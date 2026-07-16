use std::collections::{HashMap, HashSet};

use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Modifier, Style},
    widgets::Widget,
};

use crate::{gff::GffFeature, theme::Theme};

use super::ViewTransform;

use crate::cache::Strand;

fn feature_glyph(ty: &str) -> char {
    match ty {
        "exon" => '█',
        "CDS" | "start_codon" | "stop_codon" => '▓',
        "UTR" | "five_prime_UTR" | "three_prime_UTR" => '▒',
        "gene" | "pseudogene" | "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA" => '─',
        _ => '·',
    }
}

fn is_intron_backbone(ty: &str) -> bool {
    matches!(
        ty,
        "gene" | "pseudogene" | "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA"
    )
}

/// Priority order for stacking rows: genes at top, then transcripts, then sub-features.
fn feature_priority(ty: &str) -> u8 {
    match ty {
        "gene" | "pseudogene" => 0,
        "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA" => 1,
        "exon" => 2,
        "CDS" => 3,
        _ => 4,
    }
}

pub struct FeaturesTrack<'a> {
    pub features: &'a [&'a GffFeature],
    pub transform: ViewTransform,
    pub theme: Theme,
}

#[derive(Debug)]
struct FeatureBlock<'a> {
    feature_type: &'a str,
    start: u64,
    end: u64,
}

#[derive(Debug)]
struct FeatureModel<'a> {
    key: String,
    name: &'a str,
    strand: Option<Strand>,
    start: u64,
    end: u64,
    blocks: Vec<FeatureBlock<'a>>,
}

impl<'a> Widget for FeaturesTrack<'a> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height == 0 || area.width == 0 || self.features.is_empty() {
            return;
        }

        let models = build_feature_models(self.features);
        let grouped_keys: HashSet<&str> = models.iter().map(|model| model.key.as_str()).collect();

        // Lay features into rows (greedy packing, same algorithm as reads).
        let mut rows: Vec<Vec<FeatureRenderItem<'_>>> = Vec::new();
        let mut row_ends: Vec<u16> = Vec::new();

        let mut items = Vec::new();
        for feature in self.features {
            if should_skip_grouped_feature(feature, &grouped_keys) {
                continue;
            }
            items.push(FeatureRenderItem::Feature(feature));
        }
        for model in &models {
            items.push(FeatureRenderItem::Model(model));
        }
        items.sort_by(|a, b| {
            a.priority()
                .cmp(&b.priority())
                .then_with(|| a.start().cmp(&b.start()))
                .then_with(|| a.end().cmp(&b.end()))
                .then_with(|| a.sort_name().cmp(b.sort_name()))
        });

        for item in items {
            let (col_start, col_end) = self.transform.bp_range_to_cols(item.start(), item.end());
            let target = row_ends
                .iter()
                .position(|&end| col_start > end)
                .unwrap_or(row_ends.len());

            if target >= area.height as usize {
                continue;
            }
            if target == rows.len() {
                rows.push(Vec::new());
                row_ends.push(0);
            }
            rows[target].push(item);
            row_ends[target] = col_end;
        }

        for (row_idx, row) in rows.iter().enumerate() {
            let y = area.y + row_idx as u16;
            for item in row {
                match item {
                    FeatureRenderItem::Feature(feature) => {
                        render_feature(feature, y, area, &self.transform, self.theme, buf);
                    }
                    FeatureRenderItem::Model(model) => {
                        render_model(model, y, area, &self.transform, self.theme, buf);
                    }
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
enum FeatureRenderItem<'a> {
    Feature(&'a GffFeature),
    Model(&'a FeatureModel<'a>),
}

impl FeatureRenderItem<'_> {
    fn priority(&self) -> u8 {
        match self {
            FeatureRenderItem::Feature(feature) => feature_priority(&feature.feature_type),
            FeatureRenderItem::Model(_) => 1,
        }
    }

    fn start(&self) -> u64 {
        match self {
            FeatureRenderItem::Feature(feature) => feature.start,
            FeatureRenderItem::Model(model) => model.start,
        }
    }

    fn end(&self) -> u64 {
        match self {
            FeatureRenderItem::Feature(feature) => feature.end,
            FeatureRenderItem::Model(model) => model.end,
        }
    }

    fn sort_name(&self) -> &str {
        match self {
            FeatureRenderItem::Feature(feature) => feature.display_name(),
            FeatureRenderItem::Model(model) => &model.key,
        }
    }
}

fn build_feature_models<'a>(features: &'a [&'a GffFeature]) -> Vec<FeatureModel<'a>> {
    let mut id_counts: HashMap<&str, usize> = HashMap::new();
    for feature in features {
        if is_block_feature(&feature.feature_type)
            && let Some(id) = feature.id.as_deref()
        {
            *id_counts.entry(id).or_default() += 1;
        }
    }

    let mut models: HashMap<String, FeatureModel<'a>> = HashMap::new();
    for feature in features {
        if !is_block_feature(&feature.feature_type) {
            continue;
        }
        let Some(key) = grouped_block_key(feature, &id_counts) else {
            continue;
        };
        let model = models.entry(key.clone()).or_insert_with(|| FeatureModel {
            key,
            name: feature.display_name(),
            strand: feature.strand,
            start: feature.start,
            end: feature.end,
            blocks: Vec::new(),
        });
        model.start = model.start.min(feature.start);
        model.end = model.end.max(feature.end);
        if model.strand.is_none() {
            model.strand = feature.strand;
        }
        model.blocks.push(FeatureBlock {
            feature_type: &feature.feature_type,
            start: feature.start,
            end: feature.end,
        });
    }

    let mut models: Vec<_> = models
        .into_values()
        .filter(|model| model.blocks.len() > 1)
        .collect();
    for model in &mut models {
        model.blocks.sort_by_key(|block| {
            (
                block.start,
                block.end,
                block_draw_priority(block.feature_type),
            )
        });
    }
    models.sort_by_key(|model| (model.start, model.end, model.key.clone()));
    models
}

fn grouped_block_key(feature: &GffFeature, id_counts: &HashMap<&str, usize>) -> Option<String> {
    let id = feature.id.as_deref();
    if id.and_then(|id| id_counts.get(id)).copied().unwrap_or(0) > 1 {
        id.map(str::to_string)
    } else {
        feature
            .parent
            .as_deref()
            .or(id)
            .or(feature.gene_name.as_deref())
            .map(str::to_string)
    }
}

fn should_skip_grouped_feature(feature: &GffFeature, grouped_keys: &HashSet<&str>) -> bool {
    if is_block_feature(&feature.feature_type) {
        return feature
            .id
            .as_deref()
            .into_iter()
            .chain(feature.parent.as_deref())
            .chain(feature.gene_name.as_deref())
            .any(|key| grouped_keys.contains(key));
    }

    if !is_intron_backbone(&feature.feature_type) {
        return false;
    }

    feature
        .id
        .as_deref()
        .into_iter()
        .chain(feature.parent.as_deref())
        .any(|key| grouped_keys.contains(key))
}

fn block_draw_priority(ty: &str) -> u8 {
    match ty {
        "exon" => 0,
        "UTR" | "five_prime_UTR" | "three_prime_UTR" => 1,
        "CDS" | "start_codon" | "stop_codon" => 2,
        _ => 3,
    }
}

fn is_block_feature(ty: &str) -> bool {
    matches!(
        ty,
        "exon"
            | "CDS"
            | "start_codon"
            | "stop_codon"
            | "UTR"
            | "five_prime_UTR"
            | "three_prime_UTR"
    )
}

fn render_model(
    model: &FeatureModel<'_>,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    theme: Theme,
    buf: &mut Buffer,
) {
    let pseudo_feature = GffFeature {
        seqname: String::new(),
        feature_type: "transcript".to_string(),
        start: model.start,
        end: model.end,
        strand: model.strand,
        id: Some(model.key.clone()),
        name: Some(model.name.to_string()),
        parent: None,
        gene_name: None,
    };

    render_feature(&pseudo_feature, y, area, transform, theme, buf);

    for block in &model.blocks {
        render_block(block, y, area, transform, theme, buf);
    }
}

fn render_block(
    block: &FeatureBlock<'_>,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    theme: Theme,
    buf: &mut Buffer,
) {
    let (col_start, col_end) = transform.bp_range_to_cols(block.start, block.end);
    let x_start = area.x + col_start;
    let x_end = (area.x + col_end).min(area.x + area.width);
    if x_start >= x_end {
        return;
    }

    let style = Style::default()
        .fg(theme.feature_label_fg())
        .bg(theme.feature_color(block.feature_type))
        .add_modifier(Modifier::BOLD);
    let ch = feature_glyph(block.feature_type);
    for x in x_start..x_end {
        if let Some(cell) = buf.cell_mut((x, y)) {
            cell.set_char(ch).set_style(style);
        }
    }
}

fn render_feature(
    feat: &GffFeature,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    theme: Theme,
    buf: &mut Buffer,
) {
    let (col_start, col_end) = transform.bp_range_to_cols(feat.start, feat.end);
    let x_start = area.x + col_start;
    let x_end = (area.x + col_end).min(area.x + area.width);
    if x_start >= x_end {
        return;
    }

    let color = theme.feature_color(&feat.feature_type);
    let style = Style::default().fg(color);
    let width = (x_end - x_start) as usize;

    let body_ch = feature_glyph(&feat.feature_type);

    // Fill body. Transcript/gene spans act as intron/backbone lines; directional
    // ticks make strand visible without making them look like exon blocks.
    for x in x_start..x_end {
        let ch = if is_intron_backbone(&feat.feature_type) {
            intron_backbone_char(feat, x - x_start)
        } else {
            body_ch
        };
        if let Some(cell) = buf.cell_mut((x, y)) {
            cell.set_char(ch).set_style(style);
        }
    }

    // Overlay name label only on backbone spans. Exon/CDS blocks are usually
    // repeated and short, so labels obscure the shape more than they help.
    if !is_intron_backbone(&feat.feature_type) {
        return;
    }

    // Overlay name label if there's enough room (>= 3 cols)
    let name = feat.display_name();
    if width >= 3 && !name.is_empty() {
        let label: String = name.chars().take(width.saturating_sub(2)).collect();
        let label_x = x_start + 1;
        let label_style = Style::default()
            .fg(theme.feature_label_fg())
            .bg(color)
            .add_modifier(Modifier::BOLD);
        for (i, ch) in label.chars().enumerate() {
            let lx = label_x + i as u16;
            if lx < x_end
                && let Some(cell) = buf.cell_mut((lx, y))
            {
                cell.set_char(ch).set_style(label_style);
            }
        }
    }
}

fn intron_backbone_char(feat: &GffFeature, offset: u16) -> char {
    if offset % 6 != 3 {
        return '─';
    }

    match feat.strand {
        Some(Strand::Forward) => '>',
        Some(Strand::Reverse) => '<',
        None => '─',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exons_and_backbones_use_distinct_glyphs() {
        assert_eq!(feature_glyph("exon"), '█');
        assert_eq!(feature_glyph("CDS"), '▓');
        assert_eq!(feature_glyph("transcript"), '─');
        assert!(is_block_feature("exon"));
        assert!(is_intron_backbone("transcript"));
        assert!(!is_intron_backbone("exon"));
    }

    #[test]
    fn repeated_gtf_transcript_ids_form_one_model() {
        let exon_a = feature("exon", 100, 150, Some("TX1"), Some("GENE1"), Some("TPTE2"));
        let exon_b = feature("exon", 300, 350, Some("TX1"), Some("GENE1"), Some("TPTE2"));
        let features = vec![&exon_a, &exon_b];

        let models = build_feature_models(&features);

        assert_eq!(models.len(), 1);
        assert_eq!(models[0].key, "TX1");
        assert_eq!(models[0].start, 100);
        assert_eq!(models[0].end, 350);
        assert_eq!(models[0].blocks.len(), 2);
    }

    #[test]
    fn grouped_transcript_record_is_not_rendered_twice() {
        let transcript = feature(
            "transcript",
            100,
            350,
            Some("TX1"),
            Some("GENE1"),
            Some("TPTE2"),
        );
        let exon_a = feature("exon", 100, 150, Some("TX1"), Some("GENE1"), Some("TPTE2"));
        let exon_b = feature("exon", 300, 350, Some("TX1"), Some("GENE1"), Some("TPTE2"));
        let features = vec![&transcript, &exon_a, &exon_b];
        let models = build_feature_models(&features);
        let grouped_keys: HashSet<&str> = models.iter().map(|model| model.key.as_str()).collect();

        assert!(should_skip_grouped_feature(&transcript, &grouped_keys));
        assert!(should_skip_grouped_feature(&exon_a, &grouped_keys));
    }

    #[test]
    fn feature_models_have_deterministic_tie_breakers() {
        let tx2_a = feature("exon", 100, 150, Some("TX2"), Some("GENE1"), Some("TPTE2"));
        let tx2_b = feature("exon", 300, 350, Some("TX2"), Some("GENE1"), Some("TPTE2"));
        let tx1_a = feature("exon", 100, 150, Some("TX1"), Some("GENE1"), Some("TPTE2"));
        let tx1_b = feature("exon", 300, 350, Some("TX1"), Some("GENE1"), Some("TPTE2"));
        let features = vec![&tx2_a, &tx2_b, &tx1_a, &tx1_b];

        let models = build_feature_models(&features);

        assert_eq!(models.len(), 2);
        assert_eq!(models[0].key, "TX1");
        assert_eq!(models[1].key, "TX2");
    }

    #[test]
    fn gene_name_grouped_blocks_do_not_hide_gene_rows() {
        let gene = feature("gene", 100, 350, Some("GENE1"), None, Some("TPTE2"));
        let exon_a = feature("exon", 100, 150, None, None, Some("TPTE2"));
        let exon_b = feature("exon", 300, 350, None, None, Some("TPTE2"));
        let features = vec![&gene, &exon_a, &exon_b];
        let models = build_feature_models(&features);
        let grouped_keys: HashSet<&str> = models.iter().map(|model| model.key.as_str()).collect();

        assert_eq!(models.len(), 1);
        assert_eq!(models[0].key, "TPTE2");
        assert!(!should_skip_grouped_feature(&gene, &grouped_keys));
        assert!(should_skip_grouped_feature(&exon_a, &grouped_keys));
    }

    fn feature(
        feature_type: &str,
        start: u64,
        end: u64,
        id: Option<&str>,
        parent: Option<&str>,
        name: Option<&str>,
    ) -> GffFeature {
        GffFeature {
            seqname: "chr1".to_string(),
            feature_type: feature_type.to_string(),
            start,
            end,
            strand: Some(Strand::Forward),
            id: id.map(str::to_string),
            name: name.map(str::to_string),
            parent: parent.map(str::to_string),
            gene_name: name.map(str::to_string),
        }
    }
}
