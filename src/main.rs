mod app;
mod bam;
mod cache;
mod cli;
mod error;
mod events;
mod gff;
mod methylation;
mod reference;
mod region;
mod render;
mod screenshot;
mod theme;
mod ui;

use anyhow::{Context, Result};
use clap::Parser;
use crossterm::{
    event::{DisableMouseCapture, EnableMouseCapture},
    execute,
    terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode},
};
use ratatui::{Terminal, backend::CrosstermBackend};

use app::App;
use bam::BamSource;
use cli::{Args, Command};
use error::LocusError;
use gff::{GffStore, prepare_indexed_annotation};
use reference::ReferenceStore;
use region::parse_region;
use theme::Theme;

fn main() -> Result<()> {
    let args = Args::parse();

    if let Some(command) = args.command {
        return match command {
            Command::PrepareAnnotations { input, output } => {
                let prepared = prepare_indexed_annotation(&input, &output)?;
                println!(
                    "prepared {} records: {} + {}",
                    prepared.record_count,
                    prepared.output_path.display(),
                    prepared.index_path.display()
                );
                Ok(())
            }
        };
    }

    let Some(bam) = args.bam.as_ref() else {
        anyhow::bail!("missing BAM path");
    };

    let source = BamSource::open(bam).with_context(|| format!("opening {bam}"))?;

    let gff = if let Some(ref path) = args.gff {
        Some(GffStore::load(path).with_context(|| format!("loading annotation {path}"))?)
    } else {
        None
    };

    let initial_region = if let Some(ref region) = args.region {
        Some(resolve_initial_region(region, &source, gff.as_ref())?)
    } else {
        source.first_mapped_region()?
    };

    let reference = if let Some(ref path) = args.reference {
        Some(ReferenceStore::load(path).with_context(|| format!("loading reference {path}"))?)
    } else {
        None
    };

    let theme = if args.light {
        Theme::Light
    } else {
        Theme::Dark
    };
    let mut app = App::new(source, gff, reference, initial_region, theme, args.min_mapq)?;

    enable_raw_mode()?;
    let mut stdout = std::io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let size = terminal.size()?;
    app.terminal_cols = size.width;
    app.terminal_rows = size.height;

    if let Err(e) = app.refresh() {
        app.status_msg = Some(format!("{e}"));
    }

    let result = run(&mut terminal, &mut app);

    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    result
}

fn resolve_initial_region(
    query: &str,
    source: &BamSource,
    gff: Option<&GffStore>,
) -> Result<region::Region> {
    let parsed = parse_region(query)?;
    if query.contains(':') || source.contig_len(&parsed.contig).is_some() {
        return source.resolve_region(&parsed);
    }

    let Some(gff) = gff else {
        return source.resolve_region(&parsed);
    };
    let feature = gff
        .resolve_feature(query)
        .ok_or_else(|| LocusError::UnknownFeature(query.to_string()))?;
    source.resolve_region(&feature.padded_region())
}

fn run<B>(terminal: &mut Terminal<B>, app: &mut App) -> Result<()>
where
    B: ratatui::backend::Backend,
    B::Error: std::error::Error + Send + Sync + 'static,
{
    loop {
        if app.needs_fetch
            && let Err(e) = app.refresh()
        {
            app.status_msg = Some(format!("{e}"));
            app.needs_fetch = false;
        }

        terminal.draw(|frame| ui::draw(frame, app))?;

        if !events::handle_events(app)? {
            break;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn demo_source() -> BamSource {
        BamSource::open(Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam"))
            .expect("open demo BAM")
    }

    fn demo_gff() -> GffStore {
        GffStore::load(Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.gff"))
            .expect("open demo GFF")
    }

    #[test]
    fn startup_region_resolves_a_feature_name_from_annotations() {
        let source = demo_source();
        let gff = demo_gff();

        let region = resolve_initial_region("demo1", &source, Some(&gff)).expect("resolve feature");

        assert_eq!(region, region::Region::new("chrDemo", 36, 128));
    }

    #[test]
    fn startup_region_keeps_coordinate_and_contig_inputs() {
        let source = demo_source();
        let gff = demo_gff();

        assert_eq!(
            resolve_initial_region("chrDemo:50-60", &source, Some(&gff)).expect("coordinate"),
            region::Region::new("chrDemo", 49, 60)
        );
        assert_eq!(
            resolve_initial_region("chrDemo", &source, Some(&gff)).expect("contig"),
            region::Region::new("chrDemo", 0, 154)
        );
    }

    #[test]
    fn startup_region_reports_an_unknown_annotation_feature() {
        let source = demo_source();
        let gff = demo_gff();

        let error = resolve_initial_region("NO_SUCH_FEATURE", &source, Some(&gff))
            .expect_err("unknown feature should fail");

        assert!(
            error
                .to_string()
                .contains("Unknown feature: NO_SUCH_FEATURE")
        );
    }
}
