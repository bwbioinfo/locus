use crate::cache::{ModificationStrand, ModifiedBaseCall};

#[derive(Debug, Clone, PartialEq, Eq)]
struct MmGroup<'a> {
    canonical_base: u8,
    strand: ModificationStrand,
    modifications: Vec<&'a str>,
    deltas: Vec<usize>,
}

pub fn parse_modified_bases(mm: &str, ml: Option<&[u8]>, sequence: &[u8]) -> Vec<ModifiedBaseCall> {
    let mut calls = Vec::new();
    let mut probabilities = ml.unwrap_or(&[]).iter().copied();

    for group in parse_mm_groups(mm) {
        let Some(group) = group else {
            continue;
        };
        let positions = modified_read_positions(sequence, group.canonical_base, &group.deltas);
        for read_pos in positions {
            for modification in &group.modifications {
                calls.push(ModifiedBaseCall {
                    read_pos,
                    canonical_base: group.canonical_base.to_ascii_uppercase(),
                    strand: group.strand.clone(),
                    modification: (*modification).to_string(),
                    probability: probabilities.next(),
                });
            }
        }
    }

    calls
}

fn parse_mm_groups(mm: &str) -> impl Iterator<Item = Option<MmGroup<'_>>> {
    mm.split(';')
        .filter(|group| !group.trim().is_empty())
        .map(parse_mm_group)
}

fn parse_mm_group(group: &str) -> Option<MmGroup<'_>> {
    let mut chars = group.char_indices();
    let (_, canonical_base) = chars.next()?;
    if !canonical_base.is_ascii_alphabetic() {
        return None;
    }

    let (_, strand_char) = chars.next()?;
    let strand = match strand_char {
        '+' => ModificationStrand::Forward,
        '-' => ModificationStrand::Reverse,
        _ => return None,
    };

    let delta_start = group.find(',').unwrap_or(group.len());
    let mut modification_part = &group[2..delta_start];
    if modification_part.ends_with('?') || modification_part.ends_with('.') {
        modification_part = &modification_part[..modification_part.len().saturating_sub(1)];
    }
    if modification_part.is_empty() {
        return None;
    }

    let modifications = split_modifications(modification_part);
    if modifications.is_empty() {
        return None;
    }

    let deltas = if delta_start < group.len() {
        group[delta_start + 1..]
            .split(',')
            .filter(|field| !field.is_empty())
            .map(str::parse::<usize>)
            .collect::<Result<Vec<_>, _>>()
            .ok()?
    } else {
        Vec::new()
    };

    Some(MmGroup {
        canonical_base: canonical_base as u8,
        strand,
        modifications,
        deltas,
    })
}

fn split_modifications(modification_part: &str) -> Vec<&str> {
    let mut modifications = Vec::new();
    let mut start = 0;
    for (idx, ch) in modification_part.char_indices().skip(1) {
        if ch.is_ascii_uppercase() || ch == '+' || ch == '-' {
            modifications.push(&modification_part[start..idx]);
            start = idx;
        }
    }
    modifications.push(&modification_part[start..]);
    modifications
        .into_iter()
        .filter(|modification| !modification.is_empty())
        .collect()
}

fn modified_read_positions(sequence: &[u8], canonical_base: u8, deltas: &[usize]) -> Vec<usize> {
    let canonical_positions: Vec<usize> = sequence
        .iter()
        .enumerate()
        .filter_map(|(idx, &base)| base.eq_ignore_ascii_case(&canonical_base).then_some(idx))
        .collect();

    let mut positions = Vec::new();
    let mut canonical_idx: isize = -1;
    for &delta in deltas {
        canonical_idx += delta as isize + 1;
        let Ok(idx) = usize::try_from(canonical_idx) else {
            continue;
        };
        let Some(&read_pos) = canonical_positions.get(idx) else {
            continue;
        };
        positions.push(read_pos);
    }

    positions
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_cytosine_methylation_with_probabilities() {
        let calls = parse_modified_bases("C+m,0,1;", Some(&[240, 128]), b"ACCCG");

        assert_eq!(
            calls,
            vec![
                ModifiedBaseCall {
                    read_pos: 1,
                    canonical_base: b'C',
                    strand: ModificationStrand::Forward,
                    modification: "m".to_string(),
                    probability: Some(240),
                },
                ModifiedBaseCall {
                    read_pos: 3,
                    canonical_base: b'C',
                    strand: ModificationStrand::Forward,
                    modification: "m".to_string(),
                    probability: Some(128),
                },
            ]
        );
    }

    #[test]
    fn missing_ml_leaves_probabilities_empty() {
        let calls = parse_modified_bases("C+m,0;", None, b"ACG");

        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].read_pos, 1);
        assert_eq!(calls[0].probability, None);
    }

    #[test]
    fn parses_multiple_groups_and_skip_modes() {
        let calls = parse_modified_bases("C+m?,0;A+a.,1;", Some(&[200, 100]), b"ACCA");

        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].read_pos, 1);
        assert_eq!(calls[0].modification, "m");
        assert_eq!(calls[1].read_pos, 3);
        assert_eq!(calls[1].canonical_base, b'A');
        assert_eq!(calls[1].modification, "a");
    }

    #[test]
    fn parses_reverse_strand_methylation_group() {
        let calls = parse_modified_bases("C-m,1;", Some(&[180]), b"ACCCG");

        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].read_pos, 2);
        assert_eq!(calls[0].strand, ModificationStrand::Reverse);
        assert_eq!(calls[0].probability, Some(180));
    }

    #[test]
    fn empty_mm_tag_has_no_methylation_calls() {
        let calls = parse_modified_bases("", None, b"ACCCG");

        assert!(calls.is_empty());
    }

    #[test]
    fn malformed_groups_are_ignored() {
        let calls = parse_modified_bases("bad;C+m,nope;C+m,0;", Some(&[42]), b"ACG");

        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].read_pos, 1);
        assert_eq!(calls[0].probability, Some(42));
    }

    #[test]
    fn out_of_range_delta_is_ignored() {
        let calls = parse_modified_bases("C+m,10;", Some(&[42]), b"ACG");

        assert!(calls.is_empty());
    }
}
