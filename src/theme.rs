use ratatui::style::Color;

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum Theme {
    #[default]
    Dark,
    Light,
}

impl Theme {
    pub fn name(self) -> &'static str {
        match self {
            Self::Dark => "dark",
            Self::Light => "light",
        }
    }

    pub fn toggle(self) -> Self {
        match self {
            Self::Dark => Self::Light,
            Self::Light => Self::Dark,
        }
    }

    pub fn chrome_bg(self) -> Color {
        match self {
            Self::Dark => Color::DarkGray,
            Self::Light => Color::White,
        }
    }

    pub fn chrome_fg(self) -> Color {
        match self {
            Self::Dark => Color::White,
            Self::Light => Color::Black,
        }
    }

    pub fn brand_fg(self) -> Color {
        match self {
            Self::Dark => Color::Cyan,
            Self::Light => Color::Rgb(0, 76, 140),
        }
    }

    pub fn status_fg(self) -> Color {
        match self {
            Self::Dark => Color::Yellow,
            Self::Light => Color::Rgb(128, 82, 0),
        }
    }

    pub fn subtle_fg(self) -> Color {
        match self {
            Self::Dark => Color::DarkGray,
            Self::Light => Color::Rgb(105, 105, 105),
        }
    }

    pub fn low_contrast_fg(self) -> Color {
        match self {
            Self::Dark => Color::Gray,
            Self::Light => Color::Rgb(82, 82, 82),
        }
    }

    pub fn base_color(self, base: u8) -> Color {
        match (self, base.to_ascii_uppercase()) {
            (_, b'A') if self == Self::Dark => Color::Green,
            (_, b'T') if self == Self::Dark => Color::Red,
            (_, b'G') if self == Self::Dark => Color::Yellow,
            (_, b'C') if self == Self::Dark => Color::Blue,
            (Self::Light, b'A') => Color::Rgb(0, 110, 40),
            (Self::Light, b'T') => Color::Rgb(170, 0, 0),
            (Self::Light, b'G') => Color::Rgb(132, 89, 0),
            (Self::Light, b'C') => Color::Rgb(0, 76, 170),
            (Self::Dark, _) => Color::DarkGray,
            (Self::Light, _) => Color::Rgb(88, 88, 88),
        }
    }

    pub fn mismatch_fg(self) -> Color {
        match self {
            Self::Dark => Color::Black,
            Self::Light => Color::White,
        }
    }

    pub fn insertion_marker_fg(self) -> Color {
        match self {
            Self::Dark => Color::Black,
            Self::Light => Color::White,
        }
    }

    pub fn insertion_marker_bg(self) -> Color {
        match self {
            Self::Dark => Color::Magenta,
            Self::Light => Color::Rgb(132, 0, 132),
        }
    }

    pub fn methylation_high_fg(self) -> Color {
        match self {
            Self::Dark => Color::Black,
            Self::Light => Color::White,
        }
    }

    pub fn methylation_high_bg(self) -> Color {
        match self {
            Self::Dark => Color::Cyan,
            Self::Light => Color::Rgb(0, 102, 150),
        }
    }

    pub fn methylation_mid_bg(self) -> Color {
        match self {
            Self::Dark => Color::Cyan,
            Self::Light => Color::Rgb(170, 225, 242),
        }
    }

    pub fn methylation_low_fg(self) -> Color {
        match self {
            Self::Dark => Color::White,
            Self::Light => Color::Black,
        }
    }

    pub fn methylation_low_bg(self) -> Color {
        match self {
            Self::Dark => Color::DarkGray,
            Self::Light => Color::Rgb(224, 224, 224),
        }
    }

    pub fn feature_color(self, ty: &str) -> Color {
        match (self, ty) {
            (_, "gene" | "pseudogene") if self == Self::Dark => Color::Green,
            (_, "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA")
                if self == Self::Dark =>
            {
                Color::Yellow
            }
            (_, "exon") if self == Self::Dark => Color::Cyan,
            (_, "CDS" | "start_codon" | "stop_codon") if self == Self::Dark => Color::Blue,
            (_, "UTR" | "five_prime_UTR" | "three_prime_UTR") if self == Self::Dark => {
                Color::Magenta
            }
            (Self::Light, "gene" | "pseudogene") => Color::Rgb(0, 110, 40),
            (Self::Light, "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA") => {
                Color::Rgb(132, 89, 0)
            }
            (Self::Light, "exon") => Color::Rgb(0, 116, 116),
            (Self::Light, "CDS" | "start_codon" | "stop_codon") => Color::Rgb(0, 76, 170),
            (Self::Light, "UTR" | "five_prime_UTR" | "three_prime_UTR") => Color::Rgb(132, 0, 132),
            (Self::Dark, _) => Color::DarkGray,
            (Self::Light, _) => Color::Rgb(88, 88, 88),
        }
    }

    pub fn feature_label_fg(self) -> Color {
        match self {
            Self::Dark => Color::Black,
            Self::Light => Color::White,
        }
    }

    pub fn coverage_color(self, frac: f64) -> Color {
        match self {
            Self::Dark if frac > 0.8 => Color::Red,
            Self::Dark if frac > 0.5 => Color::Yellow,
            Self::Dark => Color::Cyan,
            Self::Light if frac > 0.8 => Color::Rgb(170, 0, 0),
            Self::Light if frac > 0.5 => Color::Rgb(132, 89, 0),
            Self::Light => Color::Rgb(0, 116, 116),
        }
    }

    pub fn html_background(self) -> &'static str {
        match self {
            Self::Dark => "#111",
            Self::Light => "#ffffff",
        }
    }

    pub fn html_foreground(self) -> &'static str {
        match self {
            Self::Dark => "#ddd",
            Self::Light => "#111111",
        }
    }
}
