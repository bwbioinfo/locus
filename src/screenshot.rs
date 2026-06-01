use std::{
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use ratatui::{
    Terminal,
    backend::TestBackend,
    buffer::{Buffer, Cell},
    style::{Color, Modifier},
};

use crate::{app::App, ui};

const SCREENSHOT_DIR: &str = "screenshots";

pub struct ScreenshotPaths {
    pub text: PathBuf,
    pub html: PathBuf,
}

pub fn save(app: &App) -> io::Result<ScreenshotPaths> {
    let cols = app.terminal_cols.max(1);
    let rows = app.terminal_rows.max(1);
    let backend = TestBackend::new(cols, rows);
    let mut terminal = Terminal::new(backend)?;
    terminal.draw(|frame| ui::draw(frame, app))?;

    let paths = screenshot_paths()?;
    let buffer = terminal.backend().buffer();
    write_buffer(buffer, &paths.text)?;
    write_html(buffer, &paths.html)?;
    Ok(paths)
}

fn screenshot_paths() -> io::Result<ScreenshotPaths> {
    fs::create_dir_all(SCREENSHOT_DIR)?;

    let millis = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_millis();
    let stem = format!("locus-{millis}-{}", std::process::id());

    Ok(ScreenshotPaths {
        text: Path::new(SCREENSHOT_DIR).join(format!("{stem}.txt")),
        html: Path::new(SCREENSHOT_DIR).join(format!("{stem}.html")),
    })
}

fn write_buffer(buffer: &Buffer, path: &Path) -> io::Result<()> {
    let mut file = File::create(path)?;
    for row in 0..buffer.area.height {
        let mut line = String::new();
        let mut current_style = CellStyle::default();
        let mut wrote_style = false;
        let end_col = last_visible_col(buffer, row);

        for col in 0..end_col {
            if let Some(cell) = buffer.cell((col, row)) {
                let style = CellStyle::from(cell);
                if style != current_style {
                    if wrote_style {
                        line.push_str("\x1b[0m");
                    }
                    push_sgr(&mut line, style);
                    current_style = style;
                    wrote_style = true;
                }
                line.push_str(cell.symbol());
            }
        }
        if wrote_style {
            line.push_str("\x1b[0m");
        }
        writeln!(file, "{line}")?;
    }
    Ok(())
}

fn write_html(buffer: &Buffer, path: &Path) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "<!doctype html>")?;
    writeln!(file, "<html lang=\"en\">")?;
    writeln!(file, "<head>")?;
    writeln!(file, "<meta charset=\"utf-8\">")?;
    writeln!(
        file,
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
    )?;
    writeln!(file, "<title>locus screenshot</title>")?;
    writeln!(
        file,
        "<style>body{{margin:0;background:#111;color:#ddd}}pre{{margin:0;padding:16px;font:14px/1.2 ui-monospace,SFMono-Regular,Consolas,monospace;white-space:pre;overflow:auto}}span{{display:inline}}</style>"
    )?;
    writeln!(file, "</head>")?;
    writeln!(file, "<body>")?;
    writeln!(file, "<pre>")?;

    for row in 0..buffer.area.height {
        let end_col = last_visible_col(buffer, row);
        let mut current_style = CellStyle::default();
        let mut span_open = false;

        for col in 0..end_col {
            if let Some(cell) = buffer.cell((col, row)) {
                let style = CellStyle::from(cell);
                if style != current_style {
                    if span_open {
                        write!(file, "</span>")?;
                    }
                    write!(file, "<span style=\"{}\">", css_style(style))?;
                    current_style = style;
                    span_open = true;
                }
                write_escaped_html(&mut file, cell.symbol())?;
            }
        }

        if span_open {
            write!(file, "</span>")?;
        }
        writeln!(file)?;
    }

    writeln!(file, "</pre>")?;
    writeln!(file, "</body>")?;
    writeln!(file, "</html>")?;
    Ok(())
}

fn last_visible_col(buffer: &Buffer, row: u16) -> u16 {
    (0..buffer.area.width)
        .rev()
        .find(|&col| {
            buffer
                .cell((col, row))
                .map(|cell| !cell.symbol().trim_end().is_empty())
                .unwrap_or(false)
        })
        .map(|col| col + 1)
        .unwrap_or(0)
}

fn css_style(style: CellStyle) -> String {
    let mut rules = Vec::new();

    if let Some(color) = css_color(style.fg) {
        rules.push(format!("color:{color}"));
    }
    if let Some(color) = css_color(style.bg) {
        rules.push(format!("background-color:{color}"));
    }
    if style.modifier.contains(Modifier::BOLD) {
        rules.push("font-weight:700".to_string());
    }
    if style.modifier.contains(Modifier::DIM) {
        rules.push("opacity:.65".to_string());
    }
    if style.modifier.contains(Modifier::ITALIC) {
        rules.push("font-style:italic".to_string());
    }
    if style.modifier.contains(Modifier::REVERSED) {
        rules.push("filter:invert(1)".to_string());
    }

    let mut text_decoration = Vec::new();
    if style.modifier.contains(Modifier::UNDERLINED) {
        text_decoration.push("underline");
    }
    if style.modifier.contains(Modifier::CROSSED_OUT) {
        text_decoration.push("line-through");
    }
    if !text_decoration.is_empty() {
        rules.push(format!("text-decoration:{}", text_decoration.join(" ")));
    }
    if style.modifier.contains(Modifier::HIDDEN) {
        rules.push("visibility:hidden".to_string());
    }

    rules.join(";")
}

fn css_color(color: Color) -> Option<String> {
    let value = match color {
        Color::Reset => return None,
        Color::Black => "#000000".to_string(),
        Color::Red => "#800000".to_string(),
        Color::Green => "#008000".to_string(),
        Color::Yellow => "#808000".to_string(),
        Color::Blue => "#000080".to_string(),
        Color::Magenta => "#800080".to_string(),
        Color::Cyan => "#008080".to_string(),
        Color::Gray => "#c0c0c0".to_string(),
        Color::DarkGray => "#808080".to_string(),
        Color::LightRed => "#ff0000".to_string(),
        Color::LightGreen => "#00ff00".to_string(),
        Color::LightYellow => "#ffff00".to_string(),
        Color::LightBlue => "#0000ff".to_string(),
        Color::LightMagenta => "#ff00ff".to_string(),
        Color::LightCyan => "#00ffff".to_string(),
        Color::White => "#ffffff".to_string(),
        Color::Indexed(index) => indexed_color(index),
        Color::Rgb(r, g, b) => format!("#{r:02x}{g:02x}{b:02x}"),
    };

    Some(value)
}

fn indexed_color(index: u8) -> String {
    const ANSI_16: [&str; 16] = [
        "#000000", "#800000", "#008000", "#808000", "#000080", "#800080", "#008080", "#c0c0c0",
        "#808080", "#ff0000", "#00ff00", "#ffff00", "#0000ff", "#ff00ff", "#00ffff", "#ffffff",
    ];

    match index {
        0..=15 => ANSI_16[index as usize].to_string(),
        16..=231 => {
            let value = index - 16;
            let r = value / 36;
            let g = (value % 36) / 6;
            let b = value % 6;
            format!(
                "#{:02x}{:02x}{:02x}",
                indexed_color_component(r),
                indexed_color_component(g),
                indexed_color_component(b)
            )
        }
        232..=255 => {
            let gray = 8 + (index - 232) * 10;
            format!("#{gray:02x}{gray:02x}{gray:02x}")
        }
    }
}

fn indexed_color_component(value: u8) -> u8 {
    if value == 0 { 0 } else { 55 + value * 40 }
}

fn write_escaped_html<W: Write>(writer: &mut W, text: &str) -> io::Result<()> {
    for ch in text.chars() {
        match ch {
            '&' => write!(writer, "&amp;")?,
            '<' => write!(writer, "&lt;")?,
            '>' => write!(writer, "&gt;")?,
            '"' => write!(writer, "&quot;")?,
            '\'' => write!(writer, "&#39;")?,
            _ => write!(writer, "{ch}")?,
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct CellStyle {
    fg: Color,
    bg: Color,
    modifier: Modifier,
}

impl From<&Cell> for CellStyle {
    fn from(cell: &Cell) -> Self {
        Self {
            fg: cell.fg,
            bg: cell.bg,
            modifier: cell.modifier,
        }
    }
}

fn push_sgr(line: &mut String, style: CellStyle) {
    let mut codes = Vec::new();

    push_color_code(&mut codes, style.fg, 30, 90, 38);
    push_color_code(&mut codes, style.bg, 40, 100, 48);
    push_modifier_codes(&mut codes, style.modifier);

    if !codes.is_empty() {
        line.push_str("\x1b[");
        line.push_str(&codes.join(";"));
        line.push('m');
    }
}

fn push_color_code(
    codes: &mut Vec<String>,
    color: Color,
    base: u8,
    bright_base: u8,
    extended_base: u8,
) {
    let code = match color {
        Color::Reset => return,
        Color::Black => base,
        Color::Red => base + 1,
        Color::Green => base + 2,
        Color::Yellow => base + 3,
        Color::Blue => base + 4,
        Color::Magenta => base + 5,
        Color::Cyan => base + 6,
        Color::Gray => base + 7,
        Color::DarkGray => bright_base,
        Color::LightRed => bright_base + 1,
        Color::LightGreen => bright_base + 2,
        Color::LightYellow => bright_base + 3,
        Color::LightBlue => bright_base + 4,
        Color::LightMagenta => bright_base + 5,
        Color::LightCyan => bright_base + 6,
        Color::White => bright_base + 7,
        Color::Indexed(index) => {
            codes.push(format!("{extended_base};5;{index}"));
            return;
        }
        Color::Rgb(r, g, b) => {
            codes.push(format!("{extended_base};2;{r};{g};{b}"));
            return;
        }
    };

    codes.push(code.to_string());
}

fn push_modifier_codes(codes: &mut Vec<String>, modifier: Modifier) {
    if modifier.contains(Modifier::BOLD) {
        codes.push("1".to_string());
    }
    if modifier.contains(Modifier::DIM) {
        codes.push("2".to_string());
    }
    if modifier.contains(Modifier::ITALIC) {
        codes.push("3".to_string());
    }
    if modifier.contains(Modifier::UNDERLINED) {
        codes.push("4".to_string());
    }
    if modifier.contains(Modifier::SLOW_BLINK) {
        codes.push("5".to_string());
    }
    if modifier.contains(Modifier::RAPID_BLINK) {
        codes.push("6".to_string());
    }
    if modifier.contains(Modifier::REVERSED) {
        codes.push("7".to_string());
    }
    if modifier.contains(Modifier::HIDDEN) {
        codes.push("8".to_string());
    }
    if modifier.contains(Modifier::CROSSED_OUT) {
        codes.push("9".to_string());
    }
}

#[cfg(test)]
mod tests {
    use ratatui::{
        buffer::Buffer,
        layout::Rect,
        style::{Color, Modifier, Style},
    };

    use super::*;

    #[test]
    fn writes_buffer_as_trimmed_lines() {
        let path =
            std::env::temp_dir().join(format!("locus-screenshot-test-{}.txt", std::process::id()));
        let mut buffer = Buffer::empty(Rect::new(0, 0, 4, 2));
        buffer.set_string(0, 0, "ab", Style::default());
        buffer.set_string(1, 1, "cd", Style::default());

        write_buffer(&buffer, &path).unwrap();
        let written = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(path);

        assert_eq!(written, "ab\n cd\n");
    }

    #[test]
    fn writes_buffer_with_ansi_color_by_default() {
        let path = std::env::temp_dir().join(format!(
            "locus-screenshot-color-test-{}.txt",
            std::process::id()
        ));
        let mut buffer = Buffer::empty(Rect::new(0, 0, 3, 1));
        buffer.set_string(
            0,
            0,
            "A",
            Style::default()
                .fg(Color::LightRed)
                .bg(Color::Blue)
                .add_modifier(Modifier::BOLD),
        );
        buffer.set_string(1, 0, "b", Style::default().fg(Color::Green));

        write_buffer(&buffer, &path).unwrap();
        let written = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(path);

        assert_eq!(written, "\u{1b}[91;44;1mA\u{1b}[0m\u{1b}[32mb\u{1b}[0m\n");
    }

    #[test]
    fn writes_html_with_color_styles() {
        let path = std::env::temp_dir().join(format!(
            "locus-screenshot-html-test-{}.html",
            std::process::id()
        ));
        let mut buffer = Buffer::empty(Rect::new(0, 0, 4, 1));
        buffer.set_string(
            0,
            0,
            "A",
            Style::default()
                .fg(Color::LightRed)
                .bg(Color::Blue)
                .add_modifier(Modifier::BOLD),
        );
        buffer.set_string(1, 0, "<&", Style::default().fg(Color::Rgb(1, 2, 3)));

        write_html(&buffer, &path).unwrap();
        let written = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(path);

        assert!(written.contains("<!doctype html>"));
        assert!(written.contains(
            "<span style=\"color:#ff0000;background-color:#000080;font-weight:700\">A</span>"
        ));
        assert!(written.contains("<span style=\"color:#010203\">&lt;&amp;</span>"));
    }

    #[test]
    fn indexed_color_maps_xterm_cube_and_grayscale() {
        assert_eq!(indexed_color(16), "#000000");
        assert_eq!(indexed_color(21), "#0000ff");
        assert_eq!(indexed_color(231), "#ffffff");
        assert_eq!(indexed_color(232), "#080808");
        assert_eq!(indexed_color(255), "#eeeeee");
    }
}
