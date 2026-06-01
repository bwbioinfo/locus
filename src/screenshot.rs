use std::{
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use ratatui::{Terminal, backend::TestBackend, buffer::Buffer};

use crate::{app::App, ui};

const SCREENSHOT_DIR: &str = "screenshots";

pub fn save(app: &App) -> io::Result<PathBuf> {
    let cols = app.terminal_cols.max(1);
    let rows = app.terminal_rows.max(1);
    let backend = TestBackend::new(cols, rows);
    let mut terminal = Terminal::new(backend)?;
    terminal.draw(|frame| ui::draw(frame, app))?;

    let path = screenshot_path()?;
    write_buffer(terminal.backend().buffer(), &path)?;
    Ok(path)
}

fn screenshot_path() -> io::Result<PathBuf> {
    fs::create_dir_all(SCREENSHOT_DIR)?;

    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    Ok(Path::new(SCREENSHOT_DIR).join(format!("locus-{secs}.txt")))
}

fn write_buffer(buffer: &Buffer, path: &Path) -> io::Result<()> {
    let mut file = File::create(path)?;
    for row in 0..buffer.area.height {
        let mut line = String::new();
        for col in 0..buffer.area.width {
            if let Some(cell) = buffer.cell((col, row)) {
                line.push_str(cell.symbol());
            }
        }
        writeln!(file, "{}", line.trim_end())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use ratatui::{
        buffer::Buffer,
        layout::Rect,
        style::{Color, Style},
    };

    use super::*;

    #[test]
    fn writes_buffer_as_trimmed_lines() {
        let path =
            std::env::temp_dir().join(format!("locus-screenshot-test-{}.txt", std::process::id()));
        let mut buffer = Buffer::empty(Rect::new(0, 0, 4, 2));
        buffer.set_string(0, 0, "ab", Style::default().fg(Color::White));
        buffer.set_string(1, 1, "cd", Style::default().fg(Color::White));

        write_buffer(&buffer, &path).unwrap();
        let written = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(path);

        assert_eq!(written, "ab\n cd\n");
    }
}
