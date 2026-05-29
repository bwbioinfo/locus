mod app;
mod bam;
mod cache;
mod cli;
mod error;
mod events;
mod gff;
mod region;
mod render;
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
use cli::Args;
use gff::GffStore;
use region::parse_region;

fn main() -> Result<()> {
    let args = Args::parse();

    let source = BamSource::open(&args.bam).with_context(|| format!("opening {}", args.bam))?;

    let initial_region = if let Some(ref r) = args.region {
        let parsed = parse_region(r)?;
        let resolved = source.resolve_region(&parsed)?;
        Some(resolved)
    } else {
        None
    };

    let gff = if let Some(ref path) = args.gff {
        Some(GffStore::load(path).with_context(|| format!("loading GFF {path}"))?)
    } else {
        None
    };

    let mut app = App::new(source, gff, initial_region)?;

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

fn run<B: ratatui::backend::Backend>(terminal: &mut Terminal<B>, app: &mut App) -> Result<()> {
    loop {
        if app.needs_fetch {
            if let Err(e) = app.refresh() {
                app.status_msg = Some(format!("{e}"));
                app.needs_fetch = false;
            }
        }

        terminal.draw(|frame| ui::draw(frame, app))?;

        if !events::handle_events(app)? {
            break;
        }
    }
    Ok(())
}
