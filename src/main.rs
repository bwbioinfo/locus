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
use gff::{GffStore, prepare_indexed_annotation};
use reference::ReferenceStore;
use region::parse_region;

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

    let initial_region = if let Some(ref r) = args.region {
        let parsed = parse_region(r)?;
        let resolved = source.resolve_region(&parsed)?;
        Some(resolved)
    } else {
        None
    };

    let gff = if let Some(ref path) = args.gff {
        Some(GffStore::load(path).with_context(|| format!("loading annotation {path}"))?)
    } else {
        None
    };

    let reference = if let Some(ref path) = args.reference {
        Some(ReferenceStore::load(path).with_context(|| format!("loading reference {path}"))?)
    } else {
        None
    };

    let mut app = App::new(source, gff, reference, initial_region)?;

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
