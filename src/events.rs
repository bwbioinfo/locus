use anyhow::Result;
use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyModifiers};

use crate::app::{App, Mode};

pub fn handle_events(app: &mut App) -> Result<bool> {
    if !event::poll(std::time::Duration::from_millis(100))? {
        return Ok(true);
    }

    match event::read()? {
        Event::Key(key) => dispatch(app, key)?,
        Event::Resize(cols, rows) => {
            app.terminal_cols = cols;
            app.terminal_rows = rows;
            app.relayout();
        }
        _ => {}
    }

    Ok(!app.should_quit)
}

fn dispatch(app: &mut App, key: KeyEvent) -> Result<()> {
    match &app.mode {
        Mode::Normal => handle_normal(app, key),
        Mode::GoTo => handle_goto(app, key),
        Mode::FeatureSearch => handle_feature_search(app, key),
        Mode::ContigSelect => handle_contig_select(app, key),
        Mode::Help => handle_help(app, key),
    }
}

fn handle_normal(app: &mut App, key: KeyEvent) -> Result<()> {
    let step = app.view_span().max(1) / 5;
    let big_step = app.view_span().max(1) / 2;

    match key.code {
        KeyCode::Char('q') => app.should_quit = true,
        KeyCode::Char('?') => app.show_help = !app.show_help,

        KeyCode::Char('h') | KeyCode::Left => app.pan(-(step as i64)),
        KeyCode::Char('l') | KeyCode::Right => app.pan(step as i64),
        KeyCode::Char('H') => app.pan(-(big_step as i64)),
        KeyCode::Char('L') => app.pan(big_step as i64),

        KeyCode::Char('+') | KeyCode::Char('=') | KeyCode::Up => app.zoom_in(),
        KeyCode::Char('-') | KeyCode::Down => app.zoom_out(),
        KeyCode::Char('i') => app.toggle_insertions(),

        KeyCode::Char('g') => {
            app.mode = Mode::GoTo;
            app.command_buffer.clear();
        }

        KeyCode::Char('f') if app.gff.is_some() => {
            app.mode = Mode::FeatureSearch;
            app.command_buffer.clear();
            app.feature_matches.clear();
        }

        // n / N: cycle through feature search results without re-opening the overlay
        KeyCode::Char('n') if !app.feature_matches.is_empty() => {
            let _ = app.next_feature_match();
        }
        KeyCode::Char('N') if !app.feature_matches.is_empty() => {
            let _ = app.prev_feature_match();
        }

        KeyCode::Char('c') => {
            if key.modifiers.contains(KeyModifiers::CONTROL) {
                app.should_quit = true;
            } else {
                app.mode = Mode::ContigSelect;
            }
        }

        KeyCode::Char('r') => {
            app.needs_fetch = true;
        }

        KeyCode::Char('s') => {
            app.save_screenshot();
        }

        _ => {}
    }
    Ok(())
}

fn handle_goto(app: &mut App, key: KeyEvent) -> Result<()> {
    match key.code {
        KeyCode::Enter => {
            let _ = app.confirm_goto();
        }
        KeyCode::Esc => app.cancel_input(),
        KeyCode::Backspace => {
            app.command_buffer.pop();
        }
        KeyCode::Char(c) => app.handle_goto_input(c),
        _ => {}
    }
    Ok(())
}

fn handle_feature_search(app: &mut App, key: KeyEvent) -> Result<()> {
    match key.code {
        KeyCode::Esc => {
            app.cancel_input();
        }
        KeyCode::Enter => {
            // Jump to the currently highlighted match and close overlay
            let _ = app.jump_to_current_match();
            app.command_buffer.clear();
            app.mode = Mode::Normal;
        }
        KeyCode::Tab => {
            // Cycle forward through results without closing
            let _ = app.next_feature_match();
        }
        KeyCode::BackTab => {
            let _ = app.prev_feature_match();
        }
        KeyCode::Down => {
            let _ = app.next_feature_match();
        }
        KeyCode::Up => {
            let _ = app.prev_feature_match();
        }
        KeyCode::Backspace => {
            app.command_buffer.pop();
            app.run_feature_search();
        }
        KeyCode::Char(c) => {
            app.command_buffer.push(c);
            app.run_feature_search();
        }
        _ => {}
    }
    Ok(())
}

fn handle_contig_select(app: &mut App, key: KeyEvent) -> Result<()> {
    match key.code {
        KeyCode::Esc | KeyCode::Char('q') => app.mode = Mode::Normal,
        KeyCode::Enter => {
            if let Ok(n) = app.command_buffer.trim().parse::<usize>() {
                app.select_contig(n.saturating_sub(1));
            }
            app.command_buffer.clear();
            app.mode = Mode::Normal;
        }
        KeyCode::Backspace => {
            app.command_buffer.pop();
        }
        KeyCode::Char(c) if c.is_ascii_digit() => app.command_buffer.push(c),
        _ => {}
    }
    Ok(())
}

fn handle_help(app: &mut App, key: KeyEvent) -> Result<()> {
    match key.code {
        KeyCode::Esc | KeyCode::Char('q') | KeyCode::Char('?') => {
            app.show_help = false;
            app.mode = Mode::Normal;
        }
        _ => {}
    }
    Ok(())
}
