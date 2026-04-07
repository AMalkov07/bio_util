# fasta-viewer

Interactive DNA sequence viewer with a SnapGene-like double-stranded display. Runs in your browser or VS Code.

![screenshot](https://img.shields.io/badge/viewer-HTML%2FJS-blue)

## Features

- Forward + complement strand displayed side-by-side
- **Live search** — type a DNA motif, see all matches highlighted instantly (both strands)
- **Live BLAST** — paste a sequence, BLAST it against the reference, navigate hits with prev/next
- Virtual scrolling — handles multi-megabase chromosomes smoothly
- Sequence selector for multi-record FASTA files
- Go-to-position navigation
- VS Code right-click integration for `.fasta` / `.fa` / `.fna` files

## Requirements

- Python 3.8+
- BLAST+ (`blastn`) — for the BLAST feature

## Install

```bash
git clone <this-repo> && cd fasta-viewer
bash install_visualizer.sh
```

This installs:
| What | Where |
|------|-------|
| Viewer script | `~/.local/share/fasta-viewer/visualize_sequence.py` |
| `visualize` command | `~/.local/bin/visualize` |
| VS Code extension | `~/.vscode/extensions/fasta-viewer/` (or `.vscode-server` on WSL) |

If `~/.local/bin` isn't in your PATH, the installer will tell you what to add to your shell rc file.

After install, reload VS Code: `Ctrl+Shift+P` → `Developer: Reload Window`.

## Uninstall

```bash
bash uninstall_visualizer.sh
```

Removes the viewer script, `visualize` command, and VS Code extension. The git repo itself is left for you to delete manually.

## Usage

```bash
# Start the interactive viewer (opens on http://127.0.0.1:8050)
visualize my_genome.fasta

# Use a different port
visualize my_genome.fasta --port 9000

# Pre-bake BLAST hits into a static HTML file
visualize my_genome.fasta --blast query.fa --output result.html
```

Or right-click any `.fasta` file in VS Code's sidebar → **Visualize Sequence**.

Open `http://127.0.0.1:8050` in VS Code's Simple Browser (`Ctrl+Shift+P` → `Simple Browser: Show`) or any web browser.

## Static HTML mode

Generate a standalone `.html` file (no server needed, but no live BLAST):

```bash
python ~/.local/share/fasta-viewer/visualize_sequence.py genome.fasta -o genome.html
python ~/.local/share/fasta-viewer/visualize_sequence.py genome.fasta --search GAATTC -o genome.html
python ~/.local/share/fasta-viewer/visualize_sequence.py genome.fasta --blast query.fa -o genome.html
```

## Keyboard shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+F` | Focus search bar |
| `Enter` | Next search match |
| `Shift+Enter` | Previous search match |
| `Escape` | Clear search |
