#!/usr/bin/env bash
# Install the FASTA sequence visualizer on any machine.
#
# What it does:
#   1. Copies visualize_sequence.py to ~/.local/share/fasta-viewer/
#   2. Creates the `visualize` command in ~/.local/bin/
#   3. Installs the VS Code right-click extension
#
# Usage:
#   bash install_visualizer.sh
#
# Requirements: python3, blastn (for BLAST features)

set -euo pipefail

INSTALL_DIR="$HOME/.local/share/fasta-viewer"
BIN_DIR="$HOME/.local/bin"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── 1. Copy the viewer script ──
mkdir -p "$INSTALL_DIR"
cp "$SCRIPT_DIR/visualize_sequence.py" "$INSTALL_DIR/visualize_sequence.py"
echo "Installed viewer script to $INSTALL_DIR/visualize_sequence.py"

# ── 2. Create the `visualize` wrapper ──
mkdir -p "$BIN_DIR"
cat > "$BIN_DIR/visualize" << 'WRAPPER'
#!/usr/bin/env bash
set -euo pipefail

VIEWER="$HOME/.local/share/fasta-viewer/visualize_sequence.py"

if [ $# -lt 1 ]; then
    echo "Usage: visualize <fasta> [options]"
    echo "Starts the interactive sequence viewer with BLAST support."
    echo ""
    echo "Examples:"
    echo "  visualize references/6991.fasta"
    echo "  visualize references/6991.fasta --port 9000"
    echo "  visualize references/6991.fasta --blast query.fa"
    exit 1
fi

FASTA="$1"
shift
python "$VIEWER" "$FASTA" --serve "$@"
WRAPPER
chmod +x "$BIN_DIR/visualize"
echo "Installed 'visualize' command to $BIN_DIR/visualize"

# ── 3. Install VS Code extension ──
# Detect extensions directory (local vs WSL/remote)
if [ -d "$HOME/.vscode-server/extensions" ]; then
    VSCODE_EXT_DIR="$HOME/.vscode-server/extensions"
elif [ -d "$HOME/.vscode/extensions" ]; then
    VSCODE_EXT_DIR="$HOME/.vscode/extensions"
else
    VSCODE_EXT_DIR="$HOME/.vscode/extensions"
    mkdir -p "$VSCODE_EXT_DIR"
fi

EXT_DIR="$VSCODE_EXT_DIR/fasta-viewer"
mkdir -p "$EXT_DIR"

cat > "$EXT_DIR/package.json" << 'EXTPKG'
{
  "name": "fasta-viewer",
  "displayName": "FASTA Sequence Viewer",
  "description": "Right-click a .fasta file to open the interactive sequence viewer",
  "version": "0.0.1",
  "engines": { "vscode": "^1.85.0" },
  "activationEvents": [],
  "main": "extension.js",
  "contributes": {
    "commands": [
      {
        "command": "fastaViewer.visualize",
        "title": "Visualize Sequence"
      }
    ],
    "menus": {
      "explorer/context": [
        {
          "command": "fastaViewer.visualize",
          "when": "resourceExtname == .fasta || resourceExtname == .fa || resourceExtname == .fna",
          "group": "navigation"
        }
      ]
    }
  }
}
EXTPKG

cat > "$EXT_DIR/extension.js" << 'EXTJS'
const vscode = require("vscode");

function activate(context) {
  const cmd = vscode.commands.registerCommand(
    "fastaViewer.visualize",
    (uri) => {
      if (!uri || !uri.fsPath) {
        vscode.window.showErrorMessage("No file selected.");
        return;
      }
      const filePath = uri.fsPath;
      const terminal =
        vscode.window.terminals.find((t) => t.name === "Sequence Viewer") ||
        vscode.window.createTerminal("Sequence Viewer");
      terminal.show();
      terminal.sendText(`visualize "${filePath}"`);
    }
  );
  context.subscriptions.push(cmd);
}

function deactivate() {}
module.exports = { activate, deactivate };
EXTJS

echo "Installed VS Code extension to $EXT_DIR"

# ── 4. Check PATH ──
if ! echo "$PATH" | tr ':' '\n' | grep -qx "$BIN_DIR"; then
    echo ""
    echo "NOTE: $BIN_DIR is not in your PATH. Add this to your ~/.bashrc or ~/.zshrc:"
    echo "  export PATH=\"\$HOME/.local/bin:\$PATH\""
fi

echo ""
echo "Done! Reload VS Code (Ctrl+Shift+P → 'Reload Window') for the right-click menu."
echo "Then run: visualize <your-file.fasta>"
