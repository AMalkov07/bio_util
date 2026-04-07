#!/usr/bin/env bash
# Uninstall the FASTA sequence visualizer.
#
# What it removes:
#   1. ~/.local/share/fasta-viewer/         (viewer script)
#   2. ~/.local/bin/visualize               (CLI command)
#   3. VS Code extension (fasta-viewer)
#
# Usage:
#   bash uninstall_visualizer.sh

set -euo pipefail

INSTALL_DIR="$HOME/.local/share/fasta-viewer"
BIN_CMD="$HOME/.local/bin/visualize"

echo "FASTA Viewer Uninstaller"
echo "========================"
echo ""
echo "This will remove:"
[ -d "$INSTALL_DIR" ] && echo "  - $INSTALL_DIR"
[ -f "$BIN_CMD" ]     && echo "  - $BIN_CMD"

# Find VS Code extension
VSCODE_EXT_DIR=""
for dir in "$HOME/.vscode-server/extensions/fasta-viewer" "$HOME/.vscode/extensions/fasta-viewer"; do
    if [ -d "$dir" ] || [ -L "$dir" ]; then
        VSCODE_EXT_DIR="$dir"
        echo "  - $VSCODE_EXT_DIR"
    fi
done

echo ""
read -rp "Proceed? [y/N] " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

# ── 1. Remove viewer script ──
if [ -d "$INSTALL_DIR" ]; then
    rm -rf "$INSTALL_DIR"
    echo "Removed $INSTALL_DIR"
fi

# ── 2. Remove CLI command ──
if [ -f "$BIN_CMD" ]; then
    rm "$BIN_CMD"
    echo "Removed $BIN_CMD"
fi

# ── 3. Remove VS Code extension ──
if [ -n "$VSCODE_EXT_DIR" ]; then
    rm -rf "$VSCODE_EXT_DIR"
    echo "Removed $VSCODE_EXT_DIR"
fi

echo ""
echo "Done. Reload VS Code (Ctrl+Shift+P → 'Reload Window') to remove the right-click menu."
echo ""
echo "The git repo at $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) was NOT removed."
echo "Delete it manually if you no longer need it:"
echo "  rm -rf $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
