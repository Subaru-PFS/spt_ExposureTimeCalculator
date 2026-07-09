#!/bin/bash
# PostToolUse hook: auto-format Python files edited by Claude Code with ruff + black.
# Receives hook JSON on stdin; never blocks (always exits 0).

file_path=$(python3 -c "import json,sys; print(json.load(sys.stdin).get('tool_input',{}).get('file_path',''))" 2>/dev/null)

[[ "$file_path" == *.py && -f "$file_path" ]] || exit 0

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"

run_tool() {
    local name="$1"
    shift
    if [[ -x "$PROJECT_DIR/.venv/bin/$name" ]]; then
        "$PROJECT_DIR/.venv/bin/$name" "$@"
    elif command -v "$name" >/dev/null 2>&1; then
        "$name" "$@"
    fi
}

run_tool ruff check --fix --quiet "$file_path" >/dev/null 2>&1
run_tool black --quiet "$file_path" >/dev/null 2>&1

exit 0
