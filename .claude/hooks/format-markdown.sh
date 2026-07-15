#!/bin/bash
# PostToolUse hook: auto-format Markdown files edited by Claude Code with
# prettier (general formatting incl. table-column alignment) then
# markdownlint-cli2 --fix (project style rules, e.g. MD049 asterisk emphasis).
# Order matters: prettier defaults to underscore emphasis, which markdownlint's
# fix pass then corrects back to asterisk per .markdownlint.jsonc.
# Receives hook JSON on stdin; never blocks (always exits 0).

file_path=$(python3 -c "import json,sys; print(json.load(sys.stdin).get('tool_input',{}).get('file_path',''))" 2>/dev/null)

[[ "$file_path" == *.md && -f "$file_path" ]] || exit 0

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"

run_tool() {
    local name="$1"
    shift
    if [[ -x "$PROJECT_DIR/node_modules/.bin/$name" ]]; then
        "$PROJECT_DIR/node_modules/.bin/$name" "$@"
    elif command -v "$name" >/dev/null 2>&1; then
        "$name" "$@"
    elif command -v npx >/dev/null 2>&1; then
        npx --yes "$name" "$@"
    fi
}

run_tool prettier --write "$file_path" >/dev/null 2>&1
run_tool markdownlint-cli2 --fix "$file_path" >/dev/null 2>&1

exit 0
