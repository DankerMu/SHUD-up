#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

MAP_FILE="epic/issue-map.tsv"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

require_cmd gh
require_cmd git
require_cmd find
require_cmd sort
require_cmd head

if ! gh auth status -h github.com >/dev/null 2>&1; then
  echo "ERROR: gh not authenticated (or token invalid)." >&2
  echo "Run: gh auth login -h github.com" >&2
  exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "ERROR: not inside a git repository." >&2
  exit 1
fi

mapped_url_for_file() {
  local file="$1"
  [[ -f "$MAP_FILE" ]] || return 1
  awk -F $'\t' -v f="$file" '$1 == f {print $2; found=1; exit} END {exit(found?0:1)}' "$MAP_FILE"
}

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

find epic -type f -name 'backlog-*.md' | sort >"$tmp"

touch "$MAP_FILE"

while IFS= read -r file; do
  if url_existing="$(mapped_url_for_file "$file" 2>/dev/null)"; then
    echo "SKIP: already mapped: $file -> $url_existing"
    continue
  fi
  title="$(head -n 1 "$file" | sed 's/^# *//')"
  echo "CREATE: $title"
  url="$(gh issue create --title "$title" --body-file "$file")"
  number="${url##*/}"
  printf "%s\t%s\t%s\n" "$file" "$url" "$number" >>"$MAP_FILE"
  echo "OK: $url"
done <"$tmp"

echo "Done. Issue map written to: $MAP_FILE"
