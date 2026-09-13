#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
web_source="$(mktemp)"
trap 'rm -f "$web_source"' EXIT

# Convert knitr chunk headers to standard fenced R blocks for Pandoc.
sed -E 's/^```\{r[^}]*\}$/```r/' \
  "$repo_dir/tutorials/meadow_plot_tutorial.Rmd" > "$web_source"

pandoc "$web_source" \
  --from=markdown+fenced_code_attributes \
  --to=html5 \
  --standalone \
  --section-divs \
  --toc \
  --toc-depth=3 \
  --highlight-style=zenburn \
  --template="$repo_dir/website/template.html" \
  --metadata title="Meadow Plot Tutorial" \
  --output="$repo_dir/docs/index.html"

echo "Built docs/index.html"
