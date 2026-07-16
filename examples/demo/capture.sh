#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"

"$script_dir/build.sh"

cd "$repo_root"
cargo build --quiet

mkdir -p docs/images docs/captures screenshots

before="$(mktemp)"
find screenshots -maxdepth 1 -name 'locus-*.html' -print > "$before"

(
  sleep 1
  printf 'i'
  sleep 0.3
  printf '\t'
  sleep 0.3
  printf 'm'
  sleep 0.3
  printf 's'
  sleep 0.8
  printf 'q'
) | script -qfec "stty cols 110 rows 30; ./target/debug/locus examples/demo/demo.sorted.bam --region chrDemo:45-115 --reference examples/demo/demo.fa --gff examples/demo/demo.sorted.gff.gz" /tmp/locus-demo.typescript >/dev/null

latest_html="$(comm -13 <(sort "$before") <(find screenshots -maxdepth 1 -name 'locus-*.html' -print | sort) | tail -n 1)"
rm -f "$before"

if [[ -z "$latest_html" ]]; then
  echo "no screenshot html was created" >&2
  exit 1
fi

latest_txt="${latest_html%.html}.txt"
cp -f "$latest_html" docs/captures/demo-expanded-methylation.html
cp -f "$latest_txt" docs/captures/demo-expanded-methylation.ansi.txt

if command -v chromium >/dev/null 2>&1; then
  chromium --headless --no-sandbox --disable-gpu --disable-crash-reporter \
    --window-size=1500,700 \
    --screenshot=docs/images/demo-expanded-methylation.png \
    "file://$repo_root/docs/captures/demo-expanded-methylation.html"
else
  echo "chromium not found; HTML and ANSI screenshots were updated, PNG was not regenerated" >&2
fi
