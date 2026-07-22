#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"

"$script_dir/build.sh"

cd "$repo_root"
cargo build --quiet

mkdir -p docs/images docs/captures screenshots

capture_demo() {
  local label="$1"
  shift
  local before
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
    printf 'p'
    sleep 0.3
    printf 's'
    sleep 0.8
    printf 'q'
  ) | script -qfec "stty cols 110 rows 20; ./target/debug/locus examples/demo/demo.sorted.bam --region chrDemo:45-115 --reference examples/demo/demo.fa --gff examples/demo/demo.sorted.gff.gz $*" "/tmp/locus-demo-$label.typescript" >/dev/null

  latest_html="$(comm -13 <(sort "$before") <(find screenshots -maxdepth 1 -name 'locus-*.html' -print | sort) | tail -n 1)"
  rm -f "$before"

  if [[ -z "$latest_html" ]]; then
    echo "no screenshot html was created for $label" >&2
    exit 1
  fi

  latest_txt="${latest_html%.html}.txt"
  cp -f "$latest_html" "docs/captures/demo-$label-expanded-methylation.html"
  cp -f "$latest_txt" "docs/captures/demo-$label-expanded-methylation.ansi.txt"

  if command -v chromium >/dev/null 2>&1; then
    if ! chromium --headless --no-sandbox --disable-gpu --disable-crash-reporter \
      --window-size=1100,430 \
      --screenshot="docs/images/demo-$label-expanded-methylation.png" \
      "file://$repo_root/docs/captures/demo-$label-expanded-methylation.html"; then
      echo "chromium failed; HTML and ANSI screenshots were updated, PNG was not regenerated for $label" >&2
    fi
  else
    echo "chromium not found; HTML and ANSI screenshots were updated, PNG was not regenerated" >&2
  fi
}

capture_demo dark
cp -f docs/captures/demo-dark-expanded-methylation.html docs/captures/demo-expanded-methylation.html
cp -f docs/captures/demo-dark-expanded-methylation.ansi.txt docs/captures/demo-expanded-methylation.ansi.txt
if [[ -f docs/images/demo-dark-expanded-methylation.png ]]; then
  cp -f docs/images/demo-dark-expanded-methylation.png docs/images/demo-expanded-methylation.png
fi

capture_demo light --light
