#!/usr/bin/env bash
# Print "m deg X_m" for m = 0..M for the system x1'=x2^p, x2'=x1^q.
# Usage: code/run_degrees.sh p q M   (needs Singular on PATH)
set -euo pipefail
if [ "$#" -ne 3 ]; then echo "usage: $0 p q M" >&2; exit 2; fi
command -v Singular >/dev/null || { echo "Singular not found on PATH" >&2; exit 2; }
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
p="$1"; q="$2"; M="$3"
printf '< "%s"; int m; for(m=0;m<=%s;m++){ print(string(m)+" "+string(degX(%s,%s,m))); } quit;\n' \
  "$here/degree.sing" "$M" "$p" "$q" | Singular -q
