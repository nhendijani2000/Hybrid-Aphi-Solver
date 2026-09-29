#!/bin/sh
# Build the portable CopperRodValidation.md from CopperRodValidation.src.md by
# inlining every figure as a base64 data URI.
#
# The .src.md uses relative paths and is the file to EDIT. The generated .md is
# self-contained, so it renders standalone anywhere -- emailed, copied, opened
# from any folder -- which the relative-path version does not.
#
#     sh build_doc.sh
#
# Re-run it after regenerating figures with make_plots.py, or the document will
# still show the previous run's pictures.
set -e
cd "$(dirname "$0")"
SRC=CopperRodValidation.src.md
OUT=CopperRodValidation.md
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
for f in $(grep -o '](output/plots/[^)]*)' "$SRC" | sed 's/](//; s/)$//' | sort -u); do
  base64 -w0 "$f" > "$TMP/$(basename "$f").txt"
done
awk -v d="$TMP" '
{
  line = $0
  while (match(line, /\]\(output\/plots\/[^)]*\)/)) {
    ref  = substr(line, RSTART, RLENGTH)
    path = substr(ref, 3, RLENGTH - 3)
    name = path; sub(/.*\//, "", name)
    bf = d "/" name ".txt"; data = ""
    if ((getline data < bf) > 0) {
      close(bf)
      line = substr(line,1,RSTART-1) "](data:image/png;base64," data ")" substr(line, RSTART+RLENGTH)
    } else { close(bf); break }
  }
  print line
}' "$SRC" > "$OUT"
n=$(grep -o 'data:image/png;base64,' "$OUT" | wc -l)
r=$(grep -c 'output/plots/' "$OUT" || true)
echo "wrote $OUT: $n figures embedded, $r relative refs left, $(wc -c < "$OUT") bytes"
