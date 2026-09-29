#!/bin/sh
# Build CopperRodValidation.html -- a single self-contained file that opens in
# any browser with every figure inlined, nothing beside it.
#
#     sh build_doc.sh
#
# WHY HTML AND NOT MARKDOWN. An earlier version inlined the figures into
# CopperRodValidation.md as base64 data URIs. The file was valid but unusable:
# one table row of two mesh figures made a single line 1.8 MILLION characters
# long, and editors respond by dropping rendering and showing raw base64. A
# browser handles the same content without complaint, and shows the pictures on
# a double click rather than only in a preview pane.
#
# CopperRodValidation.md keeps relative paths: small, readable, correct inside
# the repository. This script produces the portable copy.
#
# Re-run after regenerating figures with make_plots.py, or the document will
# still show the previous run's pictures.
set -e
cd "$(dirname "$0")"
SRC=page.html
OUT=CopperRodValidation.html
[ -f "$SRC" ] || { echo "error: $SRC missing (the styled page body)"; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

{
  echo '<!doctype html><html lang="en"><head><meta charset="utf-8">'
  echo '<meta name="viewport" content="width=device-width,initial-scale=1">'
  cat "$SRC"
  echo '</body></html>'
} > "$TMP/raw.html"

# fig/<alias>.png -> the real plot file
cat > "$TMP/map.txt" <<'MAP'
mesh_domain 01_mesh_domain.png
mesh_wire 02_mesh_wire.png
phi_real 03_phi_on_surface.png
phi_mag 03b_phi_magnitude_surface.png
b_domain 04_B_magnitude.png
b_zoom 04b_B_magnitude_zoom.png
e_domain_nodal 05c_E_magnitude_nodal.png
e_domain_nodal_zoom 05d_E_magnitude_nodal_zoom.png
e_wire 06_E_in_wire.png
e_wire_true 06b_E_in_wire_true_scale.png
j_vectors 07_J_vectors.png
b_domain_smoothed 08_B_magnitude_SMOOTHED.png
b_zoom_smoothed 08b_B_magnitude_zoom_SMOOTHED.png
MAP
while read alias file; do
  [ -n "$alias" ] || continue
  base64 -w0 "output/plots/$file" > "$TMP/$alias.txt"
done < "$TMP/map.txt"

awk -v d="$TMP" '
{
  line = $0
  while (match(line, /"fig\/[A-Za-z0-9_]+\.png"/)) {
    ref   = substr(line, RSTART, RLENGTH)
    alias = substr(ref, 6, RLENGTH - 10)
    bf = d "/" alias ".txt"; data = ""
    if ((getline data < bf) > 0) {
      close(bf)
      line = substr(line,1,RSTART-1) "\"data:image/png;base64," data "\"" substr(line, RSTART+RLENGTH)
    } else { close(bf); print "MISSING: " alias > "/dev/stderr"; break }
  }
  print line
}' "$TMP/raw.html" > "$OUT"

n=$(grep -o 'data:image/png;base64,' "$OUT" | wc -l)
r=$(grep -c 'src="fig/' "$OUT" || true)
echo "wrote $OUT: $n figures inlined, $r unresolved, $(wc -c < "$OUT") bytes"

# --- PDF -------------------------------------------------------------------
# Printed from the self-contained HTML with headless Edge, so the PDF picks up
# the inlined figures and needs nothing beside it either. Edge writes noisy
# task-manager and sync errors to stderr even on success; the exit code and the
# output file are what matter.
EDGE="/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
if [ -x "$EDGE" ]; then
  ABS=$(pwd -W 2>/dev/null || pwd)
  rm -f CopperRodValidation.pdf
  "$EDGE" --headless=new --disable-gpu --no-first-run --no-pdf-header-footer \
          --user-data-dir="${TMPDIR:-/tmp}/edgepdf" --virtual-time-budget=30000 \
          --print-to-pdf="$ABS/CopperRodValidation.pdf" \
          "file:///$ABS/CopperRodValidation.html" >/dev/null 2>&1 || true
  if [ -f CopperRodValidation.pdf ]; then
    echo "wrote CopperRodValidation.pdf: $(wc -c < CopperRodValidation.pdf) bytes"
  else
    echo "PDF step failed; the HTML above is still good"
  fi
else
  echo "Edge not found, skipping the PDF"
fi
#
# NOTE: the PDF step launches a browser, which some sandboxed shells block --
# the HTML is still written and the script reports the skip. Run it from an
# ordinary terminal, or print the HTML from a browser yourself, if that happens.
