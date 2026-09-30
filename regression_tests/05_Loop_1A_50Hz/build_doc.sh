#!/bin/sh
# Build LoopValidation.html and .pdf -- single self-contained files that
# open anywhere with every figure inlined, nothing beside them.
#
#     sh build_doc.sh
#
# WHY HTML AND NOT MARKDOWN FOR THE PORTABLE COPY. An earlier attempt in case 02
# inlined figures into the .md as base64 data URIs. The file was valid but
# unusable: one row of two figures made a single line 1.8 MILLION characters
# long, and editors respond by dropping rendering and showing raw base64.
#
# LoopValidation.md keeps relative fig/ paths instead -- small, readable and
# correct inside the repository, where the figures are tracked beside it. This
# script produces the portable pair.
#
# fig/e_port_pair.png is COMPOSED from two of make_plots.py's outputs by
# make_port_pair.py, so regenerating figures means running BOTH:
#     pvbatch make_plots.py output  &&  pvbatch make_port_pair.py
# Re-run after regenerating figures with make_plots.py or
# make_convergence_plot.py, or the documents keep showing the previous run's
# pictures: they are baked in.
set -e
cd "$(dirname "$0")"
SRC=page.html
OUT=LoopValidation.html
PDF=LoopValidation.pdf
[ -f "$SRC" ] || { echo "error: $SRC missing (the styled page body)"; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

{
  echo '<!doctype html><html lang="en"><head><meta charset="utf-8">'
  echo '<meta name="viewport" content="width=device-width,initial-scale=1">'
  cat "$SRC"
  echo '</body></html>'
} > "$TMP/raw.html"

# Unlike case 02 there is no alias map: make_plots.py already writes fig/<name>.png
# with the names page.html references, so the figure is found by its own path.
missing=0
for f in $(grep -o 'src="fig/[A-Za-z0-9_]*\.png"' "$SRC" | sed 's/.*fig\///; s/"//'); do
  if [ -f "fig/$f" ]; then
    base64 -w0 "fig/$f" > "$TMP/${f%.png}.txt"
  else
    echo "MISSING: fig/$f" >&2
    missing=$((missing + 1))
  fi
done
[ "$missing" -eq 0 ] || { echo "error: $missing figure(s) missing -- run make_plots.py first"; exit 1; }

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
[ "$r" -eq 0 ] || { echo "error: $r figure reference(s) did not resolve"; exit 1; }

# --- PDF -------------------------------------------------------------------
# Printed from the self-contained HTML with headless Edge, so the PDF picks up
# the inlined figures and needs nothing beside it either.
#
# THE PROFILE DIRECTORY MUST BE UNIQUE PER RUN. With a fixed --user-data-dir, a
# headless Edge that has not fully exited still holds its lockfile, and the next
# launch quietly attaches to that instance and returns WITHOUT writing the PDF.
# The script then reports "PDF step failed" over a document that is perfectly
# fine, which is a confusing way to lose half an hour. $TMP is mktemp -d with a
# trap, so it is fresh every run and cleaned up afterwards. Edge writes noisy
# task-manager and sync errors to stderr even on success; the exit code and the
# output file are what matter.
EDGE="/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
if [ -x "$EDGE" ]; then
  ABS=$(pwd -W 2>/dev/null || pwd)
  rm -f "$PDF"
  "$EDGE" --headless=new --disable-gpu --no-first-run --no-pdf-header-footer \
          --user-data-dir="$TMP/edge" --virtual-time-budget=30000 \
          --print-to-pdf="$ABS/$PDF" \
          "file:///$ABS/$OUT" >/dev/null 2>&1 || true
  if [ -f "$PDF" ]; then
    echo "wrote $PDF: $(wc -c < "$PDF") bytes"
  else
    echo "PDF step failed; the HTML above is still good"
  fi
else
  echo "Edge not found, skipping the PDF"
fi
#
# IF THE PDF STEP REPORTS A FAILURE. It launches a browser, and a sandboxed
# shell can block Edge from writing the output file. Edge still exits 0 and
# prints nothing, so the only symptom is a missing PDF -- it is not a problem
# with the document, and the HTML beside it is fine. Two ways to tell:
#
#   * it also fails on a trivial one-line HTML file, and
#   * the same command run from an ordinary terminal or from PowerShell works.
#
# PowerShell fallback, which does not go through this shell's sandbox:
#
#   $edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
#   $dir  = "<this folder>"
#   $prof = Join-Path $env:TEMP ("edgepdf_" + [guid]::NewGuid().ToString("N"))
#   Start-Process $edge -Wait -WindowStyle Hidden -ArgumentList @(
#     "--headless=new","--disable-gpu","--no-first-run","--no-pdf-header-footer",
#     "--user-data-dir=$prof","--virtual-time-budget=30000",
#     "--print-to-pdf=$dir\LoopValidation.pdf",
#     "file:///$($dir -replace '\\','/')/LoopValidation.html")
#
# A stray headless Edge left over from an earlier run is a separate failure with
# the same symptom -- see the profile-directory note above. `Get-Process msedge`
# will show it; stop it before retrying.
