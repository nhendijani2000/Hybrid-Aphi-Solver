"""Copy the solver's own post-processing pictures into fig/ for the report.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvpython.exe" collect_postprocess_figs.py [output-dir]

This is deliberately a copy and not a re-render. The sixteen PP*.png in the
report are the ones `tools/postprocess.py` produced from the manifest the solver
wrote -- that is the whole claim of the case, so the report must show those files
and not a second set drawn by a script that happens to agree with them.

It is plain Python (no paraview import), so pvpython is only a convenient
interpreter here; any Python 3 will do.
"""
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

OUT = sys.argv[1] if len(sys.argv) > 1 else "output"
FIG = "fig"
EXPECTED = 16                # PP1 .. PP16, one per [postprocess] section

if not os.path.isdir(FIG):
    os.makedirs(FIG)

names = sorted(n for n in os.listdir(OUT)
               if re.match(r"^PP\d+_.*\.png$", n))
if not names:
    raise SystemExit("no PP*.png in %s -- run tools/postprocess.py first" % OUT)

# Sort by the request number so the report's order matches the input file's.
names.sort(key=lambda n: int(re.match(r"^PP(\d+)_", n).group(1)))

for n in names:
    shutil.copy2(os.path.join(OUT, n), os.path.join(FIG, n))
    print("  fig/%s" % n)

got = set(int(re.match(r"^PP(\d+)_", n).group(1)) for n in names)
missing = sorted(set(range(1, EXPECTED + 1)) - got)
print("  copied %d of %d" % (len(names), EXPECTED))
if missing:
    # Not fatal -- the input file may legitimately have been edited, which is
    # what the case is for. But say so, because the report names each figure and
    # a silently absent one becomes a broken image in the built document.
    print("  WARNING: no picture for PP%s" % ", PP".join(str(i) for i in missing))
