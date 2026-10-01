"""Compose fig/e_nodal_vs_pertet.png -- the two renderings of |E|, side by side.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvpython.exe" make_nodal_comparison.py

Needs fig/e_xy_ring.png and fig/e_xy_ring_nodal.png, which make_plots.py writes
from the same slice and the same pinned colour range -- the ONLY difference
between them is which array the colour comes from.

WHY THIS FIGURE EXISTS. The nodal rendering is smoother everywhere and wrong in
one place, and the place is small enough that it is invisible at full-frame
scale. Describing it in prose invites the reader to take it on trust; a 3x crop
of the boundary settles it in one look. The box drawn on the full views is the
crop, so there is no question which part is being magnified.
"""
import os
from PIL import Image, ImageDraw

FIG = "fig"
BOX = (300, 180, 560, 360)      # the inner conductor/air boundary, upper left
PANEL = 560                     # width of each panel
PAD, HDR = 12, 30
INK = (15, 15, 15)


def load(name):
    p = os.path.join(FIG, name)
    if not os.path.isfile(p):
        raise SystemExit("missing %s -- run make_plots.py first" % p)
    return Image.open(p).convert("RGB")


def full(im):
    """The whole frame, with the crop marked."""
    im = im.copy()
    ImageDraw.Draw(im).rectangle(BOX, outline=(0, 0, 0), width=4)
    s = PANEL / im.size[0]
    return im.resize((PANEL, int(im.size[1] * s)), Image.LANCZOS)


def zoom(im):
    """The crop at 1:1 of the panel width, so the fringe is visible."""
    c = im.crop(BOX)
    s = PANEL / c.size[0]
    return c.resize((PANEL, int(c.size[1] * s)), Image.LANCZOS)


cell, node = load("e_xy_ring.png"), load("e_xy_ring_nodal.png")
tops = [full(cell), full(node)]
zooms = [zoom(cell), zoom(node)]

W = PANEL * 2 + PAD * 3
H = HDR + tops[0].size[1] + HDR + zooms[0].size[1] + PAD * 2
out = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(out)

for i, (im, lab) in enumerate(zip(tops, ("PER-TET  (cell data)  -- what the report uses",
                                         "NODAL  (point data)"))):
    x = PAD + i * (PANEL + PAD)
    d.text((x, 8), lab, fill=INK)
    out.paste(im, (x, HDR))

y = HDR + tops[0].size[1] + PAD
d.text((PAD, y + 6), "the boxed conductor / air boundary, magnified", fill=INK)
for i, im in enumerate(zooms):
    out.paste(im, (PAD + i * (PANEL + PAD), y + HDR))

path = os.path.join(FIG, "e_nodal_vs_pertet.png")
out.save(path)
print("wrote %s  %dx%d" % (path, out.size[0], out.size[1]))
