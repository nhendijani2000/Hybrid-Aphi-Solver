"""Compose the two port-zoom figures into ONE image.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_port_pair.py

WHY THIS EXISTS. fig/e_port_zoom.png and fig/e_port_zoom_off.png are the same
view of the same tube cross-section on the same colour scale, taken ON the cut
and one tube radius off it. Each ALONE is an unremarkable blob; the measurement
is the difference between them. Placing them adjacent in the report is not the
same as putting them in one frame, where the eye compares them without being
asked to.

Both panels are cropped to their data area and share the single colour bar
lifted from one of them, so no colormap has to be re-matched in matplotlib --
the bar is the same pixels ParaView drew.
"""
import os
import numpy as np
import vtk
from paraview.numpy_support import vtk_to_numpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = "fig"


def load(path):
    r = vtk.vtkPNGReader()
    r.SetFileName(path)
    r.Update()
    im = r.GetOutput()
    d = im.GetDimensions()
    a = vtk_to_numpy(im.GetPointData().GetScalars()).reshape(d[1], d[0], -1)[::-1]
    return a[:, :, :3].astype(np.uint8)


def split(a):
    """Separate the rendered data area from the colour bar.

    THE BAR IS DRAWN OVER THE DATA, not beside it -- the slice fills the frame
    edge to edge, so there is no white gutter to find and an automatic split
    returns nothing. The bar's position is known instead: make_plots.py puts it
    at Position 0.86 with ScalarBarLength 0.42, measured from the BOTTOM, so in
    top-down rows it spans 1 - 0.72 to 1 - 0.28 of the height.
    """
    h, w, _ = a.shape
    x_split = int(0.845 * w)
    data = a[:, :x_split]
    # Tight to the ramp itself. make_plots.py sets Position 0.30 with
    # ScalarBarLength 0.42, measured from the BOTTOM, so the coloured strip is
    # rows 1-0.72 to 1-0.30 top-down. A looser crop samples the green field
    # above and below it and puts flat green caps on the rebuilt colour bar.
    bar = a[int(0.287 * h):int(0.693 * h), x_split:]
    return data, bar


on = load(os.path.join(FIG, "e_port_zoom.png"))
off = load(os.path.join(FIG, "e_port_zoom_off.png"))
on_data, on_bar = split(on)
off_data, _ = split(off)
print("  data panels: %s and %s; bar %s"
      % (on_data.shape[:2], off_data.shape[:2],
         on_bar.shape[:2] if on_bar is not None else "none"))

fig = plt.figure(figsize=(11.0, 5.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.22], wspace=0.04)

ax1 = fig.add_subplot(gs[0, 0])
ax1.imshow(on_data)
ax1.set_title("ON the cut  ($\\theta = 0$)\nair around the conductor reaches "
              "1$-$2 V/m", fontsize=11)

ax2 = fig.add_subplot(gs[0, 1])
ax2.imshow(off_data)
ax2.set_title("one tube radius OFF the cut\nthe same air peaks near 0.2 V/m",
              fontsize=11)

for ax in (ax1, ax2):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

# A clean colour bar, drawn by matplotlib but using ParaView's OWN ramp.
#
# Cropping ParaView's bar out of the render does not work: the bar is drawn
# OVER the slice, so the crop carries the green field behind it, and a colour
# key with data showing through it is the last thing this figure needs. Instead
# the ramp is SAMPLED from that crop -- the saturated strip down its left side --
# and rebuilt as a ListedColormap. That matches the panels exactly without
# having to know or re-specify which ParaView preset produced them.
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.cm import ScalarMappable

bar_rgb = on_bar.astype(float)
sat = bar_rgb.max(axis=2) - bar_rgb.min(axis=2)
col = int(np.argmax(sat.mean(axis=0)))          # the most colourful column
strip = bar_rgb[:, max(col - 1, 0):col + 2].mean(axis=1) / 255.0
strip = strip[::-1]                              # bar runs high->low top->bottom
cmap = ListedColormap(strip)

axc = fig.add_subplot(gs[0, 2])
cb = fig.colorbar(ScalarMappable(norm=LogNorm(5e-3, 2.0), cmap=cmap),
                  cax=axc)
cb.set_label("|E|   V/m", fontsize=10)
cb.ax.tick_params(labelsize=9)

fig.suptitle("|E| at the port cross-section, same view and same log scale",
             fontsize=13, y=0.99)
fig.text(0.5, 0.015,
         "Blue disc = the conductor, uniform at $J/\\sigma$ = 8.6e-03 V/m "
         "(no skin effect at $a/\\delta$ = 0.086).  The difference between the "
         "two panels IS the port:\n"
         "$\\Phi$ jumps 0.354 mV across a zero-thickness interface on the left, "
         "and is continuous on the right.",
         ha="center", fontsize=9.5)
fig.subplots_adjust(top=0.84, bottom=0.14)
out = os.path.join(FIG, "e_port_pair.png")
fig.savefig(out, dpi=150)
print("wrote %s" % out)
