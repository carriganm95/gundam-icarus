"""Figures for gundam_fit_notes.tex.

Run from this directory:  python3 make_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

# ---------------------------------------------------------------------------
# Figure 1: one Theta_Delta2Npi event as stored in gundam_v30.root
# (knots at -1, -0.5, 0, 0.5, 1; weights symmetric about 0)
# ---------------------------------------------------------------------------
x = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
w = np.array([1.1789, 1.0895, 1.0, 1.0895, 1.1789])
xx = np.linspace(-1, 1, 401)

spline = CubicSpline(x, w, bc_type="not-a-knot")  # GUNDAM default "Spline"
graph = np.interp(xx, x, w)                       # GUNDAM "Graph"
w_lin = np.where(x < 0, 2.0 - np.interp(-x, x, w), w)  # w(-x) = 2 - w(x)

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), sharey=True)
for ax in axes:
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("dial value $x$")
axes[0].set_ylabel("event weight $w(x)$")

axes[0].plot(xx, spline(xx), color="C0")
axes[0].plot(x, w, "ko", ms=4)
axes[0].set_title("Stored knots, cubic Spline\n(slope at 0 is exactly 0)")

axes[1].plot(xx, graph, color="C1")
axes[1].plot(x, w, "ko", ms=4)
axes[1].axvspan(-1, 0, color="0.85", zorder=0)
axes[1].set_title(r"Graph (linear), bounds [0, 1]" "\n" r"(slope at 0 is $2\,[w(0.5)-1]$)")

axes[2].plot(xx, np.interp(xx, x, w_lin), color="C2")
axes[2].plot(x, w_lin, "ko", ms=4)
axes[2].set_title("Rewritten $w(-x)=2-w(x)$\n(one line through nominal)")

fig.tight_layout()
fig.savefig("fig_morph_dials.pdf")

# ---------------------------------------------------------------------------
# Figure 2: illustrative scan shapes (toy numbers, not from a fit)
# p in units of prior sigma, prior at 0
# ---------------------------------------------------------------------------
p = np.linspace(-3, 3, 601)
pen = p ** 2

cases = [
    ("(a) data constrains", 6.0 * (p - 0.4) ** 2),
    ("(b) data flat", 0.0 * p),
    ("(c) data linear (pull)", -1.6 * p),
    ("(d) data curves down", -0.93 * p ** 2),
]

fig, axes = plt.subplots(1, 4, figsize=(11, 3.1), sharey=True)
for ax, (title, stat) in zip(axes, cases):
    tot = stat + pen
    i0 = np.argmin(tot)
    tot_rel = tot - tot[i0]
    stat_rel = stat - stat[i0]
    ax.plot(p, pen, color="C3", ls="--", label="llhPenalty")
    ax.plot(p, stat_rel, color="C0", ls=":", label="llhStat")
    ax.plot(p, tot_rel, color="k", label="llh")
    ax.axhline(1, color="0.5", lw=0.8)
    inside = p[tot_rel <= 1]
    ax.axvspan(inside.min(), inside.max(), color="C2", alpha=0.15)
    ax.set_title(title)
    ax.set_xlabel(r"parameter (prior $\sigma$ units)")
    ax.set_ylim(-1, 8)
axes[0].set_ylabel(r"$\Delta$ value ($-2\ln L$ scale)")
axes[1].legend(loc="upper center", fontsize=8)
fig.tight_layout()
fig.savefig("fig_scan_shapes.pdf")

print("wrote fig_morph_dials.pdf, fig_scan_shapes.pdf")

# ---------------------------------------------------------------------------
# Figure 3: rectangle analogy -- the data fixes the product (area =
# distribution), not the split into length (cross section) x width (systematic)
# ---------------------------------------------------------------------------
from matplotlib.patches import Rectangle

A = 1.0  # fixed area = the measured distribution
fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 4.0))

# -- left: nested rectangles, same area, corners on the hyperbola l = A/w
ww = np.linspace(0.33, 3.0, 400)
axL.plot(ww, A / ww, color="0.4", ls="--", lw=1.2, zorder=1)
for (w, col) in [(0.5, "C0"), (1.0, "C1"), (2.0, "C2")]:
    axL.add_patch(Rectangle((0, 0), w, A / w, fill=False, edgecolor=col, lw=2.0))
    axL.plot(w, A / w, "o", color=col, ms=5)
axL.text(1.75, 1.75, "same area\n(distribution)", color="0.3", fontsize=9,
         ha="left", va="bottom")
axL.set_xlim(0, 2.6); axL.set_ylim(0, 2.6); axL.set_aspect("equal")
axL.set_xlabel(r"width $\propto$ systematic")
axL.set_ylabel(r"length $\propto$ cross section")
axL.set_title("Same area, different split\n(all fit the data equally)")

# -- right: data fixes the product tightly; a wide prior on the width leaves
#    the length (cross section) uncertain
axR.fill_between(ww, 0.99 * A / ww, 1.01 * A / ww, color="0.6", alpha=0.5,
                 zorder=2, label=r"measured area $\pm1\%$")
w0, sig = 1.0, 0.15
axR.axvspan(w0 * (1 - sig), w0 * (1 + sig), color="C0", alpha=0.18, zorder=1,
            label=r"width prior ($\pm15\%$)")
lo, hi = A / (w0 * (1 + sig)), A / (w0 * (1 - sig))
axR.annotate("", xy=(0.12, hi), xytext=(0.12, lo),
             arrowprops=dict(arrowstyle="<->", color="C3", lw=1.8))
axR.text(0.2, 0.5 * (lo + hi), "length\n(cross section)\nuncertainty\n$\\approx$ width prior",
         color="C3", fontsize=8.5, va="center")
axR.plot([w0 * (1 - sig), w0 * (1 + sig)], [hi, lo], "o", color="k", ms=4, zorder=3)
axR.set_xlim(0, 2.6); axR.set_ylim(0, 2.6); axR.set_aspect("equal")
axR.set_xlabel(r"width $\propto$ systematic")
axR.set_ylabel(r"length $\propto$ cross section")
axR.set_title("Data fixes the product, not the split")
axR.legend(loc="upper right", fontsize=8)

fig.tight_layout()
fig.savefig("fig_xsec_rectangle.pdf")

# ---------------------------------------------------------------------------
# Figure 4: spring analogy -- a systematic is held to nominal by a prior
# spring; the template has none, so it slides to absorb the degenerate shift
# ---------------------------------------------------------------------------
def draw_spring(ax, x0, x1, y, coils=9, amp=0.13, lw=1.6, color="C0"):
    lead = (x1 - x0) * 0.15
    cx = np.linspace(x0 + lead, x1 - lead, 2 * coils)
    cy = y + amp * np.array([1, -1] * coils)
    xs = np.concatenate([[x0, x0 + lead], cx, [x1 - lead, x1]])
    ys = np.concatenate([[y, y], cy, [y, y]])
    ax.plot(xs, ys, color=color, lw=lw)

fig, ax = plt.subplots(figsize=(9.2, 3.4))
y_sys, y_tmp = 1.0, 0.0
nominal = 0.0

# wall at nominal
ax.axvline(nominal, color="0.5", ls="--", lw=1.0)
ax.text(nominal, 1.6, "nominal", ha="center", color="0.4", fontsize=9)

# systematic: stiff prior spring, bead barely moves
d_sys = 0.3
draw_spring(ax, nominal, d_sys, y_sys, coils=11, amp=0.10, color="C0")
ax.plot(d_sys, y_sys, "o", color="C0", ms=14)
ax.text(-0.15, y_sys, "systematic", ha="right", va="center", fontsize=10)
ax.text(d_sys + 0.15, y_sys + 0.32,
        "stiff prior spring\n$\\Rightarrow$ barely moves ($\\sigma_{\\rm post}\\approx1$)",
        fontsize=8.5, va="center", color="C0")

# template: no prior spring, slides far to absorb the shift
d_tmp = 2.4
ax.plot([nominal, d_tmp], [y_tmp, y_tmp], color="0.75", ls=":", lw=1.4)
ax.plot(d_tmp, y_tmp, "o", color="C1", ms=14)
ax.text(-0.15, y_tmp, "template", ha="right", va="center", fontsize=10)
ax.text(d_tmp - 0.1, y_tmp + 0.34,
        "no prior spring\n$\\Rightarrow$ slides to absorb the shift",
        fontsize=8.5, va="center", ha="right", color="C1")

# the common degenerate shift both must accommodate
ax.annotate("", xy=(d_tmp, -0.55), xytext=(nominal, -0.55),
            arrowprops=dict(arrowstyle="->", color="C3", lw=1.8))
ax.text(0.5 * d_tmp, -0.78,
        r"required shift $\Delta$ (degenerate direction), split $\propto 1/$stiffness",
        ha="center", color="C3", fontsize=9)

ax.set_xlim(-0.9, 3.0); ax.set_ylim(-1.0, 1.8)
ax.set_xlabel("parameter displacement from nominal")
ax.set_yticks([])
ax.set_title("Stiff prior holds the systematic; the free template absorbs the shift")
fig.tight_layout()
fig.savefig("fig_xsec_spring.pdf")

print("wrote fig_xsec_rectangle.pdf, fig_xsec_spring.pdf")

# ---------------------------------------------------------------------------
# Figure 5: prior spring vs data spring -- how well a parameter is pinned is
# the balance of the prior spring (2/sigma^2) and the data spring (grows with
# event count). Template = remove the prior spring.
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3))
# (title, bead position x_eq in [0,1], post-fit half-width, data-spring coils, amp)
scenarios = [
    ("low stats / degenerate\nprior wins: $\\sigma_{\\rm post}\\approx$ prior",
     0.13, 0.21, 4, 0.14),
    ("balanced\n$\\sigma_{\\rm post}<$ prior",
     0.50, 0.15, 7, 0.10),
    ("high stats / distinguishable\ndata wins: $\\sigma_{\\rm post}\\ll$ prior",
     0.84, 0.08, 13, 0.07),
]
for ax, (title, xeq, hw, dcoils, damp) in zip(axes, scenarios):
    ax.grid(False)
    ax.axvline(0.0, color="0.55", lw=4)   # nominal / prior anchor
    ax.axvline(1.0, color="0.55", lw=4)   # data-preferred anchor
    draw_spring(ax, 0.0, xeq, 0.5, coils=10, amp=0.07, color="C0")   # prior (stiff)
    draw_spring(ax, xeq, 1.0, 0.5, coils=dcoils, amp=damp, color="C1")  # data (varies)
    ax.errorbar(xeq, 0.5, xerr=hw, fmt="o", color="k", ms=11,
                ecolor="C3", elinewidth=2.2, capsize=4, zorder=5)
    ax.set_title(title, fontsize=9)
    ax.set_xlim(-0.12, 1.12); ax.set_ylim(0, 1)
    ax.set_yticks([]); ax.set_xticks([0, 1])
    ax.set_xticklabels(["nominal", "data\nprefers"], fontsize=8.5)
# legend (drawn as proxy lines on the middle axis)
axes[1].plot([], [], color="C0", lw=1.8, label=r"prior spring ($2/\sigma_{\rm prior}^2$)")
axes[1].plot([], [], color="C1", lw=1.8, label=r"data spring ($\propto$ events)")
axes[1].plot([], [], color="C3", lw=2.2, label=r"post-fit width $\propto 1/\sqrt{k_{\rm prior}+k_{\rm data}}$")
axes[1].legend(loc="lower center", fontsize=7.3, framealpha=0.9)
fig.suptitle("A parameter is pinned by the balance of two springs "
             "(a template has no prior spring)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("fig_xsec_balance.pdf")

print("wrote fig_xsec_balance.pdf")
