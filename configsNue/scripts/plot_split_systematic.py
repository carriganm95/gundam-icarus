"""Per-parameter systematic breakdown, for the cross section and (optionally) the
pre/post-fit reconstructed distributions.

Reads the ROOT files produced by `run_frac_uncertainty.sh --split <set>`:
one subdirectory per isolated parameter (via enableOnlyParameters), plus
'template_only' (the stat baseline) and the full-set reference (e.g.
'detector_only').

For each parameter p the isolated contribution to a fractional uncertainty is
sqrt(max(0, frac_p^2 - frac_stat^2)) per bin.

Outputs (to --output-dir), for the cross section and, with --with-dist, the
post-fit and pre-fit distributions:
  <set>_split[_dist|_distprefit]_ranking.png/pdf   parameters ranked by their
      mean isolated contribution (the "which parameter dominates" plot).
  <set>_split[_dist|_distprefit]_perbin.png/pdf    per-bin step curves for the
      top parameters, with the stat baseline and the full-set reference.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_frac_uncertainty import (parse_dist_edges, load_frac_unc, isolated,
                                    XSEC_HIST, DIST_HIST, DIST_PREFIT_HIST, HIST_FILE, VAR)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from makePrettyPlots import stamp_ax, POT_DEFAULT


def summarize(arr):
    """Scalar summary for ranking: mean isolated contribution over finite bins."""
    return float(np.nanmean(arr)) if np.any(np.isfinite(arr)) else 0.0


def _step(ax, edges, frac, **kw):
    ax.step(np.append(edges[:-1], edges[-1]), np.append(frac, frac[-1]),
            where="post", **kw)


def make_split_plots(data_dir, output_dir, set_name, baseline, reference, pot, top,
                     hist_key, file_name, suffix, xlabel, qlabel):
    """Ranking + per-bin plots for one quantity (xsec / post-fit / pre-fit dist)."""
    def load(name):
        p = os.path.join(data_dir, name, file_name)
        if not os.path.exists(p):
            return None, None
        return parse_dist_edges(p, hist_key), load_frac_unc(p, hist_key)

    edges, stat = load(baseline)
    if stat is None:
        print(f"SKIP {qlabel}: baseline '{baseline}/{file_name}' not found")
        return

    exclude = {baseline}
    ref = None
    if reference:
        exclude.add(reference)
        _, ref = load(reference)

    params = {}
    for name in sorted(os.listdir(data_dir)):
        if name in exclude or not os.path.isdir(os.path.join(data_dir, name)):
            continue
        e, frac = load(name)
        if frac is None:
            continue
        if edges is None:
            edges = e
        params[name] = isolated(frac, stat)
    if not params:
        print(f"SKIP {qlabel}: no parameter variants with {file_name}")
        return

    ranked = sorted(params.items(), key=lambda kv: summarize(kv[1]), reverse=True)
    tag = f"{set_name}_split" + (f"_{suffix}" if suffix else "")

    # ---- ranking bar chart ------------------------------------------------
    names = [k for k, _ in ranked]
    vals = [summarize(v) for _, v in ranked]
    fig, ax = plt.subplots(figsize=(8, max(3.0, 0.3 * len(names) + 1.5)))
    ypos = np.arange(len(names))[::-1]
    ax.barh(ypos, vals, color="steelblue")
    for y, v in zip(ypos, vals):
        ax.text(v, y, f" {v:.1f}", va="center", ha="left", fontsize=7)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel(f"mean isolated contribution to {qlabel} unc. [%]")
    ax.set_xlim(0, max(vals) * 1.15 if vals and max(vals) > 0 else 1)
    ax.set_title(f"{set_name}: per-parameter {qlabel} contribution",
                 fontsize=11, pad=30)
    stamp_ax(ax, pot)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(output_dir, f"{tag}_ranking.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {tag}_ranking.png/pdf")

    # ---- per-bin step plot of the top-N -----------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    _step(ax, edges, stat, label="Stat only (reference)", color="grey",
          ls="--", lw=2.0, alpha=0.7)
    if ref is not None:
        _step(ax, edges, ref, label=f"{reference} (full set)", color="black",
              ls="--", lw=2.0, alpha=0.7)
    cmap = plt.cm.tab10
    ntop = min(top, len(ranked))
    for i, (name, contrib) in enumerate(ranked[:ntop]):
        _step(ax, edges, contrib, label=name, color=cmap(i % 10), lw=1.7)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fractional uncertainty [%]")
    ax.legend(fontsize=8, frameon=False, ncol=2, loc="upper right")
    stamp_ax(ax, pot)
    ax.text(0.03, 0.97, f"{set_name}: top {ntop} ({qlabel})", transform=ax.transAxes,
            va="top", ha="left", fontsize=11, fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(output_dir, f"{tag}_perbin.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {tag}_perbin.png/pdf")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--set", default="systematic", help="set name, used in titles/filenames")
    ap.add_argument("--baseline", default="template_only", help="stat-only variant")
    ap.add_argument("--reference", default=None, help="full-set variant, drawn as reference")
    ap.add_argument("--pot", default=POT_DEFAULT)
    ap.add_argument("--top", type=int, default=8, help="how many top params in the per-bin plot")
    ap.add_argument("--with-dist", action="store_true",
                    help="also make the pre/post-fit reconstructed-distribution breakdowns")
    ap.add_argument("--var", default=None, help="Analysis variable (default: $GUNDAM_VAR, else EleE)")
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # (hist_key, filename-suffix, x-axis label, quantity label)
    quantities = [(XSEC_HIST, "", VAR["true_label"], "cross section")]
    if args.with_dist:
        quantities += [
            (DIST_HIST,        "dist",       VAR["reco_label"], "post-fit distribution"),
            (DIST_PREFIT_HIST, "distprefit", VAR["reco_label"], "pre-fit distribution"),
        ]

    for hist_key, suffix, xlabel, qlabel in quantities:
        make_split_plots(args.data_dir, args.output_dir, args.set, args.baseline,
                         args.reference, args.pot, args.top,
                         hist_key, HIST_FILE[hist_key], suffix, xlabel, qlabel)


if __name__ == "__main__":
    main()
