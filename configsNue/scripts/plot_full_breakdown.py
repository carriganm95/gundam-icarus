"""Full cross-section uncertainty breakdown, with the stat-only line decomposed.

Reads the xsec.root files produced by `run_frac_uncertainty.sh --decompose`:

  group variants (full calcXsec config):
    all, template_only, xsec_only, flux_only, detector_only
  stat decomposition (calcXsec re-run on the template_only fit with one toggle):
    template_base   template post-fit (Hesse) covariance only, no throws/norms
    statthrow       + enableStatThrowInToys   (data Poisson)
    mcthrow         + enableEventMcThrow       (efficiency / MC statistics)
    fluxnorm        + "Constant flux" normalisation (fitSamples_True<var>.yaml)
    targets         + "Number of targets" normalisation

Each contribution is isolated as sqrt(max(0, a^2 - b^2)); together they add in
quadrature to Total:

  Total^2 = Template^2 + Stat^2 + EffStat^2 + FluxNorm^2 + Targets^2
          + GENIE^2 + FluxBeamline^2 + Detector^2
"""
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_frac_uncertainty import parse_dist_edges, load_frac_unc, isolated, XSEC_HIST, VAR

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from makePrettyPlots import stamp_ax, POT_DEFAULT


# label, color, dashed?
STYLE = {
    "template": ("Template (Hesse)",        "purple",       False),
    "stat":     ("Stat (data Poisson)",     "grey",         False),
    "effstat":  ("Efficiency stat (MC)",    "saddlebrown",  False),
    "fluxnorm": ("Flux norm",               "gold",         False),
    "targets":  ("Targets",                 "pink",         False),
    "genie":    ("GENIE",                   "steelblue",    False),
    "fluxbeam": ("Flux beamline",           "darkorange",   False),
    "detector": ("Detector",                "forestgreen",  False),
    "total":    ("Total",                   "black",        True),
}
# order the lines are drawn / legended
ORDER = ["total", "template", "stat", "effstat", "fluxnorm", "targets",
         "genie", "fluxbeam", "detector"]


def load(data_dir, name):
    path = os.path.join(data_dir, name, "xsec.root")
    if not os.path.exists(path):
        return None
    return load_frac_unc(path, XSEC_HIST)


def load_edges(data_dir, name):
    return parse_dist_edges(os.path.join(data_dir, name, "xsec.root"), XSEC_HIST)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--pot", default=POT_DEFAULT)
    ap.add_argument("--drop-targets", action="store_true",
                    help="omit the (usually negligible) targets line")
    ap.add_argument("--var", default=None, help="Analysis variable (default: $GUNDAM_VAR, else EleE)")
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tpl_full = load(args.data_dir, "template_only")   # full stat-only (systematic baseline)
    base = load(args.data_dir, "template_base")       # template Hesse only (decomp baseline)
    if tpl_full is None or base is None:
        sys.exit("need both 'template_only' and 'template_base' xsec.root in data-dir")

    data = {
        "total":    load(args.data_dir, "all"),
        "template": base,
        "stat":     isolated(load(args.data_dir, "statthrow"), base),
        "effstat":  isolated(load(args.data_dir, "mcthrow"),   base),
        "fluxnorm": isolated(load(args.data_dir, "fluxnorm"),  base),
        "targets":  isolated(load(args.data_dir, "targets"),   base),
        "genie":    isolated(load(args.data_dir, "xsec_only"),     tpl_full),
        "fluxbeam": isolated(load(args.data_dir, "flux_only"),     tpl_full),
        "detector": isolated(load(args.data_dir, "detector_only"), tpl_full),
    }
    edges = load_edges(args.data_dir, "template_only")

    fig, ax = plt.subplots(figsize=(10, 6.5))
    for key in ORDER:
        if key == "targets" and args.drop_targets:
            continue
        frac = data.get(key)
        if frac is None:
            print(f"WARNING: missing '{key}', skipping")
            continue
        label, color, dashed = STYLE[key]
        ax.step(np.append(edges[:-1], edges[-1]), np.append(frac, frac[-1]),
                where="post", label=label, color=color,
                ls="--" if dashed else "-", lw=2.2 if dashed else 1.8,
                alpha=0.8 if dashed else 1.0)

    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.set_xlabel(VAR["true_label"])
    ax.set_ylabel("Fractional uncertainty [%]")
    ax.legend(fontsize=10, frameon=False, ncol=2, loc="upper right")
    stamp_ax(ax, args.pot)
    ax.text(0.03, 0.97, "Decomposed contributions", transform=ax.transAxes,
            va="top", ha="left", fontsize=12, fontweight="bold")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.output_dir, f"xsec_full_breakdown.{ext}"),
                    bbox_inches="tight")
    plt.close(fig)
    print("Saved xsec_full_breakdown.png/pdf")

    # quadrature closure check
    comps = ["template", "stat", "effstat", "fluxnorm", "targets",
             "genie", "fluxbeam", "detector"]
    quad = np.sqrt(np.nansum([np.nan_to_num(data[c])**2 for c in comps if data[c] is not None], axis=0))
    tot = data["total"]
    if tot is not None:
        print("quadrature-sum vs Total (mean %): "
              f"{np.nanmean(quad):.1f} vs {np.nanmean(tot):.1f}")


if __name__ == "__main__":
    main()
