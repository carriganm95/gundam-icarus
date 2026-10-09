"""Compare fractional uncertainties across systematic subsets.

Reads xsec.root and toy_bf.root outputs from run_frac_uncertainty.sh and
produces two figures saved to --output-dir:

  xsec_frac_uncertainty.png/pdf
    Cross-section fractional uncertainty breakdown in true-energy bins.

  dist_frac_uncertainty.png/pdf
    Post-fit reco-distribution fractional uncertainty breakdown.

Both figures have:
  Left  — isolated contribution from each group (quadrature above stat)
  Right — isolated contribution of each group within the full set (N-1)

Individual contributions are quadrature differences:
  GENIE alone    = sqrt(max(0, xsec_only^2  - stat^2))
  GENIE in total = sqrt(max(0, all^2        - no_xsec^2))

Note: GUNDAM stores histogram bin edges as sequential integers; physical
energies come from XSEC_EDGES (xsec) or toy bin labels (distribution).
"""
import argparse
import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from makePrettyPlots import stamp_ax, POT_DEFAULT

import ROOT as r
r.gROOT.SetBatch(True)


XSEC_HIST        = "calcXsec/histograms/TrueEleE_XS_TH1D"
DIST_HIST        = "toyGen/histograms/signal_EleE_TH1D"
DIST_PREFIT_HIST = "toyGen/histograms/signal_EleE_pre_fit_TH1D"

HIST_FILE = {
    XSEC_HIST:        "xsec.root",
    DIST_HIST:        "toy_bf.root",
    DIST_PREFIT_HIST: "toy_pf.root",
}

# True-energy bin edges in GeV — matches binnings/binning_TrueEleE.txt (MeV/1000)
XSEC_EDGES = np.array([0.300, 0.390, 0.490, 0.600, 0.730, 0.870,
                        1.030, 1.240, 1.500, 2.000, 5.000])

COLORS = {
    "stat":     "grey",
    "genie":    "steelblue",
    "flux":     "darkorange",
    "detector": "forestgreen",
    "total":    "black",
}

LABELS_IND = {
    "stat":     "Stat only",
    "genie":    "GENIE",
    "flux":     "Flux",
    "detector": "Detector",
    "total":    "Total (reference)",
}

LABELS_N1 = {
    "stat":     "Stat only (reference)",
    "genie":    "GENIE/Total",
    "flux":     "Flux/Total",
    "detector": "Detector/Total",
    "total":    "Total",
}

LABELS_TRUE_N1 = {
    "total":    "Total (reference)",
    "stat":     "Stat only (reference)",
    "no_xsec":  "No GENIE",
    "no_flux":  "No Flux",
    "no_detector": "No Detector",
}

COLORS_TRUE_N1 = {
    "total":       "black",
    "stat":        "grey",
    "no_xsec":     "steelblue",
    "no_flux":     "darkorange",
    "no_detector": "forestgreen",
}


def parse_dist_edges(path, hist_key):
    """Parse physical reco-energy bin edges in GeV from a toy root file.

    GUNDAM stores bin edges as sequential integers; the actual energies are
    encoded in the bin labels as 'variable: [lo, hi['.
    """
    rf = r.TFile.Open(str(path))
    h = rf.Get(hist_key.split(";")[0])
    edges = []
    for i in range(1, h.GetNbinsX() + 1):
        label = h.GetXaxis().GetBinLabel(i)
        nums = label.split(":")[-1].strip().lstrip("[").split(",")
        lo = float(nums[0].strip())
        hi = float(nums[1].strip().rstrip("[]()"))
        if i == 1:
            edges.append(lo / 1000.0)
        edges.append(hi / 1000.0)
    rf.Close()
    return np.array(edges)


def load_frac_unc(path, hist_key):
    """Return fractional uncertainty (%) per bin from a GUNDAM output file."""
    with uproot.open(path) as f:
        h = f[hist_key]
        vals, _ = h.to_numpy()  # discard integer bin-index edges
        errs = h.errors()
    return np.where(vals > 0, errs / vals * 100.0, np.nan)


def load_all(data_dir, hist_key):
    """Load fractional uncertainty arrays for all variants. Returns (edges, raw_dict)."""
    variants = ["all", "template_only", "xsec_only", "flux_only",
                "detector_only", "no_xsec", "no_flux", "no_detector"]

    file_name = HIST_FILE[hist_key]
    is_xsec = hist_key == XSEC_HIST

    raw = {}
    edges = None
    for name in variants:
        fpath = os.path.join(data_dir, name, file_name)
        if not os.path.exists(fpath):
            print(f"WARNING: {fpath} not found, skipping {name}")
            continue
        raw[name] = load_frac_unc(fpath, hist_key)
        if edges is None and not is_xsec:
            edges = parse_dist_edges(fpath, hist_key)

    if is_xsec:
        edges = XSEC_EDGES
    return edges, raw


# Variants needed for each kind of plot.
IND_VARS = ["template_only", "xsec_only", "flux_only", "detector_only", "all"]
N1_VARS  = ["no_xsec", "no_flux", "no_detector"]


def have(raw, variants):
    """True if every named variant is present in raw."""
    return all(v in raw for v in variants)


def isolated(a, b):
    """Quadrature difference: isolated contribution of a above b."""
    return np.sqrt(np.maximum(0.0, a**2 - b**2))


def compute_individual(raw):
    """Isolated contribution of each group above stat (quadrature)."""
    stat  = raw["template_only"]
    total = raw["all"]
    return {
        "stat":     stat,
        "genie":    isolated(raw["xsec_only"],     stat),
        "flux":     isolated(raw["flux_only"],     stat),
        "detector": isolated(raw["detector_only"], stat),
        "total":    total,
    }


def compute_n_minus_1(raw):
    """Isolated contribution of each group within the full set (N-1)."""
    stat  = raw["template_only"]
    total = raw["all"]
    return {
        "stat":     stat,
        "genie":    isolated(total, raw["no_xsec"]),
        "flux":     isolated(total, raw["no_flux"]),
        "detector": isolated(total, raw["no_detector"]),
        "total":    total,
    }


def _draw_panel(ax, data, labels, edges):
    for key, frac in data.items():
        is_ref = key in ("stat", "total")
        ax.step(
            np.append(edges[:-1], edges[-1]),
            np.append(frac, frac[-1]),
            where="post",
            label=labels[key],
            color=COLORS[key],
            linestyle="--" if is_ref else "-",
            linewidth=2.0 if is_ref else 1.8,
            alpha=0.7 if is_ref else 1.0,
        )


def _style_ax(ax, edges, xlabel, ylabel, title, pot):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.legend(fontsize=11, frameon=False)
    stamp_ax(ax, pot)
    ax.text(0.03, 0.88, title,
            transform=ax.transAxes, fontsize=12,
            va="top", ha="left", fontweight="bold")


def make_n_minus_1_figure(raw, edges, output_path, xlabel, pot):
    """True N-1 plot: total uncertainty with each systematic group removed."""
    fig, ax = plt.subplots(figsize=(9, 6))

    draw_order = ["total", "no_xsec", "no_flux", "no_detector", "stat"]
    data = {
        "total":       raw["all"],
        "stat":        raw["template_only"],
        "no_xsec":     raw["no_xsec"],
        "no_flux":     raw["no_flux"],
        "no_detector": raw["no_detector"],
    }
    for key in draw_order:
        frac = data[key]
        is_ref = key in ("total", "stat")
        ax.step(
            np.append(edges[:-1], edges[-1]),
            np.append(frac, frac[-1]),
            where="post",
            label=LABELS_TRUE_N1[key],
            color=COLORS_TRUE_N1[key],
            linestyle="--" if is_ref else "-",
            linewidth=2.0 if is_ref else 1.8,
            alpha=0.7 if is_ref else 1.0,
        )

    _style_ax(ax, edges, xlabel, "Fractional uncertainty [%]",
              "True N−1", pot)

    for ext in ("png", "pdf"):
        fig.savefig(f"{output_path}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}.png/pdf")


def make_figure(individual, edges, output_path, xlabel, pot, n_minus_1=None):
    """Individual-contributions figure.

    With n_minus_1 given, draws the N-1 panel beside it (2-panel); otherwise
    draws the individual panel alone (single panel).
    """
    if n_minus_1 is None:
        fig, ax_ind = plt.subplots(figsize=(9, 6))
    else:
        fig, (ax_ind, ax_n1) = plt.subplots(1, 2, figsize=(16, 6))
        fig.subplots_adjust(wspace=0.3)

    _draw_panel(ax_ind, individual, LABELS_IND, edges)
    _style_ax(ax_ind, edges, xlabel, "Fractional uncertainty [%]",
              "Individual contributions", pot)

    if n_minus_1 is not None:
        _draw_panel(ax_n1, n_minus_1, LABELS_N1, edges)
        _style_ax(ax_n1, edges, xlabel, "Fractional uncertainty [%]",
                  "N−1 contributions", pot)

    for ext in ("png", "pdf"):
        fig.savefig(f"{output_path}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}.png/pdf")


def process_dataset(data_dir, output_dir, hist_key, name, xlabel, pot):
    """Load one histogram's variants and emit its figures.

    Always produces the individual-contributions figure when the individual
    variants are present. The N-1 panel and the separate true N-1 figure are
    produced only when the N-1 variants are also present (e.g. the normal run,
    but not the diagonal run where those variants are skipped).
    """
    edges, raw = load_all(data_dir, hist_key)
    if not have(raw, IND_VARS):
        missing = [v for v in IND_VARS if v not in raw]
        print(f"SKIP {name}: missing individual variants {missing}")
        return

    individual = compute_individual(raw)
    n_minus_1 = compute_n_minus_1(raw) if have(raw, N1_VARS) else None

    make_figure(individual, edges,
                os.path.join(output_dir, f"{name}_frac_uncertainty"),
                xlabel=xlabel, pot=pot, n_minus_1=n_minus_1)

    if n_minus_1 is not None:
        make_n_minus_1_figure(raw, edges,
                              os.path.join(output_dir, f"{name}_n_minus_1"),
                              xlabel=xlabel, pot=pot)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir",   required=True,
                   help="Base dir containing <variant>/xsec.root and toy_bf.root")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--pot", default=POT_DEFAULT)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    process_dataset(args.data_dir, args.output_dir, XSEC_HIST,
                    "xsec", r"$E_e$ [GeV]", args.pot)
    process_dataset(args.data_dir, args.output_dir, DIST_HIST,
                    "dist", r"Reco $E_e$ [GeV]", args.pot)
    process_dataset(args.data_dir, args.output_dir, DIST_PREFIT_HIST,
                    "dist_prefit", r"Reco $E_e$ [GeV]", args.pot)


if __name__ == "__main__":
    main()
