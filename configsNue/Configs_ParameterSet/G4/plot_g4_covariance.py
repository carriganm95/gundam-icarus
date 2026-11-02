"""Styled plots of a G4 reinteraction covariance from calc_geant4_covariance.py.

Reads <job>/<job>_{covariance,correlation,univ_yields,binning} from the output
ROOT file and writes, in the analysis style used by makePrettyPlots.py:

  CovMat_<job>.png/pdf, CorrMat_<job>.png/pdf
      ROOT matrices (kBird / kThermometer, ICARUS watermark), bins labelled by
      their edges and the axes titled with the binned variable.
  G4Dist_<job>.png/pdf   (with --config, or called from calc_geant4_covariance.py)
      Fine-binned distribution of the binned variable: nominal (CV) weights
      with the ±1σ band of the G4 universes, ratio to CV below, GUNDAM
      parameter-bin edges marked.
  UnivDist_<job>_bin<i>.png/pdf
      1D distribution of universe yield / CV in bin i.

Usage: python plot_g4_covariance.py outputs/output_CovMat_<job>.root --job <job> [--outdir plots]
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))  # configsNue
from makePrettyPlots import (_style, stamp_ax, add_watermark_root, style_axis,
                             make_pixel_canvas, ROOT_FONT)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT as r

r.gROOT.SetBatch(True)

# Binned variable -> (axis title for matplotlib, axis title for ROOT, short symbol, unit)
VARIABLES = {
    "true_leading_proton_p": ("True leading proton momentum [MeV/c]",
                              "True leading proton momentum [MeV/c]", "p", "MeV/c"),
    "true_leading_pion_ke": ("True leading pion kinetic energy [MeV]",
                             "True leading pion kinetic energy [MeV]", "KE", "MeV"),
}
# Particle whose reinteraction universes are varied, from the job name
PARTICLES = [("PiPlus", r"$\pi^{+}$", "#pi^{+}"), ("PiMinus", r"$\pi^{-}$", "#pi^{-}"),
             ("Proton", "proton", "proton")]
OPEN_EDGE = 1e4  # an upper edge at or above this is an open-ended last bin


def read_binning(path):
    var, edges = None, []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue
            if fields[0].startswith("variables"):
                var = fields[1]
                continue
            lo, hi = float(fields[0]), float(fields[1])
            if not edges:
                edges.append(lo)
            edges.append(hi)
    return var, np.array(edges)


def bin_labels(edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append(f"> {lo:g}" if hi >= OPEN_EDGE else f"{lo:g}-{hi:g}")
    return out


def particle_of(job):
    for key, mpl, root in PARTICLES:
        if key in job:
            return mpl, root
    return job, job


def save(fig, outdir, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Wrote", os.path.join(outdir, name + ".png/pdf"))


def plot_matrix(m, labels, axis_title, z_title, palette, text_fmt, zrange, outdir, name, header):
    n = len(labels)
    h = r.TH2D(name, "", n, 0, n, n, 0, n)
    h.SetDirectory(0)
    for i in range(n):
        for j in range(n):
            h.SetBinContent(i + 1, j + 1, m[i][j])
    for i, lbl in enumerate(labels, start=1):
        h.GetXaxis().SetBinLabel(i, lbl)
        h.GetYaxis().SetBinLabel(i, lbl)

    r.gStyle.SetPalette(palette)
    r.gStyle.SetOptTitle(0)
    r.gStyle.SetOptStat(0)
    r.gStyle.SetPaintTextFormat(text_fmt)
    r.gStyle.SetTextFont(42)  # cell values: regular, not bold

    frame, px, z_px = 1200, 48, 48
    left, bottom = 400, 250
    c = make_pixel_canvas(f"c_{name}", frame, frame, left, 330, bottom, 140)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        style_axis(ax, px, title_px=56)
        ax.SetTitle(axis_title)
        ax.CenterTitle()
    h.GetXaxis().SetTitleOffset(1.6)
    h.GetYaxis().SetTitleOffset(2.9)
    h.GetZaxis().SetTitle(z_title)
    style_axis(h.GetZaxis(), z_px, title_px=56, title_offset=2.0)
    h.GetZaxis().SetMaxDigits(3)
    if zrange:
        h.SetMinimum(zrange[0])
        h.SetMaximum(zrange[1])
    h.SetMarkerSize(1.6)  # cell-value text size
    h.SetMarkerColor(r.kBlack)
    h.Draw("COLZ TEXT")
    c.Update()

    pal = h.GetListOfFunctions().FindObject("palette")
    if pal:
        w = c.GetWw()
        pal.SetX1NDC((left + frame + 30) / w)
        pal.SetX2NDC((left + frame + 100) / w)
        pal.SetY1NDC(c.GetBottomMargin())
        pal.SetY2NDC(1.0 - c.GetTopMargin())

    # Which set this is, above the frame on the right (watermark is on the left)
    lat = r.TLatex()
    lat.SetNDC()
    lat.SetTextFont(ROOT_FONT)
    lat.SetTextSize(48)
    lat.SetTextAlign(31)
    lat.DrawLatex(1.0 - c.GetRightMargin(), 1.0 - c.GetTopMargin() + 0.008, header)

    c.Modified()
    c.Update()
    add_watermark_root(c, pot=None)
    for ext in ("png", "pdf"):
        c.SaveAs(os.path.join(outdir, f"{name}.{ext}"))


# Fine binning for the 1D distribution plot; chosen so every GUNDAM
# parameter-bin edge (90/210/420 MeV, 400/600/850 MeV/c) is a fine-bin edge.
FINE_BINS = {
    "true_leading_pion_ke": np.arange(0, 990 + 1, 30.0),
    "true_leading_proton_p": np.arange(0, 1500 + 1, 50.0),
}


def plot_distribution(df, all_RW, cfg, binning_file, job, outdir):
    """CV distribution of the binned variable with the ±1σ band of the G4 universes.

    df/all_RW/cfg are as in calc_geant4_covariance.py (FinalWeight = CV weight).
    The last fine bin includes overflow.
    """
    var, param_edges = read_binning(binning_file)
    xlabel_mpl, _, _, _ = VARIABLES.get(var, (var, var, var, ""))
    p_mpl, _ = particle_of(job)
    fine = FINE_BINS.get(var)
    if fine is None:
        fine = np.linspace(param_edges[0], param_edges[-2] * 2, 31)

    x = df[var].values
    w = df["FinalWeight"].values
    sel = ~np.isnan(x) & (x >= fine[0])
    print(f"[{job}] {np.sum(~sel)} of {len(x)} events have no {var} (no dial) and are not plotted")
    idx = np.clip(np.digitize(x[sel], fine) - 1, 0, len(fine) - 2)  # overflow -> last bin
    w_sel, rw_sel = w[sel], all_RW[sel]
    nb = len(fine) - 1
    cv = np.bincount(idx, weights=w_sel, minlength=nb)
    univ = np.zeros((all_RW.shape[1], nb))
    for b in range(nb):
        m = idx == b
        univ[:, b] = w_sel[m] @ rw_sel[m]
    sd = univ.std(axis=0, ddof=1)
    frac = np.divide(sd, cv, out=np.zeros_like(sd), where=cv > 0)

    # Cross-check: CV summed inside each parameter bin = the yield used for the covariance
    pe = param_edges.copy()
    pe[-1] = np.inf
    check = [cv[(fine[:-1] >= lo) & (fine[:-1] < hi)].sum() for lo, hi in zip(pe[:-1], pe[1:])]
    print(f"[{job}] CV per parameter bin from fine histogram:", np.round(check, 3))

    plt.style.use(_style.get_style())
    stepy = lambda y: np.append(y, y[-1])
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    ax.fill_between(fine, stepy(cv - sd), stepy(cv + sd), step="post", color="#4c72b0",
                    alpha=0.35, linewidth=0, label=rf"G4 universes $\pm1\sigma$ ({all_RW.shape[1]})")
    ax.step(fine, stepy(cv), where="post", color="black", linewidth=1.8, label="Nominal (CV)")
    axr.fill_between(fine, stepy(1 - frac), stepy(1 + frac), step="post", color="#4c72b0",
                     alpha=0.35, linewidth=0)
    axr.axhline(1.0, color="black", linewidth=1.2)
    for e in param_edges[1:-1]:
        for a_ in (ax, axr):
            a_.axvline(e, color="grey", linestyle="--", linewidth=1)
    ax.plot([], [], color="grey", linestyle="--", linewidth=1, label="GUNDAM parameter bins")

    ax.set_xlim(fine[0], fine[-1])
    ax.set_ylim(0, cv.max() * 1.45)
    ax.set_ylabel("Selected events (weighted)")
    ax.legend(fontsize=12, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.88))
    ax.text(0.97, 0.95, f"G4 {p_mpl} reinteraction", transform=ax.transAxes,
            ha="right", va="top", fontsize=14, fontweight="bold")
    ax.text(0.97, 0.62, "last bin includes overflow", transform=ax.transAxes,
            ha="right", va="top", fontsize=10, color="grey")
    stamp_ax(ax, pot=None)
    rmax = max(0.05, 1.3 * frac[cv > 0].max())
    axr.set_ylim(1 - rmax, 1 + rmax)
    axr.set_ylabel("Universe / CV")
    axr.set_xlabel(xlabel_mpl)
    save(fig, outdir, f"G4Dist_{job}")
    return fine, cv, univ


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rootfile")
    ap.add_argument("--job", required=True, help="JobName used in calc_geant4_covariance.py")
    ap.add_argument("--outdir", default=os.path.join(HERE, "plots"))
    ap.add_argument("--config", default=None,
                    help="calc_geant4_covariance.py JSON config; also draws the CV vs universes distribution "
                         "(re-reads the input tree)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    job = args.job

    f = r.TFile.Open(args.rootfile)
    cov_t = f.Get(f"{job}/{job}_covariance")
    corr_t = f.Get(f"{job}/{job}_correlation")
    univ_t = f.Get(f"{job}/{job}_univ_yields")
    binning = f.Get(f"{job}/{job}_binning")
    if not univ_t or not binning:
        sys.exit(f"{args.rootfile} has no universe yields for {job}; rerun calc_geant4_covariance.py")
    n = cov_t.GetNrows()
    cov = np.array([[cov_t[i][j] for j in range(n)] for i in range(n)])
    corr = np.array([[corr_t[i][j] for j in range(n)] for i in range(n)])
    univ = np.array([[univ_t[u][i] for i in range(n)] for u in range(univ_t.GetNrows())])
    binning_file = binning.GetTitle()
    var, edges = read_binning(binning_file)
    f.Close()

    xlabel_mpl, xlabel_root, _, _ = VARIABLES.get(var, (var, var, var, ""))
    p_mpl, p_root = particle_of(job)
    labels = bin_labels(edges)

    # The CV is not stored separately in the matrix; recover it from the TH1D
    f = r.TFile.Open(args.rootfile)
    h_cv = f.Get(f"{job}/{job}_cv_yield")
    cv = np.array([h_cv.GetBinContent(i + 1) for i in range(n)])
    f.Close()
    ratio = univ / cv  # (n_univ, n_bins)

    # --- 2D matrices (ROOT, makePrettyPlots style) ---
    header = f"G4 {p_root} reinteraction"
    plot_matrix(cov * 1e3, labels, xlabel_root, "Fractional covariance [#times10^{-3}]",
                r.kBird, ".2f", None, args.outdir, f"CovMat_{job}", header)
    plot_matrix(corr, labels, xlabel_root, "Correlation", r.kThermometer,
                ".2f", (-1, 1), args.outdir, f"CorrMat_{job}", header)

    # --- 1D plots (matplotlib, spineplot style) ---
    plt.style.use(_style.get_style())
    title = f"G4 {p_mpl} reinteraction"


    # Per-bin distributions
    _, _, sym, unit = VARIABLES.get(var, (var, var, var, ""))
    for i in range(n):
        rr = ratio[:, i]
        m, s = rr.mean(), rr.std(ddof=1)
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.hist(rr, bins=50, histtype="stepfilled", color="#4c72b0", alpha=0.35)
        ax.hist(rr, bins=50, histtype="step", color="#4c72b0", linewidth=1.5)
        ax.axvspan(m - s, m + s, color="C3", alpha=0.12, label=rf"$\pm1\sigma$ = {100 * s:.2f}%")
        ax.axvline(m, color="C3", linewidth=2, label=f"Mean = {m:.4f}")
        ax.axvline(1.0, color="black", linestyle="--", linewidth=1.5, label="CV")
        ax.set_xlabel("Universe yield / CV")
        ax.set_ylabel("Universes")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
        ax.legend(fontsize=11, frameon=False, loc="upper left")
        rng = f"{sym} > {edges[i]:g}" if edges[i + 1] >= OPEN_EDGE else f"{edges[i]:g} < {sym} < {edges[i + 1]:g}"
        ax.text(0.97, 0.95, f"{title}\n{rng} {unit}", transform=ax.transAxes,
                ha="right", va="top", fontsize=13, fontweight="bold")
        stamp_ax(ax, pot=None)
        save(fig, args.outdir, f"UnivDist_{job}_bin{i}")

    if args.config:
        import calc_geant4_covariance as calc
        cfg = calc.load_config(args.config)
        df, all_RW = calc.load_tree_and_universes(cfg)
        df = df.eval("FinalWeight = %s" % cfg["SampleWeightExpr"])
        plot_distribution(df, all_RW, cfg, binning_file, job, args.outdir)


if __name__ == "__main__":
    main()
