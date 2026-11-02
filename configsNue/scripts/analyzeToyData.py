"""Collect best-fit likelihoods and toy data histograms from gundamFitter outputs.

Reads every .root file under the input directory (including subdirectories,
e.g. the grid outputs/NN/ layout), and writes to the output directory:
  - toyLLH.pkl / toyLLH.csv : one row per fit with the LLH values and fit status
  - llh_<name>.root/.pdf/.png : one histogram per LLH component
  - toyHistograms.root      : canvas with every toy data histogram overlaid
                              (also saved as toyHistograms.pdf/.png)

Example:
  python3 scripts/analyzeToyData.py /pnfs/.../asimov/outputs -o toy_analysis
"""
import os
import sys
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
import uproot
import ROOT as r
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from makePrettyPlots import POT_DEFAULT, ROOT_FONT, add_watermark_root, make_pixel_canvas, style_axis

BEST_FIT_STATS = "FitterEngine/postFit/bestFitStats"
HIST_DIR = "FitterEngine/preFit/plots/histograms"
STATS_BRANCHES = [
    "toyIndex",
    "fitConverged",
    "fitStatusCode",
    "covStatusCode",
    "edmBestFit",
    "totalLikelihoodAtBestFit",
    "statLikelihoodAtBestFit",
    "penaltyLikelihoodAtBestFit",
]


def parse_args():
    parser = ArgumentParser(description="Analyze toy data")
    parser.add_argument("input_dir", help="Input directory with gundamFitter output files")
    parser.add_argument("-o", "--output_dir", default="toy_analysis", help="Output directory")
    parser.add_argument("--sample", default="signal_EleE", help="Sample name in the fitter output")
    parser.add_argument("--variable", default="reco_leading_electron_energy", help="Plotted variable")
    return parser.parse_args()


def read_fit(path, hist_path):
    """Return (dict of best-fit stats, toy data TH1D) for one fitter output."""
    path = Path(path)
    with uproot.open(path) as f:
        stats = f[BEST_FIT_STATS].arrays(STATS_BRANCHES, library="pd").iloc[0].to_dict()

    rootFile = r.TFile.Open(str(path))
    hToy = rootFile.Get(hist_path)
    if not hToy:
        rootFile.Close()
        raise KeyError(hist_path)
    hToy = hToy.Clone(f"hToy_{path.stem}")
    hToy.SetDirectory(0)  # keep the histogram after the file is closed
    rootFile.Close()
    return stats, hToy


def hist_values(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])


def toy_mean_chi2(hists, hAsimov):
    """Correlated chi2 of the toy mean vs Asimov, using the toy covariance.

    Returns {name: (chi2, ndof, pvalue)} for two tests:
      'mean':   uncertainty on the mean (cov / NToys) -> is the toy mean biased?
      'spread': full toy covariance -> is Asimov consistent with the toy spread?
    """
    toys = np.array([hist_values(h) for h in hists])
    resid = toys.mean(axis=0) - hist_values(hAsimov)
    cov = np.cov(toys, rowvar=False)
    ndof = resid.size
    out = {}
    for name, c in (("mean", cov / len(toys)), ("spread", cov)):
        chi2 = float(resid @ np.linalg.lstsq(c, resid, rcond=None)[0])
        out[name] = (chi2, ndof, float(stats.chi2.sf(chi2, ndof)))
    return out


def make_llh_hists(df, output_dir):
    """Save one canvas per LLH column (total, stat, penalty) to llh_<name>.*."""
    llh_cols = [
        ("totalLLH",   "Total LLH",   r.kAzure  + 2),
        ("statLLH",    "Stat LLH",    r.kRed    - 4),
        ("penaltyLLH", "Penalty LLH", r.kGreen  + 2),
    ]

    for col, title, color in llh_cols:
        vals = df[col].dropna().values
        if len(vals) == 0:
            continue
        vmin, vmax = float(vals.min()), float(vals.max())
        margin = 0.05 * (vmax - vmin) if vmax > vmin else 1.0
        nbins = max(20, min(50, len(vals) // 2))
        h = r.TH1D(f"h_{col}", f";{title};Toys", nbins, vmin - margin, vmax + margin)
        for v in vals:
            h.Fill(float(v))

        f1 = r.TF1(f"f1_{col}", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", vmin - margin, vmax + margin)
        f1.SetParameters(h.GetMaximum(), h.GetMean(), h.GetRMS())

        canvas = make_pixel_canvas(f"cLLH_{col}", 1200, 800, 180, 60, 120, 120)
        r.SetOwnership(canvas, False)
        canvas.cd()

        h.SetFillColorAlpha(color, 0.5)
        h.SetLineColor(color)
        style_axis(h.GetXaxis(), 30, title_px=36)
        style_axis(h.GetYaxis(), 30, title_px=36, title_offset=1.8)
        h.Draw("HIST")

        h.Fit(f1, "REMNQ")
        f1.SetLineColor(r.kBlack)
        f1.SetLineWidth(2)
        f1.Draw("SAME")

        t1 = r.TLatex()
        r.SetOwnership(t1, False)
        t1.SetNDC()
        t1.SetTextFont(ROOT_FONT)
        t1.SetTextSize(28)
        t1.DrawLatex(0.65, 0.78, f"NToys = {len(vals)}")
        t1.DrawLatex(0.65, 0.74, f"Mean = {h.GetMean():.3f}")
        t1.DrawLatex(0.65, 0.70, f"RMS = {h.GetRMS():.3f}")
        t1.DrawLatex(0.65, 0.66, f"Fit #mu = {f1.GetParameter(1):.3f}")
        # f1 is a Poisson in x/[2]: mean = [1], variance = [1]*[2] ([2] is a scale, not the width)
        fit_sigma = np.sqrt(abs(f1.GetParameter(1) * f1.GetParameter(2)))
        t1.DrawLatex(0.65, 0.62, f"Fit #sigma = {fit_sigma:.3f}")

        canvas.Update()
        add_watermark_root(canvas, pot=POT_DEFAULT)
        r.gStyle.SetImageScaling(3.)
        canvas.Modified()
        canvas.Update()

        canvas.SaveAs(os.path.join(output_dir, f"llh_{col}.png"))
        canvas.SaveAs(os.path.join(output_dir, f"llh_{col}.pdf"))
        fout = r.TFile(os.path.join(output_dir, f"llh_{col}.root"), "RECREATE")
        h.Write()
        canvas.Write()
        fout.Close()


def main():

    args = parse_args()
    r.gROOT.SetBatch(True)
    r.gStyle.SetOptStat(0)
    r.gStyle.SetOptTitle(0)

    hist_path = f"{HIST_DIR}/{args.sample}/{args.variable}/Data_TH1D"
    files = sorted(Path(args.input_dir).rglob("*.root"))
    print(f"Found {len(files)} .root files in {args.input_dir}")

    rows = []
    hists = []
    for path in files:
        try:
            stats, hToy = read_fit(path, hist_path)
        except Exception as e:  # failed jobs leave incomplete or unrelated files
            print(f"Skipping {path}: {type(e).__name__}: {str(e).splitlines()[0]}")
            continue
        stats["file"] = str(path)
        rows.append(stats)
        hists.append(hToy)

    if not rows:
        raise SystemExit("No usable fitter outputs found")

    df = pd.DataFrame(rows).rename(columns={
        "totalLikelihoodAtBestFit": "totalLLH",
        "statLikelihoodAtBestFit": "statLLH",
        "penaltyLikelihoodAtBestFit": "penaltyLLH",
    })
    print(f"Read {len(df)} fits ({int(df['fitConverged'].sum())} converged)")
    print(df[["totalLLH", "statLLH", "penaltyLLH"]].describe())

    os.makedirs(args.output_dir, exist_ok=True)
    df.to_pickle(os.path.join(args.output_dir, "toyLLH.pkl"))
    df.to_csv(os.path.join(args.output_dir, "toyLLH.csv"), index=False)

    make_llh_hists(df, args.output_dir)

    cv_path = '/exp/icarus/app/users/micarrig/gundam/gundam-icarus/configsNue/data/EleE/asimov/fit_asimov.root'
    stats, hToy = read_fit(cv_path, hist_path)

    # Overlay every toy histogram on one canvas
    x_title = hists[0].GetXaxis().GetTitle() or args.variable
    canvas = make_pixel_canvas("cToys", 1200, 800, 180, 60, 120, 120)
    r.SetOwnership(canvas, False)
    canvas.cd()
    l1 = r.TLegend(0.72, 0.62, 0.9, 0.78)
    l1.SetBorderSize(0)
    l1.SetFillStyle(0)

    toy_color = r.kAzure + 2
    ymax = 1.15 * max(h.GetMaximum() for h in hists)

    # Bin-by-bin mean of the toy data
    hMean = hists[0].Clone("hToyMean")
    hMean.Reset()
    for h in hists:
        hMean.Add(h)
    hMean.Scale(1.0 / len(hists))
    for i, h in enumerate(hists):
        h.SetFillStyle(0)  # GUNDAM saves data histograms with a solid fill
        h.SetMarkerStyle(0)
        h.SetLineColorAlpha(toy_color, 0.3 if len(hists) > 20 else 1.0)
        h.SetLineWidth(1)
        h.SetMaximum(ymax)
        h.SetMinimum(0)
        if i == 0:
            style_axis(h.GetXaxis(), 30, title_px=36)
            style_axis(h.GetYaxis(), 30, title_px=36, title_offset=1.8)
            h.GetXaxis().SetTitle(x_title)
            h.GetYaxis().SetTitle("Events")
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")

    hMean.SetFillStyle(0)
    hMean.SetMarkerStyle(0)
    hMean.SetLineColor(r.kRed + 1)
    hMean.SetLineWidth(3)
    hMean.Draw("HIST SAME")

    hToy.SetFillStyle(0)
    hToy.SetMarkerStyle(0)
    hToy.SetLineColor(r.kBlack)
    hToy.SetLineWidth(4)
    hToy.Draw("HIST SAME")

    # Opaque stand-in so the faint toy lines are readable in the legend
    hToyLeg = r.TH1D("hToyLeg", "", 1, 0, 1)
    r.SetOwnership(hToyLeg, False)
    hToyLeg.SetLineColorAlpha(toy_color, 1.0)
    hToyLeg.SetLineWidth(2)
    l1.AddEntry(hToyLeg, "Toy Data", "l")
    l1.AddEntry(hMean, "Toy Mean", "l")
    l1.AddEntry(hToy, "Asimov Data", "l")
    l1.Draw("SAME")

    t1 = r.TLatex()
    r.SetOwnership(t1, False)
    t1.SetNDC()
    t1.SetTextFont(ROOT_FONT)
    t1.SetTextSize(28)
    t1.DrawLatex(0.72, 0.56, f"NToys = {len(hists)}")

    # Toy mean vs Asimov: bias test (error on the mean) and consistency with the toy spread
    chi2s = toy_mean_chi2(hists, hToy)
    t1.SetTextSize(21)
    for name, (chi2, ndof, p) in chi2s.items():
        print(f"Toy mean vs Asimov ({name}): chi2/ndof = {chi2:.2f}/{ndof}, p = {p:.3g}")
    # Only the toy-spread test goes on the plot
    chi2, ndof, p = chi2s["spread"]
    t1.DrawLatex(0.72, 0.50, "Toy Mean vs Asimov:")
    t1.DrawLatex(0.72, 0.45, f"#chi^{{2}}/ndof = {chi2:.1f}/{ndof}, p = {p:.3f}")

    canvas.Update()
    add_watermark_root(canvas, pot=POT_DEFAULT)
    r.gStyle.SetImageScaling(3.)
    canvas.Modified()
    canvas.Update()

    outRoot = os.path.join(args.output_dir, "toyHistograms.root")
    fout = r.TFile(outRoot, "RECREATE")
    canvas.Write()
    fout.Close()
    canvas.SaveAs(os.path.join(args.output_dir, "toyHistograms.pdf"))
    canvas.SaveAs(os.path.join(args.output_dir, "toyHistograms.png"))

    print(f"Saved dataframe and canvas to {args.output_dir}/")


if __name__ == "__main__":

    main()
