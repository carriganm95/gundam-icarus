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

import pandas as pd
import uproot
import ROOT as r

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from makePrettyPlots import ROOT_FONT, add_watermark_root, make_pixel_canvas, style_axis

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
        t1.DrawLatex(0.65, 0.66, f"Fit \mu = {f1.GetParameter(1):.3f}")
        t1.DrawLatex(0.65, 0.62, f"Fit \sigma = {f1.GetParameter(2):.3f}")

        canvas.Update()
        add_watermark_root(canvas, pot=None)
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

    cv_path = '/exp/icarus/app/users/micarrig/gundam/gundam-icarus/configsNue/data/asimov/fit_asimov.root'
    stats, hToy = read_fit(cv_path, hist_path)

    # Overlay every toy histogram on one canvas
    x_title = hists[0].GetXaxis().GetTitle() or args.variable
    canvas = make_pixel_canvas("cToys", 1200, 800, 180, 60, 120, 120)
    r.SetOwnership(canvas, False)
    canvas.cd()
    l1 = r.TLegend(0.75, 0.65, 0.9, 0.75)
    l1.SetBorderSize(0)
    # l1.SetTextSize(0.02)

    t1 = r.TLatex()
    t1.SetNDC()
    t1.SetTextFont(ROOT_FONT)
    t1.SetTextSize(20)
    t1.DrawLatex(0.25, 0.7, f"NToys = {len(hists)}")

    ymax = 1.15 * max(h.GetMaximum() for h in hists)
    for i, h in enumerate(hists):
        h.SetFillStyle(0)  # GUNDAM saves data histograms with a solid fill
        h.SetMarkerStyle(0)
        h.SetLineColorAlpha(r.kAzure + 2, 0.3 if len(hists) > 20 else 1.0)
        h.SetLineWidth(1)
        h.SetMaximum(ymax)
        h.SetMinimum(0)
        if i == 0:
            style_axis(h.GetXaxis(), 30, title_px=36)
            style_axis(h.GetYaxis(), 30, title_px=36, title_offset=1.8)
            h.GetXaxis().SetTitle(x_title)
            h.GetYaxis().SetTitle("Events")
            h.Draw("HIST")
            l1.AddEntry(h, "Toy Data", "l")
        else:
            h.Draw("HIST SAME")
    hToy.SetFillStyle(0)
    hToy.SetMarkerStyle(0)
    hToy.SetLineColor(r.kBlack)
    hToy.SetLineWidth(2)
    hToy.Draw("HIST SAME")
    l1.AddEntry(hToy, "Asimov Data", "l")
    l1.Draw("SAME")

    canvas.Update()
    add_watermark_root(canvas, pot=None)
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
