import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import os
import re
import sys
import uproot
from scipy import stats

sys.path.insert(0, '/nashome/m/micarrig/icarus/nuESpine/medulla')
from spineplot.style import stamp, load_style_from_toml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_vars  # per-variable names/labels; pick with GUNDAM_VAR or --var

POT_DEFAULT = "2.5e20"
STYLE_NAME = 'nue_icarus'
WORK_IN_PROGRESS = True

# --- MC (CAF) truth cross-section reference --------------------------------
# gundamCalcXsec builds dsigma/dE for the signal from the true-signal tree
# below, normalising by the flux and target count in
# Configs_FitSampleSet/fitSamples_True<var>.yaml (selection: true_fiducial_cut).
# We reproduce that exactly here so a *fake-data* truth curve — the nominal
# truth times whatever truth reweight was injected to build the fake data —
# lands on the same absolute scale as the GUNDAM extraction. (Verified: the
# nominal reweight reproduces the Asimov GUNDAM xsec to 1.000 in every bin.)
MC_TRUTH_FILE = "/exp/icarus/data/users/micarrig/gundam_v31.root"
MC_TRUTH_TREE = "events/NuMIFull/all_signal"
# The truth branch and its unit conversion come from analysis_vars.current().
MC_TRUTH_SELECTION = "true_fiducial_cut == 1"
MC_TRUTH_NOMINAL_WEIGHT = "0.35231940 * ppfx_cv_weight"  # -> 2.5e20 POT
MC_TRUTH_FLUX = 4.41961285e+09      # (nue+nuebar)/cm^2 @ 2.5e20 POT
MC_TRUTH_NTARGETS = 4.72144364e+30  # N Ar nuclei

# Per-tag signal-truth reweight applied to the MC truth to build the fake-data
# truth reference. Expressions are evaluated with numpy against the truth-tree
# branches (use `where(cond, a, b)` for the YAML's C-style `cond ? a : b`, and
# `&`/`|` for boolean and/or). Only reweights that change the SIGNAL truth
# belong here:
#   * Q2, res/QE/MEC/DIS 0p5/1p5 are pure functions of true quantities applied to the
#     signal -> they change dsigma/dE and need the scaled reference.
#   * nc (true_category==4) and pAngle (true_category!=0) are reco-based and
#     hit <1% of the fiducial signal (which is ~99% true_category==0), so the
#     signal truth is unchanged -> reference stays nominal.
#   * gainHigh/gainLow/emb inject a different *reco* response (separate tree in
#     merged_v30.root) with no change to the truth -> reference stays nominal.
# Feed any other version at run time with --truth-reweight "<expr>".
TRUTH_REWEIGHTS = {
    "asimov": "1.0",
    "asimovGrid": "1.0",
    "10pct": "1.0",  # real data: compare to nominal GENIE
    "res0p5": "where(true_interaction_mode == 1, 0.5, 1.0)",
    "res1p5": "where(true_interaction_mode == 1, 1.5, 1.0)",
    "QE0p5": "where(true_interaction_mode == 0, 0.5, 1.0)",
    "QE1p5": "where(true_interaction_mode == 0, 1.5, 1.0)",
    "MEC0p5": "where(true_interaction_mode == 10, 0.5, 1.0)",
    "MEC1p5": "where(true_interaction_mode == 10, 1.5, 1.0)",
    "DIS0p5": "where(true_interaction_mode == 2, 0.5, 1.0)",
    "DIS1p5": "where(true_interaction_mode == 2, 1.5, 1.0)",
    "Q2": (
        "where((true_momentum_transfer < 0.0) | (true_momentum_transfer > 0.7), 1.0, "
        "1.0 - 0.7*(1.0 - (-0.6*true_momentum_transfer*(true_momentum_transfer-0.7)/0.1225 "
        "+ true_momentum_transfer*(true_momentum_transfer-0.35)/0.245))**2)"
    ),
    "nc": "1.0",
    "pAngle": "1.0",
    "gainHigh": "1.0",
    "gainLow": "1.0",
    "emb": "1.0",
}

# Short description of each tag's rescale, shown in the plot legend as
# "GENIE (<desc>)". Only tags whose truth is actually rescaled need an entry;
# everything else falls back to "GENIE (nominal)".
TRUTH_REWEIGHT_LABELS = {
    "Q2": "Q² reweight",
    "res0p5": "RES ×0.5",
    "res1p5": "RES ×1.5",
    "QE0p5": "QE ×0.5",
    "QE1p5": "QE ×1.5",
    "MEC0p5": "MEC ×0.5",
    "MEC1p5": "MEC ×1.5",
    "DIS0p5": "DIS ×0.5",
    "DIS1p5": "DIS ×1.5",
}

_style = load_style_from_toml(STYLE_NAME)
_wm_style = load_style_from_toml('nue_icarus_data') if WORK_IN_PROGRESS else _style


def stamp_ax(ax, pot=POT_DEFAULT):
    """Apply the spineplot watermark (and POT label) to a matplotlib axis. Call after limits are set."""
    stamp(ax, style_name=STYLE_NAME, pot=float(pot) if pot else None, work_in_progress=WORK_IN_PROGRESS)


def add_watermark_root(canvas, pot=POT_DEFAULT):
    """Stamp 'ICARUS Work In Progress' (top-left) and POT (top-right) just above the plot frame's top line."""
    canvas.cd()
    top_line = 1.0 - canvas.GetTopMargin()
    left = canvas.GetLeftMargin()
    right = 1.0 - canvas.GetRightMargin()

    # ROOT canvases can't use spineplot's stamp(), so mirror its label text and color.
    label = _wm_style.mark_preliminary.replace(r'$\bf{ICARUS}$', '#bf{ICARUS}')
    latex = r.TLatex()
    r.SetOwnership(latex, False)
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.SetTextColor(r.TColor.GetColor('#d67a11'))
    latex.SetTextAlign(11)  # left, bottom
    latex.DrawLatex(left, top_line + 0.005, label)
    if pot:
        mag = int(np.floor(np.log10(float(pot))))
        latex.SetTextColor(r.kBlack)
        latex.SetTextAlign(31)  # right, bottom
        latex.DrawLatex(right, top_line + 0.005, f"{float(pot) / 10**mag:.2f}#times10^{{{mag}}} POT")


BINNING_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'binnings')


def load_template_bin_labels(binning_file=None):
    """Build e.g. 'E_{e} [lo, hi] MeV' labels for each template parameter from a GUNDAM binning text file."""
    v = analysis_vars.current()
    binning_file = binning_file or os.path.join(BINNING_DIR, v['template_binning'])
    labels = []
    with open(binning_file) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 2 or fields[0].startswith('variables'):
                continue
            lo, hi = float(fields[0]), float(fields[1])
            labels.append(v['param_label'].format(lo=lo, hi=hi))
    return labels


G4_PARAM_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Configs_ParameterSet', 'G4')
# Particle name and binned variable for each G4 reinteraction set; the bin
# edges are read from the set's own binning file (parametersBinningPath).
G4_PARTICLE = {"Proton": "Proton", "PiPlus": "#pi^{+}", "PiMinus": "#pi^{-}"}
G4_VARIABLE = {"true_leading_proton_p": ("p", "MeV/c"),
               "true_leading_pion_ke": ("KE", "MeV")}


def _fmt_edge(x):
    return f"{x:g}"


def g4_bin_labels(set_name):
    """Labels like 'Proton p: 400-600 MeV/c' or '#pi^{+} KE: 90-210 MeV' for a G4 set's parameters."""
    particle = next((v for k, v in G4_PARTICLE.items() if f"G4 {k} " in set_name + " "), set_name)
    binning = None
    for fn in sorted(os.listdir(G4_PARAM_DIR)):
        if not fn.endswith('.yaml'):
            continue
        text = open(os.path.join(G4_PARAM_DIR, fn)).read()
        if re.search(rf'^name:\s*"?{re.escape(set_name)}"?\s*$', text, re.M):
            m = re.search(r'parametersBinningPath:\s*"?\$\{GUNDAM_CONFIG_DIR\}/(\S+?)"?\s*$', text, re.M)
            binning = os.path.join(os.path.dirname(os.path.abspath(__file__)), m.group(1))
            break
    if binning is None:
        return None
    labels = []
    with open(binning) as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue
            if fields[0].startswith('variables'):
                var, unit = G4_VARIABLE.get(fields[1], (fields[1], ""))
                continue
            lo, hi = float(fields[0]), float(fields[1])
            # Open-ended last bin (upper edge is a large sentinel value).
            rng = f"> {_fmt_edge(lo)}" if hi >= 1e4 else f"{_fmt_edge(lo)}-{_fmt_edge(hi)}"
            labels.append(f"{particle} {var}: {rng} {unit}".strip())
    return labels


def pretty_param_label(label, template_labels, template_mode=False):
    """Drop the '<parameter set>/' prefix and '#N_' index; name template parameters by their true-energy bin."""
    is_template = template_mode or label.startswith('Template Parameter')
    short = label.split('/')[-1]
    match = re.match(r'#(\d+)', short)
    if is_template and match and int(match.group(1)) < len(template_labels):
        return template_labels[int(match.group(1))]
    # Eigen-decomposed sets (e.g. G4) only have '#N_eigen' names: keep the set
    # name so the entries are distinguishable, e.g. 'G4 Proton eigen 0'.
    if match and short.endswith('_eigen') and '/' in label:
        set_short = label.split('/')[0].replace(' Reinteraction', '')
        return f"{set_short} eigen {match.group(1)}"
    return re.sub(r'^#\d+_', '', short)


def relabel_axis(axis, template_labels, template_mode=False):
    """Rename every bin label on an axis with pretty_param_label and return the new labels."""
    for i in range(1, axis.GetNbins() + 1):
        axis.SetBinLabel(i, pretty_param_label(axis.GetBinLabel(i), template_labels, template_mode))
    return [axis.GetBinLabel(i) for i in range(1, axis.GetNbins() + 1)]


# Layout for the ROOT plots is done in pixels (font precision 3, i.e. font 43)
# so label sizes and margins don't depend on the canvas aspect ratio.
ROOT_FONT = 43
_measure_canvas = None


def text_width_px(text, px):
    """Rendered width (px) of a TLatex string at text size px, measured on a scratch canvas."""
    global _measure_canvas
    if _measure_canvas is None:
        _measure_canvas = r.TCanvas("c_measure", "c_measure", 1000, 1000)
        r.SetOwnership(_measure_canvas, False)
    _measure_canvas.cd()
    t = r.TLatex(0, 0, text)
    t.SetTextFont(ROOT_FONT)
    t.SetTextSize(px)
    return t.GetXsize() * _measure_canvas.GetWw()


def max_text_width_px(labels, px):
    # Small safety factor: the TLatex measurement runs a little short at small sizes.
    return 1.1 * max([text_width_px(lbl, px) for lbl in labels] + [0])


def label_px(frame_px, n_bins, cap_frac=0.07):
    """Largest label size (px) that fits the per-bin spacing once rotated 45 degrees, capped as a fraction of the frame."""
    return int(min(cap_frac * frame_px, 0.6 * frame_px / n_bins))


def rotated_label_height(labels, px):
    """Vertical space (px) needed below an axis for labels rotated 45 degrees."""
    return int(0.71 * (max_text_width_px(labels, px) + px))


def rotated_label_overhang(labels, px, bin_px):
    """How far (px) the first rotated label reaches left of the frame edge."""
    return int(0.71 * (text_width_px(labels[0], px) + px) - bin_px / 2) if labels else 0


def style_axis(axis, px, title_px=None, title_offset=None, rotate=False):
    axis.SetLabelFont(ROOT_FONT)
    axis.SetLabelSize(px)
    if title_px:
        axis.SetTitleFont(ROOT_FONT)
        axis.SetTitleSize(title_px)
    if title_offset is not None:
        axis.SetTitleOffset(title_offset)
    if rotate:
        # Rotate 45 degrees, right-aligned so each label ends under its bin.
        for i in range(1, axis.GetNbins() + 1):
            axis.ChangeLabel(i, 45, -1, 33)


def make_pixel_canvas(name, frame_w, frame_h, left, right, bottom, top):
    """Canvas whose frame is exactly frame_w x frame_h pixels, with the given pixel margins."""
    width, height = left + frame_w + right, bottom + frame_h + top
    c = r.TCanvas(name, name, width, height)
    r.SetOwnership(c, False)
    c.SetLeftMargin(left / width)
    c.SetRightMargin(right / width)
    c.SetBottomMargin(bottom / height)
    c.SetTopMargin(top / height)
    return c


def plotMatrix(tree, path, outputDir, plot_name, z_title, colorMap="kBird", pot=POT_DEFAULT,
               hist=None, labels=None):
    """Draw a parameter matrix. Pass hist (and labels) to draw an already-built TH2D instead of reading path."""

    h = hist if hist is not None else tree.Get(path)
    # Detach from the source file's directory so TFile::Close() at the end
    # of main() doesn't delete it out from under the canvas that's still
    # displaying it (SetOwnership below only guards against Python's GC,
    # not ROOT's own directory-based auto-cleanup of TH1-derived objects).
    if hasattr(h, "SetDirectory"):
        h.SetDirectory(0)

    r.gStyle.SetPalette(getattr(r, colorMap))
    r.gStyle.SetOptTitle(0)
    r.gStyle.SetOptStat(0)

    if labels is not None:
        for i, lbl in enumerate(labels, start=1):
            h.GetXaxis().SetBinLabel(i, lbl)
            h.GetYaxis().SetBinLabel(i, lbl)
        x_labels = y_labels = list(labels)
    else:
        template_labels = load_template_bin_labels()
        template_mode = 'Template Parameter' in path
        x_labels = relabel_axis(h.GetXaxis(), template_labels, template_mode)
        y_labels = relabel_axis(h.GetYaxis(), template_labels, template_mode)

    # Square frame; margins sized to the labels and colour bar so there's no dead space.
    frame = max(1600, 28 * h.GetNbinsX())
    px = label_px(frame, h.GetNbinsX(), cap_frac=0.04)
    z_px = 60
    left = max(int(max_text_width_px(y_labels, px)) + 40,
               rotated_label_overhang(x_labels, px, frame / len(x_labels)) + 20)
    bottom = rotated_label_height(x_labels, px) + 40
    right = 400
    top = 120
    c = make_pixel_canvas(f"c_{plot_name}", frame, frame, left, right, bottom, top)

    h.SetTitle("")
    h.GetXaxis().SetTitle("")
    h.GetYaxis().SetTitle("")
    style_axis(h.GetXaxis(), px, rotate=True)
    style_axis(h.GetYaxis(), px)
    h.GetZaxis().SetTitle(z_title)
    style_axis(h.GetZaxis(), z_px, title_px=int(1.2 * z_px), title_offset=1.6)
    h.GetZaxis().SetMaxDigits(3)
    h.Draw("COLZ")
    c.Update()

    # Put the palette just right of the frame, leaving room for its labels and title.
    palette = h.GetListOfFunctions().FindObject("palette")
    if palette:
        width = c.GetWw()
        palette.SetX1NDC((left + frame + 30) / width)
        palette.SetX2NDC((left + frame + 110) / width)
        palette.SetY1NDC(c.GetBottomMargin())
        palette.SetY2NDC(1.0 - c.GetTopMargin())
    c.Modified()
    c.Update()
    r.gStyle.SetImageScaling(3.)

    add_watermark_root(c, pot)
    c.SaveAs(f"{outputDir}/png/{plot_name}.png")
    c.SaveAs(f"{outputDir}/pkl/{plot_name}.root")

def plotCovarianceMatrix(tree, outputDir, plot_name="covarianceMatrix", path="FitterEngine/postFit/Hesse/hessian/postfitCovariance_TH2D", pot=POT_DEFAULT, hist=None, labels=None):
    plotMatrix(tree, path, outputDir, plot_name, "Covariance", "kBird", pot=pot, hist=hist, labels=labels)

def plotCorrelationMatrix(tree, outputDir, plot_name="correlationMatrix", path="FitterEngine/postFit/Hesse/hessian/postfitCorrelation_TH2D", pot=POT_DEFAULT, hist=None, labels=None):
    plotMatrix(tree, path, outputDir, plot_name, "Correlation", "kThermometer", pot=pot, hist=hist, labels=labels)

def plotPulls(tree, path, outputDir, plot_name, pot=POT_DEFAULT, hists=None, labels=None):
    """Redraw GUNDAM's pre/post-fit parameter canvas on a fresh canvas: no title, rotated renamed labels.

    Pass hists=(prefit, postfit, line) and labels to draw already-built histograms instead of reading path.
    """

    r.gStyle.SetOptTitle(0)
    r.gStyle.SetOptStat(0)

    # Clone everything BEFORE drawing, and tell cppyy that ROOT (not
    # Python's refcounting) owns them, since they'll be attached to a
    # canvas below. Without this, Python's GC and the canvas's own
    # primitive-list cleanup can both try to free the same object,
    # causing "double free or corruption" at exit.
    def clone_owned(obj):
        clone = obj.Clone()
        r.SetOwnership(clone, False)
        # Same TFile-directory auto-registration issue as in plotMatrix:
        # detach histogram-like clones so the source file's Close() call
        # can't delete them while they're still drawn on our canvas.
        if hasattr(clone, "SetDirectory"):
            clone.SetDirectory(0)
        return clone

    if hists is not None:
        h_prefit, h_postfit, h_line = (clone_owned(h) if h else None for h in hists)
        for i, lbl in enumerate(labels, start=1):
            h_prefit.GetXaxis().SetBinLabel(i, lbl)
        labels = list(labels)
    else:
        c_par = tree.Get(path)
        r.SetOwnership(c_par, False)
        h_prefit = clone_owned(c_par.GetPrimitive("preFitErrors"))
        h_postfit = clone_owned(c_par.GetPrimitive("postFitErrors"))
        h_line = c_par.GetPrimitive("preFitErrorHistLine")
        h_line = clone_owned(h_line) if h_line else None

        template_labels = load_template_bin_labels()
        labels = relabel_axis(h_prefit.GetXaxis(), template_labels, 'Template Parameter' in path)

    n_bins = h_prefit.GetNbinsX()
    frame_w = max(1600, 60 * n_bins)
    frame_h = 1000
    px = label_px(frame_w, n_bins, cap_frac=0.035)
    left = max(260, rotated_label_overhang(labels, px, frame_w / n_bins) + 20)
    bottom = rotated_label_height(labels, px) + 40
    c = make_pixel_canvas(f"c_{plot_name}", frame_w, frame_h, left, 60, bottom, 130)
    c.SetGridx()
    c.SetGridy()

    h_prefit.SetTitle("")
    h_prefit.GetXaxis().SetTitle("")
    style_axis(h_prefit.GetXaxis(), px, rotate=True)
    # GUNDAM's marker/line sizes were set for a much smaller canvas.
    h_postfit.SetMarkerSize(3.5)
    h_postfit.SetLineWidth(10)
    plot_cap_px, legend_cap_px = 18, 6
    r.gStyle.SetEndErrorSize(plot_cap_px)
    # Template canvases fill the pre-fit band white; match the other pull plots instead.
    h_prefit.SetFillColor(623)
    h_prefit.SetLineColor(r.kRed + 1)
    h_prefit.SetLineWidth(5)
    if h_line:
        h_line.SetLineColor(r.kRed + 1)
        h_line.SetLineWidth(5)
    style_axis(h_prefit.GetYaxis(), 50, title_px=60, title_offset=1.6)

    h_prefit.Draw("E2")
    if h_line:
        h_line.Draw("SAME")
    h_postfit.Draw("E1 X0 SAME")
    c.RedrawAxis()

    legend = r.TLegend(0.72, 0.80, 0.98, 0.98 - 130 / c.GetWh())
    r.SetOwnership(legend, False)
    legend.SetTextFont(ROOT_FONT)
    legend.SetTextSize(50)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_prefit, "Pre-fit values", "fl")
    legend.AddEntry(h_postfit, "Post-fit values", "ep")
    # The cap size is a global style setting read at paint time, so shrink it only
    # while the legend is painted (TExecs run in primitive order).
    exec_small = r.TExec("exec_legend_caps", f"gStyle->SetEndErrorSize({legend_cap_px});")
    exec_reset = r.TExec("exec_reset_caps", f"gStyle->SetEndErrorSize({plot_cap_px});")
    for ex in (exec_small, exec_reset):
        r.SetOwnership(ex, False)
    exec_small.Draw()
    legend.Draw()
    exec_reset.Draw()

    c.Modified()
    c.Update()
    add_watermark_root(c, pot)
    c.SaveAs(f"{outputDir}/png/{plot_name}.png")
    c.SaveAs(f"{outputDir}/pkl/{plot_name}.root")


# Output-file prefixes for the parameter sets that had hand-written names
# before the per-set loop, so existing plot file names don't change. Any other
# set gets a camelCase prefix from its name (parameter_set_prefix).
PARAMETER_SET_PREFIXES = {
    "Cross section Systematics": "crossSection",
    "Flux Systematics Beamline": "flux",
    "Detector Systematics":      "detectorSystematics",
}
PULLS_SUFFIX = {"flux": "ParameterFitPulls"}  # historical name: fluxParameterFitPulls


def parameter_set_prefix(set_name):
    if set_name.startswith("Template Parameter"):
        return "templateParameter"
    if set_name in PARAMETER_SET_PREFIXES:
        return PARAMETER_SET_PREFIXES[set_name]
    words = [w for w in re.split(r"[^A-Za-z0-9]+", set_name) if w]
    # "G4 PiPlus Reinteraction" -> "g4PiPlusReinteraction"
    return words[0].lower() + "".join(w[0].upper() + w[1:] for w in words[1:])


def list_parameter_sets(fit_file):
    """Names of the parameter sets with post-fit Hesse errors in the fit file."""
    d = fit_file.Get("FitterEngine/postFit/Hesse/errors")
    if not d:
        return []
    return [k.GetName() for k in d.GetListOfKeys()]


def plot_g4_combined(fit_file, g4_sets, outputDir, pot=POT_DEFAULT):
    """One pulls plot and one covariance/correlation matrix for all G4 reinteraction sets together.

    The matrices are block-diagonal: each block is that set's own post-fit matrix.
    The fit stores the G4 sets only eigen-decomposed in the total covariance, so
    there is no physical-basis cross-set block to fill in.
    """
    labels, blocks, pulls = [], [], []
    for s in g4_sets:
        base = f"FitterEngine/postFit/Hesse/errors/{s}"
        cov = fit_file.Get(f"{base}/matrices/Covariance_TH2D")
        corr = fit_file.Get(f"{base}/matrices/Correlation_TH2D")
        c = fit_file.Get(f"{base}/values/fitConstraints_TCanvas")
        n = cov.GetNbinsX()
        lbls = g4_bin_labels(s) or [f"{s} #{i}" for i in range(n)]
        if len(lbls) != n:
            print(f"Warning: {s} has {n} parameters but {len(lbls)} binning labels; using indices")
            lbls = [f"{s} #{i}" for i in range(n)]
        labels += lbls
        blocks.append((cov, corr))
        pulls.append(tuple(c.GetPrimitive(k) for k in ("preFitErrors", "postFitErrors", "preFitErrorHistLine")))

    N = len(labels)
    h_cov = r.TH2D("g4_cov", "", N, 0, N, N, 0, N)
    h_corr = r.TH2D("g4_corr", "", N, 0, N, N, 0, N)
    off = 0
    for cov, corr in blocks:
        n = cov.GetNbinsX()
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                h_cov.SetBinContent(off + i, off + j, cov.GetBinContent(i, j))
                h_corr.SetBinContent(off + i, off + j, corr.GetBinContent(i, j))
        off += n
    for h in (h_cov, h_corr):
        h.SetDirectory(0)
        r.SetOwnership(h, False)
    h_corr.GetZaxis().SetRangeUser(-1, 1)

    def concat(name, parts):
        if any(p is None for p in parts):
            return None
        h = r.TH1D(name, "", N, 0, N)
        h.SetDirectory(0)
        r.SetOwnership(h, False)
        k = 1
        for p in parts:
            for i in range(1, p.GetNbinsX() + 1):
                h.SetBinContent(k, p.GetBinContent(i))
                h.SetBinError(k, p.GetBinError(i))
                k += 1
        # Keep the source style (colours, markers, fill).
        for attr in ("LineColor", "LineWidth", "LineStyle", "MarkerColor", "MarkerStyle",
                     "MarkerSize", "FillColor", "FillStyle"):
            getattr(h, f"Set{attr}")(getattr(parts[0], f"Get{attr}")())
        return h

    h_pre = concat("g4_pre", [p[0] for p in pulls])
    h_post = concat("g4_post", [p[1] for p in pulls])
    h_line = concat("g4_line", [p[2] for p in pulls])
    if h_line:
        # The prior-value line is drawn as a plain line: no error bars.
        for i in range(1, N + 1):
            h_line.SetBinError(i, 0)
    # GUNDAM sets the y range per set; use one that covers every set.
    lo = min(h_pre.GetBinContent(i) - h_pre.GetBinError(i) for i in range(1, N + 1))
    hi = max(h_pre.GetBinContent(i) + h_pre.GetBinError(i) for i in range(1, N + 1))
    lo = min(lo, min(h_post.GetBinContent(i) - h_post.GetBinError(i) for i in range(1, N + 1)))
    hi = max(hi, max(h_post.GetBinContent(i) + h_post.GetBinError(i) for i in range(1, N + 1)))
    pad = 0.15 * (hi - lo)
    h_pre.GetYaxis().SetRangeUser(lo - pad, hi + pad)
    h_pre.GetYaxis().SetTitle(pulls[0][0].GetYaxis().GetTitle() or "Parameter value")

    plotCovarianceMatrix(None, outputDir, plot_name="g4CovarianceMatrix", pot=pot, hist=h_cov, labels=labels)
    plotCorrelationMatrix(None, outputDir, plot_name="g4CorrelationMatrix", pot=pot, hist=h_corr, labels=labels)
    plotPulls(None, None, outputDir, "g4FitPulls", pot=pot, hists=(h_pre, h_post, h_line), labels=labels)


def draw_stacked_with_band(ax, stack_data, errors, title, x_vals, widths, labels, category_colors, data_y=None, data_yerr=None):
    n_bins_local = len(x_vals)
    bottom = np.zeros(n_bins_local)
    totals_per_cat = stack_data.sum(axis=1)
    order = np.argsort(totals_per_cat)  # draw from smallest to largest contribution

    for i in order:
        ax.bar(
            x_vals,
            stack_data[i],
            bottom=bottom,
            width=widths,
            label=labels[i],
            color=category_colors[i],
            align='center'
        )
        bottom += stack_data[i]

    # total prediction and band
    x_edges = [x_vals[0] - widths[0] / 2.0]
    x_edges.extend([x + w / 2.0 for x, w in zip(x_vals, widths)])
    total_pred = np.sum(stack_data, axis=0)
    err = np.asarray(errors, dtype=float)
    lower = np.maximum(0.0, total_pred - err)
    upper = total_pred + err
    ax.fill_between(
        x_edges,
        np.r_[lower, lower[-1]],
        np.r_[upper, upper[-1]],
        step='post',
        facecolor='0.8',
        edgecolor='gray',
        hatch='///',
        alpha=0.5,
        label='Uncertainty'
    )

    # data points
    if data_y is not None:
        y = np.asarray(data_y, dtype=float)
        yerr = np.asarray(data_yerr, dtype=float) if data_yerr is not None else np.zeros_like(y)
        ax.errorbar(
            x_vals,
            y,
            yerr=yerr,
            fmt='o',
            color='black',
            ecolor='black',
            elinewidth=1.5,
            capsize=0,
            markersize=4,
            label='Data',
            zorder=10
        )

    ax.set_ylabel("Counts")
    ax.text(0.02, 0.97, title, transform=ax.transAxes, ha='left', va='top', fontsize=14, fontweight='bold')
    return x_edges, total_pred, err

def draw_ratio(ax, total_pred, err, data_y, data_yerr, x_edges, x_vals, center=1.0, ypad=(0.5, 1.5), xlabel=None):
    tot = np.asarray(total_pred, dtype=float)
    data_y = np.asarray(data_y, dtype=float)
    data_yerr = np.asarray(data_yerr, dtype=float)

    denom = tot.copy()
    denom[denom == 0.0] = np.nan
    ratio = data_y / denom
    ratio_err = data_yerr / denom

    lower = np.maximum(0.0, (tot - err) / denom)
    upper = (tot + err) / denom

    ax.fill_between(
        x_edges,
        np.r_[lower, lower[-1]],
        np.r_[upper, upper[-1]],
        step='post',
        facecolor='0.8',
        edgecolor='none',
        alpha=0.6,
        zorder=1,
    )
    ax.errorbar(x_vals, ratio, yerr=ratio_err, fmt='o', color='black', zorder=10)
    ax.axhline(center, color='gray', linestyle='--', linewidth=1)
    ax.set_ylim(ypad)
    ax.set_xlabel(xlabel or analysis_vars.current()['reco_axis_title'])
    ax.set_ylabel("Data / Pred")

def get_mc_hist(file_obj, category_index, mode, sample=None, branch=None):
    """Get category histogram from ROOT file with path fallbacks across file variants."""
    v = analysis_vars.current()
    s, b = sample or v['reco_sample'], branch or v['reco_branch']
    candidates = []
    if mode == "pre":
        candidates = [
            f"toyGen/plots/histograms/{s} (pre-fit)/{b}/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/{s}/{b}/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/{s}_pre_fit/{b}/category/{category_index}/MC_TH1D",
        ]
    else:
        candidates = [
            f"toyGen/plots/histograms/{s} (post-fit)/{b}/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/{s}/{b}/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/{s}_post_fit/{b}/category/{category_index}/MC_TH1D",
        ]

    for path in candidates:
        obj = file_obj.Get(path)
        if obj and hasattr(obj, "GetBinContent"):
            return obj
    print("Warning no MC_TH1D found for category", category_index, "in", mode, "fit file. Tried paths:", candidates)
    return None

def getHistErrors(fPre, fPost, sample=None, branch=None):
    """Pre/post-fit errors, per-category MC and data for one sample (default: the signal sample)."""
    v = analysis_vars.current()
    s, b = sample or v['reco_sample'], branch or v['reco_branch']
    h_errors_pre = fPre.Get(f"toyGen/histograms/{s}_pre_fit_TH1D")  # errors
    h_errors_post = fPost.Get(f"toyGen/histograms/{s}_TH1D")  # errors

    h_data = fPost.Get(f"toyGen/plots/histograms/{s}/{b}/Data_TH1D")  # data

    postFitErrors = []
    preFitErrors = []
    x_vals = []
    widths = []
    for i in range(1, h_errors_pre.GetNbinsX() + 1):
        bounds = [float(x) for x in h_errors_pre.GetXaxis().GetBinLabel(i).split(':')[-1].strip().replace("[", "").replace("]", "").split(',')]
        width = bounds[1] - bounds[0]
        widths.append(width)
        preFitErrors.append(h_errors_pre.GetBinError(i) * width)
        postFitErrors.append(h_errors_post.GetBinError(i) * width)
        x_vals.append((bounds[0] + bounds[1]) / 2)

    # extract MC per-category histograms (fallback paths handled by get_mc_hist)
    # categories 5 and 10 are both cosmics and are summed into one "Cosmics" entry
    h_mc_pre = []
    h_mc_post = []
    for cat_group in [[0], [1], [2], [3], [4], [5, 10]]:
        n_cat_bins = len(x_vals)
        pre_bin_content = [0.0] * n_cat_bins
        post_bin_content = [0.0] * n_cat_bins
        for cat in cat_group:
            h_data_pre = get_mc_hist(fPre, cat, "pre", s, b)
            h_data_post = get_mc_hist(fPost, cat, "post", s, b)
            if h_data_pre is not None:
                n_cat_bins = h_data_pre.GetNbinsX()
            elif h_data_post is not None:
                n_cat_bins = h_data_post.GetNbinsX()
            for j in range(n_cat_bins):
                pre_bin_content[j] += h_data_pre.GetBinContent(j + 1) if h_data_pre is not None else 0.0
                post_bin_content[j] += h_data_post.GetBinContent(j + 1) if h_data_post is not None else 0.0
        h_mc_pre.append(pre_bin_content)
        h_mc_post.append(post_bin_content)

    data = []
    data_err = []
    for i in range(1, h_data.GetNbinsX() + 1):
        data.append(h_data.GetBinContent(i))
        data_err.append(h_data.GetBinError(i))

    # labels and colors matching notebook
    labels = [r"$\nu_e$ Signal", r"$\nu_e$ OOFV", r"$\nu_\mu$", r"$\nu$ OOFV", r"$\nu$ NC", 'Cosmics']
    category_colors = ['#011959', '#1C5A62', '#687B3E', '#D29343', '#FDB7BC', '#ca80a0ff']

    n_bins = len(x_vals)

    return x_vals, widths, preFitErrors, postFitErrors, np.array(h_mc_pre), np.array(h_mc_post), data, data_err, n_bins, labels, category_colors

def plot_pre_post_comparison(prefit_file, postfit_file, outputDir, pot=POT_DEFAULT,
                             sample=None, branch=None, xlabel=None, plot_name="pre_post_fit_comparison"):
    """Create side-by-side pre/post stacked plots with ratio panels and save image.

    Defaults to the signal sample; pass sample/branch/xlabel/plot_name for another
    sample (e.g. the NC pi0 sideband).
    """
    try:
        res = getHistErrors(prefit_file, postfit_file, sample, branch)
        x_vals, widths, preFitErrors, postFitErrors, h_mc_pre, h_mc_post, data, data_err, n_bins, labels, category_colors = res

        # set globals expected by existing helper functions
        globals()['x_vals'] = x_vals
        globals()['widths'] = widths
        globals()['n_bins'] = n_bins
        globals()['labels'] = labels
        globals()['category_colors'] = category_colors

        data_pre = np.array(h_mc_pre)
        data_post = np.array(h_mc_post)

        fig, axes = plt.subplots(2, 2, figsize=(20, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex='col')
        ax1_main, ax2_main = axes[0, 0], axes[0, 1]
        ax1_ratio, ax2_ratio = axes[1, 0], axes[1, 1]

        edges1, tot1, err1 = draw_stacked_with_band(ax1_main, data_pre, preFitErrors, "Pre-fit", x_vals, widths, labels, category_colors, data, data_err)
        edges2, tot2, err2 = draw_stacked_with_band(ax2_main, data_post, postFitErrors, "Post-fit", x_vals, widths, labels, category_colors, data, data_err)

        draw_ratio(ax1_ratio, tot1, err1, data, data_err, edges1, x_vals, xlabel=xlabel)
        draw_ratio(ax2_ratio, tot2, err2, data, data_err, edges2, x_vals, xlabel=xlabel)

        # shared y-limits
        ymax1 = np.nanmax(tot1 + err1) if len(tot1) and len(err1) else 0.0
        ymax2 = np.nanmax(tot2 + err2) if len(tot2) and len(err2) else 0.0
        shared_ymax = max(ymax1, ymax2)
        pad = 0.25
        shared_ymax = max(1.0, shared_ymax * (1.0 + pad))
        ax1_main.set_ylim(0.0, shared_ymax)
        ax2_main.set_ylim(0.0, shared_ymax)
        # Put the legend over whichever side of the plot has the lower prediction
        # (EleE peaks at low energy, EleAngle at forward cos theta).
        top = np.maximum(tot1 + err1, tot2 + err2)
        third = max(1, len(top) // 3)
        if np.nanmax(top[:third]) > np.nanmax(top[-third:]):
            legend_kw = dict(loc='upper right', ncol=2)
        else:
            legend_kw = dict(loc='upper left', bbox_to_anchor=(0.0, 0.92), ncol=2)
        for ax in (ax1_main, ax2_main):
            ax.legend(**legend_kw)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        stamp_ax(ax1_main, pot)
        stamp_ax(ax2_main, pot)

        fig.savefig(f"{outputDir}/png/{plot_name}.png")
        with open(f"{outputDir}/pkl/{plot_name}.pkl", "wb") as f:
            pickle.dump(fig, f)
        plt.close(fig)
    except Exception as e:
        print(f"Warning: failed to make pre/post fit comparison ({plot_name}):", e)

def compute_flux_integral(nuisflat_file, bin_width_divided=False, is_gibuu=False):
    """Integrate a nuisflat file's FlatTree_FLUX into a flux/POT/cm^2 scalar."""
    with uproot.open(nuisflat_file) as f:
        counts, bin_edges = f['FlatTree_FLUX'].to_numpy()
        integral = np.sum(counts * np.diff(bin_edges)) if bin_width_divided else np.sum(counts)
    # GiBUU flux is already given in /cm^2; everything else needs m^2->cm^2 and MeV->GeV
    return integral if is_gibuu else integral * (1e-4 / 9.99e8)


def load_nuisflat_dataframe(file_path, branches, sf_var_name='FinalWeight', sf_var_expr='InputWeight * fScaleFactor'):
    """Load a NUISANCE-flattened GENIE tree, adding the derived cross-section weight column."""
    with uproot.open(file_path) as f:
        df = f['FlatTree_VARS'].arrays(branches, library='pd')
    return df.eval(f'{sf_var_name} = {sf_var_expr}')


def compute_genie_xsec(df, variable_name, edges, signal_selection_expr, weight_col='FinalWeight',
                        n_files=1, per_nucleon=False, n_nucleon=40):
    """Differential GENIE cross section straight from a NUISANCE flat tree's pre-normalized weights."""
    df_sel = df.query(signal_selection_expr)
    bin_widths = np.diff(edges)

    counts = np.histogram(df_sel[variable_name], bins=edges, weights=df_sel[weight_col])[0]
    counts_err2 = np.histogram(df_sel[variable_name], bins=edges, weights=df_sel[weight_col] ** 2)[0]

    # NUISANCE's fScaleFactor already yields a per-nucleon cross section;
    # scale up to per-target (Ar) unless the caller wants per-nucleon.
    scale = 1.0 if per_nucleon else n_nucleon
    xsec = counts / bin_widths / n_files * scale
    xsec_err = np.sqrt(counts_err2) / bin_widths / n_files * scale

    return xsec, xsec_err


def combine_flux_averaged_xsec(xsecs, xsec_errs, flux_integrals):
    """Combine per-species cross sections into a flux-weighted average."""
    flux_sum = sum(flux_integrals)
    xsec_sum = sum(x * f for x, f in zip(xsecs, flux_integrals))
    xsec_err2_sum = sum((e * f) ** 2 for e, f in zip(xsec_errs, flux_integrals))
    return xsec_sum / flux_sum, np.sqrt(xsec_err2_sum) / flux_sum


_REWEIGHT_NAMESPACE = {
    "where": np.where, "sin": np.sin, "cos": np.cos, "tan": np.tan,
    "arcsin": np.arcsin, "arccos": np.arccos, "arctan": np.arctan,
    "asin": np.arcsin, "acos": np.arccos, "atan": np.arctan,
    "sqrt": np.sqrt, "abs": np.abs, "exp": np.exp, "log": np.log,
    "minimum": np.minimum, "maximum": np.maximum, "clip": np.clip,
    "power": np.power, "pow": np.power, "pi": np.pi, "e": np.e,
}


def _eval_tree_expr(expr, tree):
    """Evaluate a numpy expression over a TTree's branches.

    Any identifier in `expr` that names a branch is loaded and bound; the rest
    come from a small whitelist of numpy functions/constants. This lets the
    fake-data reweight strings mirror the YAML weight formulas (with `where`
    standing in for the C-style ternary). No Python builtins are exposed.
    """
    identifiers = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr))
    available = set(tree.keys())
    needed = sorted(identifiers & available)
    namespace = dict(_REWEIGHT_NAMESPACE)
    if needed:
        namespace.update(tree.arrays(needed, library="np"))
    unknown = identifiers - available - set(_REWEIGHT_NAMESPACE)
    if unknown:
        raise ValueError(
            f"reweight expression references unknown names {sorted(unknown)}; "
            f"only truth-tree branches and {sorted(_REWEIGHT_NAMESPACE)} are allowed"
        )
    return eval(expr, {"__builtins__": {}}, namespace)  # noqa: S307 (whitelisted namespace)


def compute_mc_truth_xsec(edges, reweight_expr="1.0",
                          flux=MC_TRUTH_FLUX, ntargets=MC_TRUTH_NTARGETS):
    """dsigma/dE from the CAF signal-truth tree, optionally scaled by a reweight.

    Reproduces gundamCalcXsec's signal truth (same tree, weight, selection and
    flux x target normalisation), times `reweight_expr`. With the default
    reweight this is the nominal truth; with a fake-data reweight it is the true
    cross section of the injected fake data — the correct closure reference.
    Returns (xsec, xsec_stat_err) per unit of the plotted variable per Ar on
    `edges` (GeV for EleE, cos theta for EleAngle).
    """
    v = analysis_vars.current()
    with uproot.open(MC_TRUTH_FILE) as f:
        tree = f[MC_TRUTH_TREE]
        var = tree[v['truth_branch']].array(library="np") * v['unit_scale']  # e.g. MeV -> GeV
        nominal_w = np.asarray(_eval_tree_expr(MC_TRUTH_NOMINAL_WEIGHT, tree), dtype=float)
        selection = np.asarray(_eval_tree_expr(MC_TRUTH_SELECTION, tree)).astype(bool)
        reweight = np.broadcast_to(
            np.asarray(_eval_tree_expr(reweight_expr, tree), dtype=float), var.shape)

    weight = nominal_w * reweight
    bin_widths = np.diff(edges)
    counts = np.histogram(var[selection], bins=edges, weights=weight[selection])[0]
    counts_err2 = np.histogram(var[selection], bins=edges, weights=weight[selection] ** 2)[0]

    norm = bin_widths * flux * ntargets
    return counts / norm, np.sqrt(counts_err2) / norm


def parse_gundam_bin_edges(axis_labels, unit_scale=1e-3):
    """Parse physical bin edges from GUNDAM '[lo, hi[' axis labels.

    GUNDAM stores TH1 bin edges as sequential integers; the physical ranges live
    in the bin labels, e.g. 'true_leading_electron_energy: [300, 450['. The
    default unit_scale converts the stored MeV to GeV.
    """
    edges = []
    for i, label in enumerate(axis_labels):
        nums = label.split(":")[-1].strip().lstrip("[").split(",")
        lo = float(nums[0].strip())
        hi = float(nums[1].strip().rstrip("[]()"))
        if i == 0:
            edges.append(lo * unit_scale)
        edges.append(hi * unit_scale)
    return np.array(edges)


def load_gundam_xsec(root_file, hist_name, unit_scale=1.0):
    """Read a GUNDAM cross-section TH1D, its stat errors, and physical bin edges (GeV).

    Edges come from the bin labels (not to_numpy's integer indices), so callers
    adapt to any (re)binning of truth space.
    """
    with uproot.open(root_file) as f:
        h_xsec = f[hist_name]
        counts = h_xsec.values()
        errs = h_xsec.errors()
        labels = h_xsec.axis(0).labels()
    edges = parse_gundam_bin_edges(labels, unit_scale=analysis_vars.current()['unit_scale'])
    return counts * unit_scale, errs * unit_scale, edges


def load_gundam_xsec_cov(root_file, cov_name="calcXsec/matrices/covarianceMatrix_TH2D",
                         unit_scale=1.0):
    """Read the calcXsec xsec covariance matrix, scaled to match load_gundam_xsec values.

    The matrix has units of (xsec)^2, so the value unit_scale enters squared. Returns
    None if the matrix is absent (older outputs) so callers can fall back to the diagonal.
    """
    try:
        with uproot.open(root_file) as f:
            cov = f[cov_name].values()
        return np.asarray(cov, dtype=float) * (unit_scale ** 2)
    except Exception as e:
        print("Warning: xsec covariance unavailable, using diagonal errors:", e)
        return None


def _cov_chi2(residual, cov):
    """Correlated chi2 = r^T C^-1 r, with ndof and p-value. Falls back to a
    pseudo-inverse if C is singular. Returns (chi2, ndof, pvalue)."""
    residual = np.asarray(residual, dtype=float)
    ndof = residual.size
    try:
        solved = np.linalg.solve(cov, residual)
    except np.linalg.LinAlgError:
        solved = np.linalg.pinv(cov) @ residual
    chi2_val = float(residual @ solved)
    return chi2_val, ndof, float(stats.chi2.sf(chi2_val, ndof))


def _chi2_pvalue(vals_a, err_a, vals_b, err_b):
    """Uncorrelated-bin chi2 and p-value between two binned distributions with independent errors."""
    vals_a = np.asarray(vals_a, dtype=float)
    vals_b = np.asarray(vals_b, dtype=float)
    sigma2 = np.asarray(err_a, dtype=float) ** 2 + np.asarray(err_b, dtype=float) ** 2
    mask = np.isfinite(vals_a) & np.isfinite(vals_b) & np.isfinite(sigma2) & (sigma2 > 0)
    ndof = int(np.sum(mask))
    if ndof == 0:
        return np.nan, 0, np.nan
    chi2_val = float(np.sum((vals_b[mask] - vals_a[mask]) ** 2 / sigma2[mask]))
    pvalue = float(stats.chi2.sf(chi2_val, ndof))
    return chi2_val, ndof, pvalue


def plot_xsec_comparison(edges, ref_xsec, ref_xsec_err, gundam_xsec, gundam_xsec_err,
                          xlabel, ylabel, ref_label='GENIE', ylim=None,
                          ratio_ylim=None, outpath=None, pot=POT_DEFAULT, fontsize=18,
                          extra_xsec=None, extra_xsec_err=None, extra_label=None,
                          chi2_text=None):
    """Draw a two-panel cross-section comparison: a reference (the closure target) vs GUNDAM, with a ratio panel.

    `ref_xsec` is the primary reference the chi2 and ratio are computed against.
    An optional `extra_xsec` (e.g. an independent GENIE generator prediction) is
    overlaid for context only and does not enter the chi2.

    Any '{scale}' in ylabel is replaced by the common power of ten, which is divided out of
    the plotted values so matplotlib's offset text doesn't collide with the watermark.
    """
    bin_widths = np.diff(edges)
    bin_centers = (edges[:-1] + edges[1:]) / 2

    exponent = int(np.floor(np.log10(np.nanmax(np.abs(ref_xsec)))))
    unit = 10.0 ** exponent
    ref_xsec, ref_xsec_err = ref_xsec / unit, ref_xsec_err / unit
    gundam_xsec, gundam_xsec_err = gundam_xsec / unit, gundam_xsec_err / unit
    if extra_xsec is not None:
        extra_xsec = np.asarray(extra_xsec) / unit
    ylabel = ylabel.replace('{scale}', rf'$10^{{{exponent}}}$')
    if ylim:
        ylim = (ylim[0] / unit, ylim[1] / unit)

    fig, (ax_xsec, ax_ratio) = plt.subplots(
        2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]}, sharex=True,
    )
    fig.subplots_adjust(hspace=0.05)

    ax_xsec.hist(edges[:-1], bins=edges, weights=ref_xsec, label=ref_label, histtype='step', alpha=1.0, color='blue', linewidth=2)
    if extra_xsec is not None:
        ax_xsec.hist(edges[:-1], bins=edges, weights=extra_xsec,
                      label=extra_label or 'reference', histtype='step',
                      alpha=0.7, color='green', linewidth=2, linestyle='--')
    ax_xsec.errorbar(bin_centers, gundam_xsec, yerr=gundam_xsec_err, xerr=bin_widths / 2,
                      fmt='o', capsize=5, label='GUNDAM', color='black', markersize=6)
    # Legend and chi2 text go on whichever side of the plot is lower
    # (EleE peaks at low energy, EleAngle at forward cos theta).
    span = edges[-1] - edges[0]
    tops = np.nan_to_num(np.maximum(ref_xsec, gundam_xsec + gundam_xsec_err))
    left_top = np.max(tops[bin_centers < edges[0] + 0.4 * span], initial=0.0)
    right_top = np.max(tops[bin_centers > edges[-1] - 0.4 * span], initial=0.0)
    text_on_left = right_top > left_top
    text_x, text_ha = (0.03, 'left') if text_on_left else (0.97, 'right')
    ax_xsec.legend(loc='upper left' if text_on_left else 'upper right', fontsize=fontsize)
    ax_xsec.set_ylabel(ylabel)
    if ylim:
        ax_xsec.set_ylim(*ylim)
    else:
        top = np.nanmax([np.nanmax(ref_xsec), np.nanmax(gundam_xsec + gundam_xsec_err)])
        ax_xsec.set_ylim(0, top * 1.15)

    if chi2_text is None:
        chi2_val, ndof, pvalue = _chi2_pvalue(ref_xsec, ref_xsec_err, gundam_xsec, gundam_xsec_err)
        chi2_ndof = chi2_val / ndof if ndof > 0 else np.nan
        chi2_text = (r'$\chi^2$/ndof = %.1f/%d = %.2f' % (chi2_val, ndof, chi2_ndof)
                     + '\n' + r'p = %.3f' % (pvalue))
    ax_xsec.text(
        text_x, 0.74, chi2_text,
        transform=ax_xsec.transAxes,
        fontsize=fontsize, ha=text_ha, va='top',
    )

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(ref_xsec > 0, gundam_xsec / ref_xsec, np.nan)
        ratio_err = np.where(
            ref_xsec > 0,
            ratio * np.sqrt(
                (gundam_xsec_err / np.where(gundam_xsec > 0, gundam_xsec, np.nan)) ** 2
                + (ref_xsec_err / ref_xsec) ** 2
            ),
            np.nan,
        )

    ax_ratio.errorbar(bin_centers, ratio, yerr=ratio_err, xerr=bin_widths / 2, fmt='o', capsize=5, color='k')
    ax_ratio.axhline(1, color='black', linestyle='--', linewidth=1)
    ax_ratio.set_ylabel(f'GUNDAM / {ref_label}')
    ax_ratio.set_xlabel(xlabel)
    if ratio_ylim:
        ax_ratio.set_ylim(*ratio_ylim)
    else:
        finite = np.isfinite(ratio) & np.isfinite(ratio_err)
        max_dev = np.nanmax(np.abs(ratio[finite] - 1) + ratio_err[finite]) if np.any(finite) else 0.3
        pad = max(max_dev * 1.2, 0.1)
        ax_ratio.set_ylim(1 - pad, 1 + pad)

    stamp_ax(ax_xsec, pot)

    if outpath:
        fig.savefig(outpath, bbox_inches='tight')

    return fig, (ax_xsec, ax_ratio)


def _external_genie_xsec(xsec_edges):
    """Flux-averaged GENIE prediction from the independent NUISANCE flat trees.

    This is a separate generator sample (not the CAF MC the fit runs on), so it
    serves only as a context overlay, never as the closure target. Returns
    (xsec, xsec_err) or (None, None) if the flat trees are unavailable.
    """
    nue_flat_trees = {
        'nue': '/exp/icarus/data/users/micarrig/genie/nue_samples_v3_04_00/p12_1000180400_CC_v3_4_0_AR23_20i_00_000.flat.root',
        'nuebar': '/exp/icarus/data/users/micarrig/genie/nue_samples_v3_04_00/p-12_1000180400_CC_v3_4_0_AR23_20i_00_000.flat.root',
    }
    genie_var = analysis_vars.current()['genie_branch']
    nue_branches = ['ICARUS_1eInc_IsSignal', genie_var, 'InputWeight', 'fScaleFactor']
    try:
        genie_xsecs, genie_xsec_errs, flux_integrals = [], [], []
        for file_path in nue_flat_trees.values():
            df = load_nuisflat_dataframe(file_path, nue_branches)
            xsec, xsec_err = compute_genie_xsec(
                df, genie_var, xsec_edges,
                signal_selection_expr='(ICARUS_1eInc_IsSignal==True)',
            )
            genie_xsecs.append(xsec)
            genie_xsec_errs.append(xsec_err)
            flux_integrals.append(compute_flux_integral(file_path))
        return combine_flux_averaged_xsec(genie_xsecs, genie_xsec_errs, flux_integrals)
    except Exception as e:
        print("Warning: external GENIE overlay unavailable:", e)
        return None, None


def closure_chi2(gundam_xsec_file, truth_reweight="1.0", xsec=None, truth=None):
    """Closure chi2 of the GUNDAM xsec against the (reweighted) CAF truth.

    The correlated chi2 uses the full calcXsec covariance (xsec bins are strongly
    correlated) plus the truth's MC-stat error on the diagonal, which is negligible
    but keeps the two estimates independent. Uses the variable selected in
    analysis_vars. `xsec`/`truth` let a caller that has already loaded them skip
    the reload. Returns a dict: chi2_corr, p_corr (None if no covariance),
    chi2_uncorr, p_uncorr, ndof.
    """
    v = analysis_vars.current()
    xsec_scale = 1.0 / v['unit_scale']
    if xsec is None:
        xsec = load_gundam_xsec(gundam_xsec_file, v['xsec_hist'], unit_scale=xsec_scale)
    gundam_xsec, gundam_xsec_err, edges = xsec
    if truth is None:
        truth = compute_mc_truth_xsec(edges, truth_reweight)
    truth_xsec, truth_xsec_err = truth

    chi2_d, ndof_d, p_d = _chi2_pvalue(truth_xsec, truth_xsec_err, gundam_xsec, gundam_xsec_err)
    out = dict(chi2_corr=None, p_corr=None, chi2_uncorr=chi2_d, p_uncorr=p_d, ndof=ndof_d)
    cov = load_gundam_xsec_cov(gundam_xsec_file, unit_scale=xsec_scale)
    if cov is not None:
        cov = cov + np.diag(np.asarray(truth_xsec_err, dtype=float) ** 2)
        chi2_c, ndof_c, p_c = _cov_chi2(gundam_xsec - truth_xsec, cov)
        out.update(chi2_corr=chi2_c, p_corr=p_c, ndof=ndof_c)
    return out


def plot_genie_gundam_comparison(outputDir, gundam_xsec_file, pot=POT_DEFAULT,
                                 truth_reweight="1.0", truth_label=None):
    """GUNDAM extracted vs the fake-data truth cross section (the closure test).

    The primary reference is the CAF signal truth scaled by `truth_reweight` —
    i.e. the true cross section of the injected fake data, which is what GUNDAM
    must reproduce for the study to close. The independent GENIE flat-tree
    prediction is overlaid for context only. `truth_reweight` is a numpy
    expression over the truth-tree branches (default "1.0" = nominal truth).
    """
    try:
        v = analysis_vars.current()
        gundam_xsec_hist = v['xsec_hist']
        # GUNDAM's xsec is per branch unit (per MeV for EleE); convert to per plot unit.
        xsec_scale = 1.0 / v['unit_scale']

        # Analysis binning is read from the GUNDAM xsec histogram's bin labels, so
        # the truth comparison adapts to any (re)binning of truth space.
        gundam_xsec, gundam_xsec_err, xsec_edges = load_gundam_xsec(
            gundam_xsec_file, gundam_xsec_hist, unit_scale=xsec_scale)

        # Closure reference: the fake data's own true cross section.
        truth_xsec, truth_xsec_err = compute_mc_truth_xsec(xsec_edges, truth_reweight)
        is_nominal = str(truth_reweight).strip() in ("1.0", "1", "")
        ref_label = truth_label or ("GENIE (nominal)" if is_nominal else "GENIE (reweighted)")

        res = closure_chi2(gundam_xsec_file, truth_reweight,
                           xsec=(gundam_xsec, gundam_xsec_err, xsec_edges),
                           truth=(truth_xsec, truth_xsec_err))
        if res['chi2_corr'] is not None:
            chi2_text = (
                r'$\chi^2$/ndof = %.1f/%d = %.2f (corr.)'
                % (res['chi2_corr'], res['ndof'], res['chi2_corr'] / res['ndof'])
                + '\n' + r'p = %.3f' % res['p_corr']
                + '\n' + r'$\chi^2_{\mathrm{uncorr}}$ = %.1f/%d' % (res['chi2_uncorr'], res['ndof'])
            )
        else:
            chi2_text = None  # plot_xsec_comparison falls back to the diagonal chi2

        # Independent GENIE generator prediction, overlaid for context only.
        # Suppress it when it would just duplicate a nominal primary curve.
        genie_xsec, _ = _external_genie_xsec(xsec_edges)
        extra_label = 'GENIE (nominal)'
        if genie_xsec is None or ref_label == extra_label:
            genie_xsec, extra_label = None, None

        fig, _ = plot_xsec_comparison(
            xsec_edges, truth_xsec, truth_xsec_err, gundam_xsec, gundam_xsec_err,
            xlabel=v['label'],
            ylabel=v['xsec_ylabel'],
            ref_label=ref_label,
            extra_xsec=genie_xsec,
            extra_label=extra_label,
            chi2_text=chi2_text,
            outpath=f"{outputDir}/png/genie_gundam_xsec_comparison.png",
            pot=pot,
        )
        with open(f"{outputDir}/pkl/genie_gundam_xsec_comparison.pkl", "wb") as f:
            pickle.dump(fig, f)
        plt.close(fig)
    except Exception as e:
        print("Warning: failed to make truth/GUNDAM cross-section comparison:", e)


def main(tag='asimov', pot=POT_DEFAULT, data_dir=None, output_dir=None,
         truth_reweight=None, truth_label=None):

    r.gROOT.SetBatch(True)
    plt.style.use(_style.get_style())

    # Directories default to the tag-based layout, but can be overridden (e.g.
    # when run_fake_study.sh is given a custom output directory).
    _base = '/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue'
    # Each variable has its own tree: data/<var>/<tag>, as run_fake_study.sh does.
    var_name = analysis_vars.current_name()
    data_dir = data_dir or f'{_base}/data/{var_name}/{tag}'
    outputDir = output_dir or f'{_base}/plots/{var_name}/{tag}'

    fitFile = f'{data_dir}/fit_{tag}.root'
    # xsecFile = f'{data_dir}/nue_asimovXSec.root'
    prefitFile = f'{data_dir}/toy_pf_{tag}.root'
    postfitFile = f'{data_dir}/toy_bf_{tag}.root'
    xsecFile = f'{data_dir}/xsec_{tag}.root'

    fit_file = r.TFile(fitFile)
    # xsec_file = r.TFile(xsecFile)
    prefit_file = r.TFile(prefitFile)
    postfit_file = r.TFile(postfitFile)
    xsec_file = r.TFile(xsecFile)

    os.makedirs(f"{outputDir}/png", exist_ok=True)
    os.makedirs(f"{outputDir}/pkl", exist_ok=True)

    # Total
    plotCovarianceMatrix(fit_file, outputDir, pot=pot)
    plotCorrelationMatrix(fit_file, outputDir, pot=pot)

    # Per-parameter-set pulls + covariance + correlation. Loops over every set
    # the fit stored, so new sets are picked up automatically. The G4
    # reinteraction sets are drawn together in one set of plots instead.
    all_sets = list_parameter_sets(fit_file)
    g4_sets = [s for s in all_sets if s.startswith("G4 ")]
    for set_name in all_sets:
        if set_name in g4_sets:
            continue
        prefix = parameter_set_prefix(set_name)
        base = f"FitterEngine/postFit/Hesse/errors/{set_name}"
        names = {"pulls": f"{prefix}{PULLS_SUFFIX.get(prefix, 'FitPulls')}",
                 "cov":   f"{prefix}CovarianceMatrix",
                 "corr":  f"{prefix}CorrelationMatrix"}
        try:
            plotCovarianceMatrix(fit_file, outputDir, plot_name=names["cov"], path=f"{base}/matrices/Covariance_TH2D", pot=pot)
            plotCorrelationMatrix(fit_file, outputDir, plot_name=names["corr"], path=f"{base}/matrices/Correlation_TH2D", pot=pot)
            plotPulls(fit_file, f"{base}/values/fitConstraints_TCanvas", outputDir, names["pulls"], pot=pot)
        except Exception as e:
            print(f"Warning: failed to make plots for parameter set '{set_name}':", e)
    if g4_sets:
        try:
            plot_g4_combined(fit_file, g4_sets, outputDir, pot=pot)
        except Exception as e:
            print("Warning: failed to make the combined G4 plots:", e)

    # Pre/Post fit stacked comparison (moved to helper)
    plot_pre_post_comparison(prefit_file, postfit_file, outputDir, pot=pot)

    # Same comparison for the NC pi0 sideband, when the fit has one.
    v = analysis_vars.current()
    sb = v.get('sideband_sample')
    if sb and postfit_file.Get(f"toyGen/histograms/{sb}_TH1D"):
        plot_pre_post_comparison(prefit_file, postfit_file, outputDir, pot=pot,
                                 sample=sb, branch=v['sideband_branch'],
                                 xlabel=v['sideband_axis_title'],
                                 plot_name="sideband_pre_post_fit_comparison")

    # GUNDAM extracted vs the fake-data truth cross section (closure test).
    # The reweight defaults to this tag's known injection; override with
    # --truth-reweight to feed a specific version.
    if truth_reweight is None:
        truth_reweight = TRUTH_REWEIGHTS.get(tag, "1.0")
        if tag not in TRUTH_REWEIGHTS:
            print(f"Note: no truth reweight registered for tag '{tag}'; "
                  f"using nominal truth. Pass --truth-reweight to override.")
    if truth_label is None and tag in TRUTH_REWEIGHT_LABELS:
        truth_label = f"GENIE ({TRUTH_REWEIGHT_LABELS[tag]})"
    plot_genie_gundam_comparison(outputDir, xsecFile, pot=pot,
                                 truth_reweight=truth_reweight, truth_label=truth_label)

    # Deliberately skip TFile.Close()/normal interpreter teardown: all PNGs
    # are already on disk by this point, the files were opened read-only, and
    # ROOT/cppyy's own exit-time cleanup of the many canvases and cloned
    # primitives created above triggers a "double free or corruption" abort
    # that has nothing to do with the actual output. os._exit skips that
    # cleanup entirely instead of chasing every remaining ownership edge case.


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GUNDAM pretty plots")
    parser.add_argument("--tag", default="asimov", help="Dataset tag to use for input/output paths")
    parser.add_argument("--pot", default=POT_DEFAULT, help="POT value to display on the plots")
    parser.add_argument("--data-dir", default=None,
                        help="Directory containing the input ROOT files (default: data/<var>/<tag>)")
    parser.add_argument("--output-dir", default=None,
                        help="Directory to write plots into (default: plots/<var>/<tag>)")
    parser.add_argument("--truth-reweight", default=None,
                        help="numpy expression over the truth-tree branches giving the "
                             "fake-data signal reweight for the xsec closure reference "
                             "(e.g. \"where(true_interaction_mode==1, 0.5, 1.0)\"). "
                             "Default: this tag's registered reweight, else nominal.")
    parser.add_argument("--truth-label", default=None,
                        help="Legend label for the truth reference curve.")
    parser.add_argument("--var", default=None, choices=sorted(analysis_vars.VARIABLES),
                        help="Analysis variable (default: $GUNDAM_VAR, else EleE)")
    args = parser.parse_args()
    if args.var:
        analysis_vars.set_var(args.var)

    main(tag=args.tag, pot=args.pot, data_dir=args.data_dir, output_dir=args.output_dir,
         truth_reweight=args.truth_reweight, truth_label=args.truth_label)
    os._exit(0)
