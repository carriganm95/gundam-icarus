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

POT_DEFAULT = "2.5e20"
STYLE_NAME = 'nue_icarus'
WORK_IN_PROGRESS = True

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


TEMPLATE_BINNING_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'binnings', 'binning_TrueEleE.txt')


def load_template_bin_labels(binning_file=TEMPLATE_BINNING_FILE):
    """Build 'E_{e} [lo, hi] MeV' labels for each template parameter from a GUNDAM binning text file."""
    labels = []
    with open(binning_file) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 2 or fields[0].startswith('variables'):
                continue
            lo, hi = float(fields[0]), float(fields[1])
            labels.append(f"E_{{e}} [{lo:g}, {hi:g}] MeV")
    return labels


def pretty_param_label(label, template_labels, template_mode=False):
    """Drop the '<parameter set>/' prefix and '#N_' index; name template parameters by their true-energy bin."""
    is_template = template_mode or label.startswith('Template Parameter')
    short = label.split('/')[-1]
    match = re.match(r'#(\d+)', short)
    if is_template and match and int(match.group(1)) < len(template_labels):
        return template_labels[int(match.group(1))]
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


def plotMatrix(tree, path, outputDir, plot_name, z_title, colorMap="kBird", pot=POT_DEFAULT):

    h = tree.Get(path)
    # Detach from the source file's directory so TFile::Close() at the end
    # of main() doesn't delete it out from under the canvas that's still
    # displaying it (SetOwnership below only guards against Python's GC,
    # not ROOT's own directory-based auto-cleanup of TH1-derived objects).
    if hasattr(h, "SetDirectory"):
        h.SetDirectory(0)

    r.gStyle.SetPalette(getattr(r, colorMap))
    r.gStyle.SetOptTitle(0)
    r.gStyle.SetOptStat(0)

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

def plotCovarianceMatrix(tree, outputDir, plot_name="covarianceMatrix", path="FitterEngine/postFit/Hesse/hessian/postfitCovariance_TH2D", pot=POT_DEFAULT):
    plotMatrix(tree, path, outputDir, plot_name, "Covariance", "kBird", pot=pot)

def plotCorrelationMatrix(tree, outputDir, plot_name="correlationMatrix", path="FitterEngine/postFit/Hesse/hessian/postfitCorrelation_TH2D", pot=POT_DEFAULT):
    plotMatrix(tree, path, outputDir, plot_name, "Correlation", "kThermometer", pot=pot)

def plotPulls(tree, path, outputDir, plot_name, pot=POT_DEFAULT):
    """Redraw GUNDAM's pre/post-fit parameter canvas on a fresh canvas: no title, rotated renamed labels."""

    c_par = tree.Get(path)
    r.SetOwnership(c_par, False)
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

def draw_ratio(ax, total_pred, err, data_y, data_yerr, x_edges, x_vals, center=1.0, ypad=(0.5, 1.5)):
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
    ax.set_xlabel("Reco Leading Electron Energy [MeV]")
    ax.set_ylabel("Data / Pred")

def get_mc_hist(file_obj, category_index, mode):
    """Get category histogram from ROOT file with path fallbacks across file variants."""
    candidates = []
    if mode == "pre":
        candidates = [
            f"toyGen/plots/histograms/signal_EleE (pre-fit)/reco_leading_electron_energy/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/signal_EleE/reco_leading_electron_energy/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/signal_EleE_pre_fit/reco_leading_electron_energy/category/{category_index}/MC_TH1D",
        ]
    else:
        candidates = [
            f"toyGen/plots/histograms/signal_EleE (post-fit)/reco_leading_electron_energy/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/signal_EleE/reco_leading_electron_energy/category/{category_index}/MC_TH1D",
            f"toyGen/plots/histograms/signal_EleE_post_fit/reco_leading_electron_energy/category/{category_index}/MC_TH1D",
        ]

    for path in candidates:
        obj = file_obj.Get(path)
        if obj and hasattr(obj, "GetBinContent"):
            return obj
    print("Warning no MC_TH1D found for category", category_index, "in", mode, "fit file. Tried paths:", candidates)
    return None

def getHistErrors(fPre, fPost):
    h_errors_pre = fPre.Get("toyGen/histograms/signal_EleE_pre_fit_TH1D")  # errors
    h_errors_post = fPost.Get("toyGen/histograms/signal_EleE_TH1D")  # errors

    h_data = fPost.Get("toyGen/plots/histograms/signal_EleE/reco_leading_electron_energy/Data_TH1D")  # data

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
            h_data_pre = get_mc_hist(fPre, cat, "pre")
            h_data_post = get_mc_hist(fPost, cat, "post")
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

def plot_pre_post_comparison(prefit_file, postfit_file, outputDir, pot=POT_DEFAULT):
    """Create side-by-side pre/post stacked plots with ratio panels and save image."""
    try:
        res = getHistErrors(prefit_file, postfit_file)
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

        draw_ratio(ax1_ratio, tot1, err1, data, data_err, edges1, x_vals)
        draw_ratio(ax2_ratio, tot2, err2, data, data_err, edges2, x_vals)

        # shared y-limits
        ymax1 = np.nanmax(tot1 + err1) if len(tot1) and len(err1) else 0.0
        ymax2 = np.nanmax(tot2 + err2) if len(tot2) and len(err2) else 0.0
        shared_ymax = max(ymax1, ymax2)
        pad = 0.12
        shared_ymax = max(1.0, shared_ymax * (1.0 + pad))
        ax1_main.set_ylim(0.0, shared_ymax)
        ax2_main.set_ylim(0.0, shared_ymax)
        ax1_main.legend(loc='upper right')
        ax2_main.legend(loc='upper right')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        stamp_ax(ax1_main, pot)
        stamp_ax(ax2_main, pot)

        fig.savefig(f"{outputDir}/png/pre_post_fit_comparison.png")
        with open(f"{outputDir}/pkl/pre_post_fit_comparison.pkl", "wb") as f:
            pickle.dump(fig, f)
        plt.close(fig)
    except Exception as e:
        print("Warning: failed to make pre/post fit comparison:", e)

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


def load_gundam_xsec(root_file, hist_name, unit_scale=1.0):
    """Read a GUNDAM cross-section TH1D and its stat errors."""
    with uproot.open(root_file) as f:
        h_xsec = f[hist_name]
        counts, edges = h_xsec.to_numpy()
        errs = h_xsec.errors()
    return counts * unit_scale, errs * unit_scale, edges


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
                          ratio_ylim=None, outpath=None, pot=POT_DEFAULT, fontsize=18):
    """Draw a two-panel cross-section comparison: a reference (e.g. GENIE) vs GUNDAM, with a ratio panel.

    Any '{scale}' in ylabel is replaced by the common power of ten, which is divided out of
    the plotted values so matplotlib's offset text doesn't collide with the watermark.
    """
    bin_widths = np.diff(edges)
    bin_centers = (edges[:-1] + edges[1:]) / 2

    exponent = int(np.floor(np.log10(np.nanmax(np.abs(ref_xsec)))))
    unit = 10.0 ** exponent
    ref_xsec, ref_xsec_err = ref_xsec / unit, ref_xsec_err / unit
    gundam_xsec, gundam_xsec_err = gundam_xsec / unit, gundam_xsec_err / unit
    ylabel = ylabel.replace('{scale}', rf'$10^{{{exponent}}}$')
    if ylim:
        ylim = (ylim[0] / unit, ylim[1] / unit)

    fig, (ax_xsec, ax_ratio) = plt.subplots(
        2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]}, sharex=True,
    )
    fig.subplots_adjust(hspace=0.05)

    ax_xsec.hist(edges[:-1], bins=edges, weights=ref_xsec, label=ref_label, histtype='step', alpha=1.0, color='blue', linewidth=2)
    ax_xsec.errorbar(bin_centers, gundam_xsec, yerr=gundam_xsec_err, xerr=bin_widths / 2,
                      fmt='o', capsize=5, label='GUNDAM', color='black', markersize=6)
    ax_xsec.legend(loc='upper right', fontsize=fontsize)
    ax_xsec.set_ylabel(ylabel)
    if ylim:
        ax_xsec.set_ylim(*ylim)
    else:
        top = np.nanmax([np.nanmax(ref_xsec), np.nanmax(gundam_xsec + gundam_xsec_err)])
        ax_xsec.set_ylim(0, top * 1.15)

    chi2_val, ndof, pvalue = _chi2_pvalue(ref_xsec, ref_xsec_err, gundam_xsec, gundam_xsec_err)
    chi2_ndof = chi2_val / ndof if ndof > 0 else np.nan
    ax_xsec.text(
        0.97, 0.78,
        r'$\chi^2$/ndof = %.1f/%d = %.2f' % (chi2_val, ndof, chi2_ndof) + '\n' + r'p = %.3f' % (pvalue),
        transform=ax_xsec.transAxes,
        fontsize=fontsize, ha='right', va='top',
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


def plot_genie_gundam_comparison(outputDir, gundam_xsec_file, pot=POT_DEFAULT):
    """Create the GENIE (NUISANCE flat tree) vs GUNDAM cross-section comparison plot and save it."""
    try:
        gundam_xsec_hist = "calcXsec/histograms/TrueEleE_XS_TH1D"

        # NUISANCE-flattened GENIE samples (nue, nuebar)
        nue_flat_trees = {
            'nue': '/exp/icarus/data/users/micarrig/genie/nue_samples_v3_04_00/p12_1000180400_CC_v3_4_0_AR23_20i_00_000.flat.root',
            'nuebar': '/exp/icarus/data/users/micarrig/genie/nue_samples_v3_04_00/p-12_1000180400_CC_v3_4_0_AR23_20i_00_000.flat.root',
        }

        # Same analysis binning used by the GUNDAM template (GeV)
        xsec_edges = np.array([0.300, 0.390, 0.490, 0.600, 0.730, 0.870, 1.030, 1.240, 1.500, 2.000, 5.000])

        nue_branches = ['ICARUS_1eInc_IsSignal', 'ICARUS_1eInc_electronEnergy', 'InputWeight', 'fScaleFactor']

        genie_xsecs, genie_xsec_errs, flux_integrals = [], [], []
        for file_path in nue_flat_trees.values():
            df = load_nuisflat_dataframe(file_path, nue_branches)
            xsec, xsec_err = compute_genie_xsec(
                df, 'ICARUS_1eInc_electronEnergy', xsec_edges,
                signal_selection_expr='(ICARUS_1eInc_IsSignal==True)',
            )
            genie_xsecs.append(xsec)
            genie_xsec_errs.append(xsec_err)
            flux_integrals.append(compute_flux_integral(file_path))

        genie_xsec, genie_xsec_err = combine_flux_averaged_xsec(genie_xsecs, genie_xsec_errs, flux_integrals)

        gundam_xsec, gundam_xsec_err, _ = load_gundam_xsec(gundam_xsec_file, gundam_xsec_hist, unit_scale=1000) # units scale to convert MeV->GeV

        fig, _ = plot_xsec_comparison(
            xsec_edges, genie_xsec, genie_xsec_err, gundam_xsec, gundam_xsec_err,
            xlabel=r'$E_{e}$ (GeV)',
            ylabel=r'$d\sigma/dE$ [{scale} cm$^{2}$ / GeV / Ar]',
            outpath=f"{outputDir}/png/genie_gundam_xsec_comparison.png",
            pot=pot,
        )
        with open(f"{outputDir}/pkl/genie_gundam_xsec_comparison.pkl", "wb") as f:
            pickle.dump(fig, f)
        plt.close(fig)
    except Exception as e:
        print("Warning: failed to make GENIE/GUNDAM cross-section comparison:", e)


def main(tag='asimov', pot=POT_DEFAULT):

    r.gROOT.SetBatch(True)
    plt.style.use(_style.get_style())

    outputDir = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/plots/{tag}'
    
    fitFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/fit_{tag}.root'
    # xsecFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/nue_asimovXSec.root'
    prefitFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/toy_pf_{tag}.root'
    postfitFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/toy_bf_{tag}.root'
    xsecFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/xsec_{tag}.root'

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

    # Template Parameters
    try:
        plotCovarianceMatrix(fit_file, outputDir, plot_name="templateParameterCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Template Parameter true_leading_electron_energy/matrices/Covariance_TH2D", pot=pot)
        plotCorrelationMatrix(fit_file, outputDir, plot_name="templateParameterCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Template Parameter true_leading_electron_energy/matrices/Correlation_TH2D", pot=pot)
        plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Template Parameter true_leading_electron_energy/values/fitConstraints_TCanvas", outputDir, "templateParameterFitPulls", pot=pot)
    except Exception as e:
        print("Warning: failed to make template parameter plots:", e)


    # Cross Section Parameters
    try:
        plotCovarianceMatrix(fit_file, outputDir, plot_name="crossSectionCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Cross section Systematics/matrices/Covariance_TH2D", pot=pot)
        plotCorrelationMatrix(fit_file, outputDir, plot_name="crossSectionCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Cross section Systematics/matrices/Correlation_TH2D", pot=pot)
        plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Cross section Systematics/values/fitConstraints_TCanvas", outputDir, "crossSectionFitPulls", pot=pot)
    except Exception as e:
        print("Warning: failed to make cross section parameter plots:", e)

    # Flux Parameters
    try:
        plotCovarianceMatrix(fit_file, outputDir, plot_name="fluxCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Flux Systematics Beamline/matrices/Covariance_TH2D", pot=pot)
        plotCorrelationMatrix(fit_file, outputDir, plot_name="fluxCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Flux Systematics Beamline/matrices/Correlation_TH2D", pot=pot)
        plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Flux Systematics Beamline/values/fitConstraints_TCanvas", outputDir, "fluxParameterFitPulls", pot=pot)
    except Exception as e:
        print("Warning: failed to make flux parameter plots:", e)

    # Detector Systematics
    try:
        plotCovarianceMatrix(fit_file, outputDir, plot_name="detectorSystematicsCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Detector Systematics/matrices/Covariance_TH2D", pot=pot)
        plotCorrelationMatrix(fit_file, outputDir, plot_name="detectorSystematicsCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Detector Systematics/matrices/Correlation_TH2D", pot=pot)
        plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Detector Systematics/values/fitConstraints_TCanvas", outputDir, "detectorSystematicsFitPulls", pot=pot)
    except Exception as e:
        print("Warning: failed to make detector systematics plots:", e)

    # Pre/Post fit stacked comparison (moved to helper)
    plot_pre_post_comparison(prefit_file, postfit_file, outputDir, pot=pot)

    # GENIE (CAF truth) vs GUNDAM cross-section comparison
    plot_genie_gundam_comparison(outputDir, xsecFile, pot=pot)

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
    args = parser.parse_args()

    main(tag=args.tag, pot=args.pot)
    os._exit(0)
