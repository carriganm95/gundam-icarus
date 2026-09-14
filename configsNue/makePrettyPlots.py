import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import os

def plotMatrix(tree, path, outputDir, title, plot_name, colorMap="kBird"):

    c = r.TCanvas(f"c_{plot_name}", f"c_{plot_name}", 2800, 2000)
    r.SetOwnership(c, False)
    h = tree.Get(path)
    # Detach from the source file's directory so TFile::Close() at the end
    # of main() doesn't delete it out from under the canvas that's still
    # displaying it (SetOwnership above only guards against Python's GC,
    # not ROOT's own directory-based auto-cleanup of TH1-derived objects).
    if hasattr(h, "SetDirectory"):
        h.SetDirectory(0)

    r.gStyle.SetPalette(getattr(r, colorMap))

    h.GetXaxis().SetLabelSize(0.01)
    h.GetYaxis().SetLabelSize(0.01)
    h.Draw("COLZ")
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.2)
    c.SetRightMargin(0.13)
    c.SetRightMargin(0.13)
    c.Draw()
    r.gStyle.SetImageScaling(3.)

    h.SetTitle(f"{title};;")
    c.SaveAs(f"{outputDir}/{plot_name}.png")

def plotCovarianceMatrix(tree, outputDir, title="Covariance Matrix", plot_name="covarianceMatrix", path="FitterEngine/postFit/Hesse/hessian/postfitCovariance_TH2D"):
    plotMatrix(tree, path, outputDir, title, plot_name, "kBird")

def plotCorrelationMatrix(tree, outputDir, title="Correlation Matrix", plot_name="correlationMatrix", path="FitterEngine/postFit/Hesse/hessian/postfitCorrelation_TH2D"):
    plotMatrix(tree, path, outputDir, title, plot_name, "kThermometer")

def plotPulls(tree, path, outputDir, plot_name):

    c_par = tree.Get(path)
    r.SetOwnership(c_par, False)
    c_par.SetCanvasSize(2000, 1200)

    if 'Template' in path:
        primitives = list(c_par.GetListOfPrimitives())  # make a Python list copy

        # Identify your objects by name or type
        h_prefit = c_par.GetPrimitive("preFitErrors")
        h_postfit = c_par.GetPrimitive("postFitErrors")
        h_line = c_par.GetPrimitive("preFitErrorHistLine")
        tpave = c_par.GetPrimitive("TPave")  # or get by index/type

        # Clone all objects BEFORE drawing, and tell cppyy that ROOT (not
        # Python's refcounting) owns them, since they'll be attached to a
        # canvas below. Without this, Python's GC and the canvas's own
        # primitive-list cleanup can both try to free the same object,
        # causing "double free or corruption" at exit.
        objects_to_draw = []

        def clone_owned(obj):
            clone = obj.Clone()
            r.SetOwnership(clone, False)
            # Same TFile-directory auto-registration issue as in plotMatrix:
            # detach histogram-like clones so the source file's Close() call
            # can't delete them while they're still drawn on our canvas.
            if hasattr(clone, "SetDirectory"):
                clone.SetDirectory(0)
            return clone

        if h_prefit:
            objects_to_draw.append(clone_owned(h_prefit))
        if h_postfit:
            objects_to_draw.append(clone_owned(h_postfit))
        if h_line:
            objects_to_draw.append(clone_owned(h_line))
        if tpave:
            objects_to_draw.append(clone_owned(tpave))

        for i in range(2, 12):
            objects_to_draw.append(clone_owned(primitives[i]))

        # Now create a new canvas
        c = r.TCanvas(f"c_{plot_name}", f"c_{plot_name}", 2000, 1200)
        r.SetOwnership(c, False)

        # Draw the cloned prefit histogram (not the original, which is
        # still owned by c_par/the source file) as the base frame.
        objects_to_draw[0].Draw()

        # Draw the remaining objects in the correct order
        for obj in objects_to_draw[1:][::-1]:
            obj.Draw("same")


        c.Modified()
        c.Draw()
        c.SaveAs(f"{outputDir}/{plot_name}.png")

    else:
        c_par.Draw()
        c_par.SaveAs(f"{outputDir}/{plot_name}.png")

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
    ax.set_title(title)
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
    raise RuntimeError(f"Could not find MC_TH1D for category {category_index} in {mode}-fit file")

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
    h_mc_pre = []
    h_mc_post = []
    # assume 5 categories as in the notebook
    for cat in range(0, 5):
        pre_bin_content = []
        post_bin_content = []
        h_data_pre = get_mc_hist(fPre, cat, "pre")
        h_data_post = get_mc_hist(fPost, cat, "post")
        for j in range(h_data_pre.GetNbinsX()):
            pre_bin_content.append(h_data_pre.GetBinContent(j + 1))
            post_bin_content.append(h_data_post.GetBinContent(j + 1))
        h_mc_pre.append(pre_bin_content)
        h_mc_post.append(post_bin_content)

    data = []
    data_err = []
    for i in range(1, h_data.GetNbinsX() + 1):
        data.append(h_data.GetBinContent(i))
        data_err.append(h_data.GetBinError(i))

    # labels and colors matching notebook
    labels = [r"$\nu_e$ Signal", r"$\nu_e$ OOFV", r"$\nu_\mu$", r"$\nu$ OOFV", r"$\nu$ NC"]
    category_colors = ['#011959', '#1C5A62', '#687B3E', '#D29343', '#FDB7BC']

    n_bins = len(x_vals)

    return x_vals, widths, preFitErrors, postFitErrors, np.array(h_mc_pre), np.array(h_mc_post), data, data_err, n_bins, labels, category_colors

def plot_pre_post_comparison(prefit_file, postfit_file, outputDir):
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
        plt.tight_layout()

        outpath = f"{outputDir}/pre_post_fit_comparison.png"
        fig.savefig(outpath)
        plt.close(fig)
    except Exception as e:
        print("Warning: failed to make pre/post fit comparison:", e)

def main(tag='asimov'):

    r.gROOT.SetBatch(True)

    outputDir = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/plots/{tag}'
    
    fitFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/fit_{tag}.root'
    # xsecFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/nue_asimovXSec.root'
    prefitFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/toy_pf_{tag}.root'
    postfitFile = f'/nashome/m/micarrig/icarus/gundam/gundam-icarus/configsNue/data/{tag}/toy_bf_{tag}.root'

    fit_file = r.TFile(fitFile)
    # xsec_file = r.TFile(xsecFile)
    prefit_file = r.TFile(prefitFile)
    postfit_file = r.TFile(postfitFile)

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Total
    plotCovarianceMatrix(fit_file, outputDir)
    plotCorrelationMatrix(fit_file, outputDir)
    
    # Template Parameters
    plotCovarianceMatrix(fit_file, outputDir, title="Template Parameter Covariance Matrix", plot_name="templateParameterCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Template Parameter true_leading_electron_energy/matrices/Covariance_TH2D")
    plotCorrelationMatrix(fit_file, outputDir, title="Template Parameter Correlation Matrix", plot_name="templateParameterCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Template Parameter true_leading_electron_energy/matrices/Correlation_TH2D")
    plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Template Parameter true_leading_electron_energy/values/fitConstraints_TCanvas", outputDir, "templateParameterFitPulls")

    # Cross Section Parameters
    plotCovarianceMatrix(fit_file, outputDir, title="Cross Section Covariance Matrix", plot_name="crossSectionCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Cross section Systematics/matrices/Covariance_TH2D")
    plotCorrelationMatrix(fit_file, outputDir, title="Cross Section Correlation Matrix", plot_name="crossSectionCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Cross section Systematics/matrices/Correlation_TH2D")
    plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Cross section Systematics/values/fitConstraints_TCanvas", outputDir, "crossSectionFitPulls")

    # Flux Parameters
    plotCovarianceMatrix(fit_file, outputDir, title="Flux Covariance Matrix", plot_name="fluxCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Flux Systematics Beamline/matrices/Covariance_TH2D")
    plotCorrelationMatrix(fit_file, outputDir, title="Flux Correlation Matrix", plot_name="fluxCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Flux Systematics Beamline/matrices/Correlation_TH2D")
    plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Flux Systematics Beamline/values/fitConstraints_TCanvas", outputDir, "fluxParameterFitPulls")

    # Detector Systematics
    plotCovarianceMatrix(fit_file, outputDir, title="Detector Systematics Covariance Matrix", plot_name="detectorSystematicsCovarianceMatrix", path="FitterEngine/postFit/Hesse/errors/Detector Systematics/matrices/Covariance_TH2D")
    plotCorrelationMatrix(fit_file, outputDir, title="Detector Systematics Correlation Matrix", plot_name="detectorSystematicsCorrelationMatrix", path="FitterEngine/postFit/Hesse/errors/Detector Systematics/matrices/Correlation_TH2D")
    plotPulls(fit_file, "FitterEngine/postFit/Hesse/errors/Detector Systematics/values/fitConstraints_TCanvas", outputDir, "detectorSystematicsFitPulls")

    # Pre/Post fit stacked comparison (moved to helper)
    plot_pre_post_comparison(prefit_file, postfit_file, outputDir)

    # Deliberately skip TFile.Close()/normal interpreter teardown: all PNGs
    # are already on disk by this point, the files were opened read-only, and
    # ROOT/cppyy's own exit-time cleanup of the many canvases and cloned
    # primitives created above triggers a "double free or corruption" abort
    # that has nothing to do with the actual output. os._exit skips that
    # cleanup entirely instead of chasing every remaining ownership edge case.


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GUNDAM pretty plots")
    parser.add_argument("--tag", default="asimov", help="Dataset tag to use for input/output paths")
    args = parser.parse_args()

    main(tag=args.tag)
    os._exit(0)