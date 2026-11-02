#!/usr/bin/env python3
"""
calc_geant4_covariance.py

Computes a fractional bin-to-bin covariance matrix for one GEANT4 FSI
"fate" reweight dial (proton, pi+, or pi-), following the same method as
Kiyoung's CovMat_kjung/CalcCov_test.ipynb, but adapted to the actual
layout found in icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx_G4.root:

  - The event-level analysis variables live in one tree (MainTree,
    e.g. events/full/selected).
  - The per-event GEANT4 reweight universes live in a SEPARATE tree in
    the same file (DialTree, e.g. events/full/selected_multisimTree),
    as a plain std::vector<double> branch of length NUniv per event
    (e.g. reinteractions_proton_Geant4) -- NOT a TClonesArray of TGraph
    objects like Kiyoung's file, so no .member('fY') gymnastics needed.
  - Verified directly against the file (2026-09-10): MainTree and
    DialTree have the same number of entries and are row-for-row
    aligned by (Run, Subrun, Evt). This script re-checks that alignment
    every time it runs rather than assuming it holds -- if it doesn't,
    it stops instead of silently mismatching events to universes.

Method, same as Kiyoung's:
  1. Bin the selected events (nominal, unweighted-by-dial) using a
     binning file (see BinParser.py).
  2. For every universe u, recompute the per-bin weighted yield using
     that universe's per-event reweight value.
  3. Covariance is the fractional bin-to-bin covariance of the universe
     yields around the CV (mean-over-universes) yield:
         Cov[i][j] = sum_u (N_u[i]-Nbar[i])(N_u[j]-Nbar[j]) / (NUniv-1)
                     / Nbar[i] / Nbar[j]
     This is Kiyoung's normalization convention (fractional, divided by
     the universe MEAN, not the CV/nominal count) -- carried over as-is,
     not re-derived independently. If GUNDAM expects the raw (non
     -fractional) covariance of Eq. 12 in main.pdf instead, this needs to
     change before it's usable as a parameter-set input; flag that back
     to me once you know which convention the GEANT4 ParameterSet in
     Configs_ParameterSet/GEANT4 actually consumes.

Usage:
    python3 calc_geant4_covariance.py config_G4Proton.json

Requires: uproot, numpy, pandas, matplotlib, and PyROOT (for writing the
TMatrixDSym output). If PyROOT is not available in your environment, tell
me and I'll change the output step to something that doesn't need it
(e.g. .npz).
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import uproot

import BinParser


def load_config(config_path):
    with open(config_path) as f:
        cfg = json.load(f)

    required_keys = [
        "MainTree", "DialTree", "BinningFile", "JobName", "DialName", "NUniv",
        "SelectionCutExpr", "SampleSelectionExpr", "SampleWeightExpr",
        "OutputFileName",
    ]
    missing = [k for k in required_keys if k not in cfg]
    if missing:
        raise KeyError("Config %s is missing required key(s): %s" % (config_path, missing))
    for sub_key in ["FilePath", "TreeKey"]:
        if sub_key not in cfg["MainTree"]:
            raise KeyError("Config %s: MainTree is missing '%s'" % (config_path, sub_key))
    if "TreeKey" not in cfg["DialTree"]:
        raise KeyError("Config %s: DialTree is missing 'TreeKey'" % (config_path,))

    return cfg


def load_tree_and_universes(cfg):
    file_path = cfg["MainTree"]["FilePath"]
    tree_key = cfg["MainTree"]["TreeKey"]
    # DialTree can override FilePath if the universes ever live in a
    # different file; defaults to the same file as MainTree, which is the
    # layout actually found in icarus_..._ppfx_G4.root.
    dial_file_path = cfg["DialTree"].get("FilePath", file_path)
    dial_tree_key = cfg["DialTree"]["TreeKey"]
    dial_name = cfg["DialName"]
    n_univ_expected = cfg["NUniv"]

    print("- InputFileName: %s" % file_path)
    print("- TreeKey: %s" % tree_key)
    print("- DialTreeKey: %s (file: %s)" % (dial_tree_key, dial_file_path))
    print("- DialName: %s" % dial_name)

    with uproot.open(file_path) as root_file:
        tree = root_file[tree_key]
        df = tree.arrays(tree.keys(), library="pd")

    # Guarantee df's index is the plain row position (0..N-1). Everything
    # below indexes the universe array by this same position, so this
    # assumption has to hold exactly, not "probably hold".
    if not np.array_equal(df.index.values, np.arange(len(df))):
        raise RuntimeError(
            "DataFrame index is not a plain 0..N-1 range after reading "
            "MainTree. The per-event universe lookup below assumes "
            "positional indexing; call df = df.reset_index(drop=True) "
            "before using it, or tell me if the tree read is expected to "
            "reorder/filter rows."
        )

    with uproot.open(dial_file_path) as dial_root_file:
        dial_tree = dial_root_file[dial_tree_key]

        if dial_name not in dial_tree.keys():
            raise KeyError(
                "Branch '%s' (DialName) not found in DialTree %s. Branches present: %s"
                % (dial_name, dial_tree_key, dial_tree.keys())
            )

        # Cross-check MainTree and DialTree are row-for-row the same events
        # -- required because the covariance code below indexes both by
        # plain row position, not by any join key. Verified true for
        # events/full/selected vs events/full/selected_multisimTree in
        # icarus_..._ppfx_G4.root on 2026-09-10; re-checked here every run
        # rather than trusted from that one check.
        id_branches = [b for b in ["Run", "Subrun", "Evt"] if b in df.columns and b in dial_tree.keys()]
        if id_branches:
            dial_ids = dial_tree.arrays(id_branches, library="pd")
            if len(dial_ids) != len(df):
                raise RuntimeError(
                    "MainTree (%d entries) and DialTree (%d entries) have "
                    "different lengths -- they are not the same set of "
                    "events, can't align by position." % (len(df), len(dial_ids))
                )
            for b in id_branches:
                if not np.array_equal(df[b].values, dial_ids[b].values):
                    raise RuntimeError(
                        "MainTree and DialTree disagree on branch '%s' at "
                        "some row -- they are not row-for-row aligned. "
                        "Do not trust positional indexing here; the join "
                        "needs to be done explicitly by (Run, Subrun, Evt) "
                        "instead of assumed." % b
                    )
        else:
            print(
                "[WARN] Could not find Run/Subrun/Evt in both trees to "
                "verify row alignment -- proceeding on the unverified "
                "assumption that MainTree and DialTree are row-for-row "
                "the same events. Tell me if that's not guaranteed."
            )

        arr = dial_tree[dial_name].array(library="np")
        # Each arr[i] is already a plain 1D array of NUniv floats (a
        # std::vector<double> branch) -- unlike Kiyoung's file, there is no
        # TClonesArray/TGraph wrapper to unpack here.
        n_events = len(arr)
        n_univ_observed = len(arr[0])
        lengths = np.array([len(x) for x in arr])
        if not np.all(lengths == n_univ_observed):
            raise ValueError(
                "Branch '%s' does not store the same number of universes "
                "for every event (found lengths %s..%s) -- can't stack "
                "into a rectangular array without deciding how to handle "
                "that first." % (dial_name, lengths.min(), lengths.max())
            )
        if n_univ_observed != n_univ_expected:
            raise ValueError(
                "Config says NUniv=%d but branch '%s' actually stores %d "
                "universes per event. Fix NUniv in the config to match the "
                "file, don't silently truncate/pad."
                % (n_univ_expected, dial_name, n_univ_observed)
            )

        all_RW = np.stack(arr).astype(np.float64)  # (n_events, n_univ)

    print("DataFrame shape:", df.shape)
    print("Universe reweight array shape:", all_RW.shape)
    return df, all_RW


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="Path to the per-dial JSON config (see config_*_TEMPLATE.json)")
    parser.add_argument("--outdir", default=None, help="Override output directory (default: alongside the config, in ./outputs and ./plots)")
    parser.add_argument("--no-plots", action="store_true", help="Skip the matplotlib covariance/correlation plots")
    args = parser.parse_args()

    cfg = load_config(args.config)

    workdir = args.outdir or os.path.dirname(os.path.abspath(args.config)) or "."
    out_dir = os.path.join(workdir, "outputs")
    plot_dir = os.path.join(workdir, "plots")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    df, all_RW = load_tree_and_universes(cfg)

    # Nominal (CV) event weight
    sample_weight_expr = cfg["SampleWeightExpr"]
    print("- Sample weight: %s" % sample_weight_expr)
    df = df.eval("FinalWeight = %s" % sample_weight_expr)

    # Binning
    binning_file = cfg["BinningFile"]
    print("- Binning from %s:" % binning_file)
    bin_parser = BinParser.BinParser(binning_file)
    n_bins = len(bin_parser.CutExprs)
    for i_cut, cut_expr in enumerate(bin_parser.CutExprs):
        print("bin %d: " % i_cut, cut_expr)

    selection_cut_expr = cfg["SelectionCutExpr"]
    sample_selection_expr = cfg["SampleSelectionExpr"]

    # CV (nominal, unweighted-by-dial) per-bin yield -- sanity reference only
    h_nom = np.zeros(n_bins)
    for i_cut, cut_expr in enumerate(bin_parser.CutExprs):
        exprs = [e for e in [selection_cut_expr, sample_selection_expr, cut_expr] if e.strip()]
        this_cut_expr = " & ".join("(%s)" % e for e in exprs)
        h_nom[i_cut] = np.sum(df.query(this_cut_expr)["FinalWeight"])
    print("CV per-bin yield:", h_nom)

    # Per-universe per-bin yield
    job_name = cfg["JobName"]
    dial_name = cfg["DialName"]
    n_univ = all_RW.shape[1]

    h_univs = np.zeros((n_univ, n_bins))
    for i_cut, cut_expr in enumerate(bin_parser.CutExprs):
        exprs = [e for e in [selection_cut_expr, sample_selection_expr, cut_expr] if e.strip()]
        this_cut_expr = " & ".join("(%s)" % e for e in exprs)

        df_selected = df.query(this_cut_expr)
        idx = df_selected.index.values
        rw_selected = all_RW[idx]  # (n_selected, n_univ)
        weights = df_selected["FinalWeight"].values  # (n_selected,)

        # h_univs[:, i_cut] = sum over selected events of weight * per-universe reweight
        h_univs[:, i_cut] = weights @ rw_selected  # (n_univ,)

    # Fractional covariance / correlation (Kiyoung's convention -- see module docstring)
    h_cov = np.zeros((n_bins, n_bins))
    h_mean = np.mean(h_univs, axis=0)
    for i in range(n_bins):
        for j in range(n_bins):
            h_cov[i][j] = (
                np.sum((h_univs[:, i] - h_mean[i]) * (h_univs[:, j] - h_mean[j]))
                / (n_univ - 1) / h_mean[i] / h_mean[j]
            )

    h_corr = np.zeros((n_bins, n_bins))
    for i in range(n_bins):
        for j in range(n_bins):
            h_corr[i][j] = h_cov[i][j] / np.sqrt(h_cov[i, i]) / np.sqrt(h_cov[j, j])

    print("Universe-mean / CV ratio per bin (should be close to 1):", h_mean / h_nom)

    # --- Write ROOT output (TMatrixDSym) ---
    output_file_name = cfg["OutputFileName"]
    try:
        import ROOT
    except ImportError:
        print(
            "\n[WARN] PyROOT is not importable in this environment -- skipping the "
            ".root write. Saving covariance/correlation as .npz instead so nothing "
            "is lost; tell me if you need the TMatrixDSym ROOT file specifically "
            "and I'll adjust this.\n"
        )
        npz_path = os.path.join(out_dir, output_file_name.replace(".root", ".npz"))
        np.savez(npz_path, covariance=h_cov, correlation=h_corr, cv_yield=h_nom, mean_yield=h_mean)
        print("Wrote", npz_path)
    else:
        out_path = os.path.join(out_dir, output_file_name)
        f_out = ROOT.TFile(out_path, "RECREATE")

        h_cov_tmat = ROOT.TMatrixDSym(n_bins)
        for i in range(n_bins):
            for j in range(n_bins):
                h_cov_tmat[i][j] = h_cov[i, j]
        h_cov_tmat.SetTol(1e-20)
        print("Covariance matrix determinant:", h_cov_tmat.Determinant())

        h_corr_tmat = ROOT.TMatrixDSym(n_bins)
        for i in range(n_bins):
            for j in range(n_bins):
                h_corr_tmat[i][j] = h_corr[i, j]
        h_corr_tmat.SetTol(1e-20)

        f_out.mkdir(job_name)
        f_out.cd(job_name)
        h_cov_tmat.Write("%s_covariance" % job_name)
        h_corr_tmat.Write("%s_correlation" % job_name)

        # Per-bin universe distributions (yield / CV), the inputs to the covariance
        for i in range(n_bins):
            ratio = h_univs[:, i] / h_nom[i]
            half = max(abs(ratio.min() - 1), abs(ratio.max() - 1), 1e-4) * 1.1
            h1 = ROOT.TH1D("%s_univ_bin%d" % (job_name, i),
                           "%s;universe yield / CV;universes" % bin_parser.CutExprs[i],
                           50, 1 - half, 1 + half)
            for r in ratio:
                h1.Fill(r)
            h1.Write()
        h_cv = ROOT.TH1D("%s_cv_yield" % job_name, "CV yield;bin;weighted events", n_bins, 0, n_bins)
        for i in range(n_bins):
            h_cv.SetBinContent(i + 1, h_nom[i])
        h_cv.Write()

        # Raw per-universe yields (n_univ x n_bins), read by plot_g4_covariance.py
        m_univ = ROOT.TMatrixD(n_univ, n_bins)
        for u in range(n_univ):
            for i in range(n_bins):
                m_univ[u][i] = h_univs[u, i]
        m_univ.Write("%s_univ_yields" % job_name)
        ROOT.TNamed("%s_binning" % job_name, binning_file).Write()
        f_out.cd()
        f_out.Close()
        print("Wrote", out_path)

    # --- Plots (styled, via plot_g4_covariance.py) ---
    if not args.no_plots and "out_path" in locals():
        import subprocess
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_g4_covariance.py")
        subprocess.run([sys.executable, script, out_path, "--job", job_name, "--outdir", plot_dir], check=True)
        # CV vs G4-universe distribution of the binned variable (uses the events in memory)
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import plot_g4_covariance
        plot_g4_covariance.plot_distribution(df, all_RW, cfg, binning_file, job_name, plot_dir)

if __name__ == "__main__":
    main()
