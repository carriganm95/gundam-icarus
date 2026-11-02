"""Per-variable settings for the nueCC GUNDAM analysis.

Everything that changes between the electron-energy (EleE) and the
electron-beam opening-angle (EleAngle) measurements lives here, so the
pipeline and plotting scripts can switch variable with one setting:

    GUNDAM_VAR=EleAngle ./run_fake_study.sh asimov

The variable is chosen by the GUNDAM_VAR environment variable (default EleE).
Plot scripts also take --var, which sets GUNDAM_VAR for the rest of the run.
"""

import os

DEFAULT_VAR = "EleE"

VARIABLES = {
    "EleE": {
        # GUNDAM sample names (set in Configs_FitSampleSet/fitSamples_*.yaml)
        "reco_sample": "signal_EleE",
        "xsec_sample": "TrueEleE_XS",
        "reco_branch": "reco_leading_electron_energy",
        "truth_branch": "true_leading_electron_energy",
        "template_set": "Template Parameter true_leading_electron_energy",
        "template_binning": "binning_TrueEleE.txt",
        # Branch values are in MeV; plots are in GeV.
        "unit_scale": 1e-3,
        "label": r"$E_{e}$ (GeV)",
        "reco_label": r"Reco $E_e$ [GeV]",
        "true_label": r"$E_e$ [GeV]",
        "reco_axis_title": "Reco Leading Electron Energy [MeV]",
        "xsec_ylabel": r"$d\sigma/dE$ [{scale} cm$^{2}$ / GeV / Ar]",
        # Template parameter labels: lo/hi are in branch units.
        "param_label": "E_{{e}} [{lo:g}, {hi:g}] MeV",
        # Independent GENIE (NUISANCE flat tree) overlay branch, in plot units.
        "genie_branch": "ICARUS_1eInc_electronEnergy",
        # NC pi0 sideband sample (fitSamples_EleE.yaml); reco_visible_energy is in GeV.
        "sideband_sample": "sideband_NCpi0",
        "sideband_branch": "reco_visible_energy",
        "sideband_axis_title": "Reco Visible Energy [GeV]",
    },
    "EleAngle": {
        "reco_sample": "signal_EleAngle",
        "xsec_sample": "TrueEleAngle_XS",
        "reco_branch": "reco_ele_beam_open_angle",
        "truth_branch": "true_ele_beam_open_angle",
        "template_set": "Template Parameter true_ele_beam_open_angle",
        "template_binning": "binning_TrueEleAngle.txt",
        # Branch is already cos(theta): no unit conversion.
        "unit_scale": 1.0,
        "label": r"$\cos\theta_{e}$",
        "reco_label": r"Reco $\cos\theta_e$",
        "true_label": r"$\cos\theta_e$",
        "reco_axis_title": r"Reco $\cos\theta_e$",
        "xsec_ylabel": r"$d\sigma/d\cos\theta_e$ [{scale} cm$^{2}$ / Ar]",
        "param_label": "cos#theta_{{e}} [{lo:g}, {hi:g}]",
        "genie_branch": "ICARUS_1eInc_cosElectronAngle",
        # NC pi0 sideband sample (fitSamples_EleAngle.yaml); reco_visible_energy is in GeV.
        "sideband_sample": "sideband_NCpi0",
        "sideband_branch": "reco_visible_energy",
        "sideband_axis_title": "Reco Visible Energy [GeV]",
    },
}


def set_var(name):
    """Select the analysis variable for the rest of this process (and child processes)."""
    if name not in VARIABLES:
        raise ValueError(f"unknown GUNDAM_VAR '{name}'; expected one of {sorted(VARIABLES)}")
    os.environ["GUNDAM_VAR"] = name


def apply_cli_override(argv=None):
    """Honour a '--var X' / '--var=X' on the command line before anything reads current().

    Scripts whose module-level constants depend on the variable call this right
    after importing analysis_vars (argparse still declares --var for help/validation).
    """
    import sys
    argv = sys.argv[1:] if argv is None else argv
    for i, a in enumerate(argv):
        if a == "--var" and i + 1 < len(argv):
            set_var(argv[i + 1])
        elif a.startswith("--var="):
            set_var(a.split("=", 1)[1])


def current_name():
    name = os.environ.get("GUNDAM_VAR", DEFAULT_VAR) or DEFAULT_VAR
    if name not in VARIABLES:
        raise ValueError(f"unknown GUNDAM_VAR '{name}'; expected one of {sorted(VARIABLES)}")
    return name


def current():
    """Settings dict for the selected variable, plus the derived histogram paths."""
    v = dict(VARIABLES[current_name()])
    v["name"] = current_name()
    v["xsec_hist"] = f"calcXsec/histograms/{v['xsec_sample']}_TH1D"
    v["dist_hist"] = f"toyGen/histograms/{v['reco_sample']}_TH1D"
    v["dist_prefit_hist"] = f"toyGen/histograms/{v['reco_sample']}_pre_fit_TH1D"
    return v
