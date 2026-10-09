#!/usr/bin/env bash
# Run fitter + calcXsec + toy generator for each systematic subset (individual
# and N-1), then generate fractional uncertainty comparison plots.

set -euo pipefail

# Parse flags: optional --diagonal, --zero-flux-norm, positional POT value
diagonal=false
zero_flux_norm=false
pot="2.5e20"
for arg in "$@"; do
    case "$arg" in
        --diagonal)                   diagonal=true ;;
        --zero-flux-norm|--noflux)    zero_flux_norm=true ;;
        *)                            pot="$arg" ;;
    esac
done

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
export GUNDAM_CONFIG_DIR="${script_dir}"

fitter_base="${script_dir}/fakeStudyConfigs/config_Fitter_FakeData_EleE_asimov.yaml"
xsec_base="${script_dir}/config_CalcXSec_FakeData_EleE.yaml"
toy_base="${script_dir}/fakeStudyConfigs/config_Toy_FakeData_EleE_asimov.yaml"
fit_sample_file="${script_dir}/Configs_FitSampleSet/fitSamples_EleE.yaml"
xs_sample_file="${script_dir}/Configs_FitSampleSet/fitSamples_TrueEleE.yaml"

if $diagonal; then
    # Use true energy as the "reco" variable — response matrix is diagonal by construction.
    binning_ele="binning_TrueEleE.txt"
    binning_nominal="binning_TrueEleE.txt"
    output_base="${script_dir}/data/frac_uncertainty_asimov_diagonal"
    plot_dir="${script_dir}/plots/frac_uncertainty_asimov_diagonal"
else
    binning_ele="binning_EleE.txt"
    binning_nominal="binning_EleENominal.txt"
    output_base="${script_dir}/data/frac_uncertainty_asimov"
    plot_dir="${script_dir}/plots/frac_uncertainty_asimov"
fi

# --zero-flux-norm: set the "Constant flux" normalisation stdDev to 0 in the xsec
# sample def, removing the 9.13% flux-normalisation uncertainty baked into every
# cross-section bin by gundamCalcXsec. Only affects calcXsec output, not the fit.
if $zero_flux_norm; then
    output_base="${output_base}_noflux"
    plot_dir="${plot_dir}_noflux"
fi

mkdir -p "${plot_dir}"

# fitSamples_EleE.yaml is modified in-place to switch binning.
# In normal mode: fitter uses 20-bin EleE, toy uses 10-bin EleENominal (matches calcXsec).
# In diagonal mode: both use binning_TrueEleE.txt (true energy = reco energy → diagonal matrix).
set_binning() {
  python3 - "${fit_sample_file}" "$1" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
bfile = sys.argv[2]
path.write_text(re.sub(r'^(\s+binning: ").*(")', r'\1${GUNDAM_CONFIG_DIR}/binnings/' + bfile + r'\2',
                path.read_text(), flags=re.M))
PY
}

set_ele_binning()     { set_binning "${binning_ele}"; }
set_nominal_binning() { set_binning "${binning_nominal}"; }

# Individual-contribution variants always run. The N-1 variants are skipped in
# diagonal mode: with a diagonal response the free per-bin templates are exactly
# degenerate with any reweighting systematic, so those fits are ill-posed (Migrad
# diverges to NaN). The individual variants are all that's needed for the
# individual-contributions plots.
VARIANTS=(all template_only xsec_only flux_only detector_only)
if ! $diagonal; then
    VARIANTS+=(no_xsec no_flux no_detector)
fi

tmp_fitter=""
tmp_xsec=""
tmp_toy=""
xs_backup=""
cleanup() {
    [[ -n "${tmp_fitter}" ]] && rm -f "${tmp_fitter}"
    [[ -n "${tmp_xsec}" ]]   && rm -f "${tmp_xsec}"
    [[ -n "${tmp_toy}" ]]    && rm -f "${tmp_toy}"
    set_ele_binning  # restore fitter binning on exit
    # restore the original xsec sample def if we zeroed the flux norm
    [[ -n "${xs_backup}" ]] && mv -f "${xs_backup}" "${xs_sample_file}"
}
trap cleanup EXIT ERR

set_ele_binning  # ensure clean state before starting

if $zero_flux_norm; then
    xs_backup=$(mktemp --suffix=.yaml)
    cp "${xs_sample_file}" "${xs_backup}"
    python3 - "${xs_sample_file}" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
# Zero the stdDev of the "Constant flux" normalise parameter (first stdDev after it).
text, n = re.subn(r'(name:\s*"Constant flux".*?stdDev:\s*)[0-9.eE+\-]+',
                  r'\g<1>0.0', path.read_text(), count=1, flags=re.S)
assert n == 1, "could not find 'Constant flux' stdDev to zero"
path.write_text(text)
print("Zeroed 'Constant flux' stdDev in", path.name)
PY
fi

for variant in "${VARIANTS[@]}"; do
    echo "============================================"
    echo " Running variant: ${variant}"
    echo "============================================"

    out_dir="${output_base}/${variant}"
    mkdir -p "${out_dir}"

    tmp_fitter=$(mktemp --suffix=.yaml)
    tmp_xsec=$(mktemp --suffix=.yaml)
    tmp_toy=$(mktemp --suffix=.yaml)

    sed "s|Configs_ParameterSetList/parameterSetList_EleE|frac_uncertainty/parameterSetList_${variant}|" \
        "${fitter_base}" > "${tmp_fitter}"
    sed "s|Configs_ParameterSetList/parameterSetList_EleE|frac_uncertainty/parameterSetList_${variant}|" \
        "${xsec_base}" > "${tmp_xsec}"
    sed "s|Configs_ParameterSetList/parameterSetList_EleE|frac_uncertainty/parameterSetList_${variant}|" \
        "${toy_base}" > "${tmp_toy}"

    # Process each variant tolerantly: a failed fit (e.g. a degenerate variant
    # where Migrad diverges to NaN) is logged and skipped, never aborts the run.
    variant_ok=true

    if ! gundamFitter -c "${tmp_fitter}" -o "${out_dir}/fit.root" --asimov; then
        echo "WARNING: gundamFitter failed for '${variant}' -- skipping this variant"
        variant_ok=false
    fi

    if ${variant_ok}; then
        if ! gundamCalcXsec -c "${tmp_xsec}" -f "${out_dir}/fit.root" \
                -o "${out_dir}/xsec.root" -n 1000 --use-bf-as-xsec; then
            echo "WARNING: gundamCalcXsec failed for '${variant}' -- skipping this variant"
            variant_ok=false
        fi
    fi

    if ${variant_ok}; then
        set_nominal_binning  # switch to 10-bin nominal for toy generator (matches calcXsec)
        gundamToyGenerator -c "${tmp_toy}" -f "${out_dir}/fit.root" \
            -n 1000 --use-bf --save-data -o "${out_dir}/toy_bf.root" \
            || echo "WARNING: toy (--use-bf) failed for '${variant}'"
        gundamToyGenerator -c "${tmp_toy}" -f "${out_dir}/fit.root" \
            -n 1000 --use-prefit --save-data -o "${out_dir}/toy_pf.root" \
            || echo "WARNING: toy (--use-prefit) failed for '${variant}'"
    fi

    rm -f "${tmp_fitter}" "${tmp_xsec}" "${tmp_toy}"
    tmp_fitter=""
    tmp_xsec=""
    tmp_toy=""

    set_ele_binning  # restore for next variant's fitter
done

echo "============================================"
echo " Generating comparison plots"
echo "============================================"

python3 "${script_dir}/scripts/plot_frac_uncertainty.py" \
    --data-dir "${output_base}" \
    --output-dir "${plot_dir}" \
    --pot "${pot}"

echo "Done. Plots saved to ${plot_dir}"
