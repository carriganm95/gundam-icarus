#!/usr/bin/env bash

set -euo pipefail

tag="${1:-asimov}"
pot="${2:-2.5e20}"
# Number of toys for the toy generator and calcXsec steps. Defaults to 1000;
# override with e.g. NTOYS=10 for a quick pipeline test.
n_toys="${NTOYS:-1000}"
output="${3:-${OUTPUT:-}}"

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
export GUNDAM_CONFIG_DIR="${script_dir}"
# Analysis variable: EleE (electron energy, default) or EleAngle (cos theta_e).
# Picks the *_${var}* configs/binnings and is inherited by makePrettyPlots.py.
var="${GUNDAM_VAR:-EleE}"
export GUNDAM_VAR="${var}"
fit_sample_file="${script_dir}/Configs_FitSampleSet/fitSamples_${var}.yaml"
# Each variable has its own tree: data/<var>/<tag> and plots/<var>/<tag>.
if [[ -n "${output}" ]]; then
  output_dir="${output}"
  plot_dir="${output}/plots"
else
  output_dir="${script_dir}/data/${var}/${tag}"
  plot_dir="${script_dir}/plots/${var}/${tag}"
fi
# PLOT_DIR overrides where the plots go (run_data_10pct.sh uses plots/<var>/data_10pct).
plot_dir="${PLOT_DIR:-${plot_dir}}"
config_dir="${script_dir}/fakeStudyConfigs"
# Per-step configs default to the tag-based fake-study configs; override with
# FIT_CONFIG / TOY_CONFIG / XSEC_CONFIG (run_data_10pct.sh overrides XSEC_CONFIG).
fit_config="${FIT_CONFIG:-${config_dir}/config_Fitter_FakeData_${var}_${tag}.yaml}"
toy_config="${TOY_CONFIG:-${config_dir}/config_Toy_FakeData_${var}_${tag}.yaml}"
xsec_config="${XSEC_CONFIG:-${script_dir}/config_CalcXSec_FakeData_${var}.yaml}"

mkdir -p "${output_dir}" "${plot_dir}"

if [[ ! -f "${fit_sample_file}" ]]; then
  echo "Missing fit sample file: ${fit_sample_file}" >&2
  exit 1
fi

restore_default_binning() {
  python3 - "${fit_sample_file}" "binning_${var}.txt" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
text = re.sub(
    r'^    binning: "(?!.*Sideband).*"$',  # signal sample only; sideband keeps its binning
    '    binning: "${GUNDAM_CONFIG_DIR}/binnings/' + sys.argv[2] + '"',
    text,
    flags=re.M,
)
# Sideband: fine reco binning for the fit, coarse nominal binning for toys/xsec.
text = re.sub(
    r'binnings/binning_SidebandVisE(Nominal)?\.txt',
    'binnings/binning_SidebandVisE.txt',
    text,
)
path.write_text(text)
PY
}

set_nominal_binning() {
  python3 - "${fit_sample_file}" "binning_${var}Nominal.txt" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
text = re.sub(
    r'^    binning: "(?!.*Sideband).*"$',  # signal sample only; sideband keeps its binning
    '    binning: "${GUNDAM_CONFIG_DIR}/binnings/' + sys.argv[2] + '"',
    text,
    flags=re.M,
)
# Sideband: fine reco binning for the fit, coarse nominal binning for toys/xsec.
text = re.sub(
    r'binnings/binning_SidebandVisE(Nominal)?\.txt',
    'binnings/binning_SidebandVisENominal.txt',
    text,
)
path.write_text(text)
PY
}

# NO_SIDEBAND=1 turns the NC pi0 sideband sample off for this run only (for
# with/without-sideband comparisons). It is switched back on by the EXIT trap.
set_sideband_enabled() {
  python3 - "${fit_sample_file}" "$1" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
text, n = re.subn(
    r'(- name: "sideband_NCpi0"\s*\n\s*isEnabled: )(true|false)',
    r'\g<1>' + sys.argv[2],
    path.read_text(),
)
if n == 0:
    sys.exit("NO_SIDEBAND: no sideband_NCpi0 sample found in " + str(path))
path.write_text(text)
PY
}

cleanup() {
  restore_default_binning
  if [[ "${NO_SIDEBAND:-}" == "1" ]]; then
    set_sideband_enabled true
  fi
}

trap cleanup EXIT

restore_default_binning
if [[ "${NO_SIDEBAND:-}" == "1" ]]; then
  echo "NO_SIDEBAND=1: fitting without the NC pi0 sideband sample"
  set_sideband_enabled false
fi

fitter_args=(
  -c "${fit_config}"
  -o "${output_dir}/fit_${tag}.root"
)

if [[ "${tag}" == "asimov" || "${ASIMOV:-}" == "1" ]]; then
  fitter_args+=(--asimov)
fi

gundamFitter "${fitter_args[@]}"

set_nominal_binning

gundamToyGenerator \
  -c "${toy_config}" \
  -f "${output_dir}/fit_${tag}.root" \
  -n "${n_toys}" \
  --use-bf \
  --save-data \
  -o "${output_dir}/toy_bf_${tag}.root"

gundamToyGenerator \
  -c "${toy_config}" \
  -f "${output_dir}/fit_${tag}.root" \
  -n "${n_toys}" \
  --use-prefit \
  --save-data \
  -o "${output_dir}/toy_pf_${tag}.root"

gundamCalcXsec \
  -c "${xsec_config}" \
  -f "${output_dir}/fit_${tag}.root" \
  -o "${output_dir}/xsec_${tag}.root" \
  -n "${n_toys}" \
  --use-bf-as-xsec

python3 "${script_dir}/makePrettyPlots.py" --tag "${tag}" --pot "${pot}" \
  --data-dir "${output_dir}" --output-dir "${plot_dir}"
