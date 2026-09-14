#!/usr/bin/env bash

set -euo pipefail

tag="${1:-asimov}"

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
export GUNDAM_CONFIG_DIR="${script_dir}"
fit_sample_file="${script_dir}/Configs_FitSampleSet/fitSamples_EleE.yaml"
output_dir="${script_dir}/data/${tag}"
plot_dir="${script_dir}/plots/${tag}"

mkdir -p "${output_dir}" "${plot_dir}"

if [[ ! -f "${fit_sample_file}" ]]; then
  echo "Missing fit sample file: ${fit_sample_file}" >&2
  exit 1
fi

restore_default_binning() {
  python3 - "${fit_sample_file}" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
text = re.sub(
    r'^    binning: ".*"$',
    '    binning: "${GUNDAM_CONFIG_DIR}/binnings/binning_EleE.txt"',
    text,
    flags=re.M,
)
path.write_text(text)
PY
}

set_nominal_binning() {
  python3 - "${fit_sample_file}" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
text = re.sub(
    r'^    binning: ".*"$',
    '    binning: "${GUNDAM_CONFIG_DIR}/binnings/binning_EleENominal.txt"',
    text,
    flags=re.M,
)
path.write_text(text)
PY
}

cleanup() {
  restore_default_binning
}

trap cleanup EXIT

restore_default_binning

fitter_args=(
  -c "${script_dir}/config_Fitter_FakeData_EleE.yaml"
  -o "${output_dir}/fit_${tag}.root"
)

if [[ "${tag}" == "asimov" ]]; then
  fitter_args+=(--asimov)
fi

gundamFitter "${fitter_args[@]}"

set_nominal_binning

gundamToyGenerator \
  -c "${script_dir}/config_Toy_FakeData_EleE.yaml" \
  -f "${output_dir}/fit_${tag}.root" \
  -n 1000 \
  --use-bf \
  --save-data \
  -o "${output_dir}/toy_bf_${tag}.root"

gundamToyGenerator \
  -c "${script_dir}/config_Toy_FakeData_EleE.yaml" \
  -f "${output_dir}/fit_${tag}.root" \
  -n 1000 \
  --use-prefit \
  --save-data \
  -o "${output_dir}/toy_pf_${tag}.root"

python3 "${script_dir}/makePrettyPlots.py" --tag "${tag}"