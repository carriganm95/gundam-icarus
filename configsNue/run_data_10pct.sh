#!/usr/bin/env bash
# Full pipeline (fit -> toys -> xsec -> plots) on the 10% real-data sample.
#
# Fit/toys use fakeStudyConfigs/config_{Fitter,Toy}_FakeData_<var>_10pct.yaml
# (<var> = GUNDAM_VAR: EleE by default, or EleAngle)
# (events/data/selected, MC scaled to 2.36124e19 POT via
# dataSetListConfig_nueCCInclusiveData.yaml). calcXsec uses
# config_CalcXSec_Data_EleE_10pct.yaml, whose truth MC
# (dataSetListConfig_TrueEleE_dev.yaml) and flux (fitSamples_TrueEleE_10pct.yaml)
# are both normalised to 2.36124e19 POT.
#
# Usage: ./run_data_10pct.sh                # ROOT files -> data/<var>/data_10pct, plots -> plots/<var>/data_10pct
#        ./run_data_10pct.sh <dir>          # ROOT files -> <dir>, plots -> <dir>/plots
#   NTOYS=10 ./run_data_10pct.sh <dir>      # quick pipeline test

set -euo pipefail

script_dir="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
output="${1:-}"
pot="2.36124e19"

# GUNDAM_VAR=EleAngle runs the opening-angle analysis instead of EleE.
var="${GUNDAM_VAR:-EleE}"
export GUNDAM_VAR="${var}"
export XSEC_CONFIG="${script_dir}/config_CalcXSec_Data_${var}_10pct.yaml"

# Same data/<var>/<name> + plots/<var>/<name> layout as the fake-data studies.
name="data_10pct"
if [[ -z "${output}" ]]; then
  output="${script_dir}/data/${var}/${name}"
  export PLOT_DIR="${script_dir}/plots/${var}/${name}"
fi

exec bash "${script_dir}/run_fake_study.sh" 10pct "${pot}" "${output}"
