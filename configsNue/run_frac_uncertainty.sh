#!/usr/bin/env bash
# Run fitter + calcXsec + toy generator for each systematic subset (individual
# and N-1), then generate fractional uncertainty comparison plots.

set -euo pipefail

# Parse flags: optional --diagonal, --zero-flux-norm, --output <dir>, positional POT.
# A custom output dir can also be given via the OUTPUT env var.
diagonal=false
zero_flux_norm=false
no_mc_throw=false
decompose=false
split_dist=false
pot="2.5e20"
output="${OUTPUT:-}"
# --split <set> isolates each individual parameter of one set (detector|genie|flux)
# in turn, instead of the group variants, via enableOnlyParameters.
# --split-dist additionally throws toys per variant so the split plots also show
# the pre/post-fit reconstructed-distribution breakdown (slower: +2 toys/variant).
# --decompose additionally splits the stat-only line into Template/Stat/EffStat/
# FluxNorm by re-running calcXsec on the template_only fit with each throw toggled.
split_set=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --diagonal)                   diagonal=true ;;
        --zero-flux-norm|--noflux)    zero_flux_norm=true ;;
        --no-mc-throw|--nomcthrow)    no_mc_throw=true ;;
        --decompose)                  decompose=true ;;
        --split-dist)                 split_dist=true ;;
        --output)                     output="$2"; shift ;;
        --output=*)                   output="${1#*=}" ;;
        --split)                      split_set="$2"; shift ;;
        --split=*)                    split_set="${1#*=}" ;;
        *)                            pot="$1" ;;
    esac
    shift
done

# Catch a bare word (e.g. a set name) that fell through to the positional POT,
# which otherwise reaches matplotlib as a non-numeric POT and crashes the plots.
if ! [[ "${pot}" =~ ^[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
    echo "ERROR: '${pot}' is not a numeric POT value." >&2
    echo "       Flags need the leading '--'. Did you mean:  --split ${pot}  ?" >&2
    exit 1
fi
if $split_dist && [[ -z "${split_set}" ]]; then
    echo "ERROR: --split-dist requires --split <detector|genie|flux>" >&2
    exit 1
fi

if $decompose; then
    [[ -n "${split_set}" ]] && { echo "ERROR: --decompose cannot be combined with --split" >&2; exit 1; }
    $zero_flux_norm && { echo "ERROR: --decompose manages the flux norm itself; drop --zero-flux-norm" >&2; exit 1; }
    $no_mc_throw && { echo "ERROR: --decompose manages the MC throw itself; drop --no-mc-throw" >&2; exit 1; }
fi

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
export GUNDAM_CONFIG_DIR="${script_dir}"

# Analysis variable: EleE (electron energy, default) or EleAngle (cos theta_e).
# Inherited by the plot scripts through the environment.
var="${GUNDAM_VAR:-EleE}"
export GUNDAM_VAR="${var}"

fitter_base="${script_dir}/fakeStudyConfigs/config_Fitter_FakeData_${var}_asimov.yaml"
xsec_base="${script_dir}/config_CalcXSec_FakeData_${var}.yaml"
toy_base="${script_dir}/fakeStudyConfigs/config_Toy_FakeData_${var}_asimov.yaml"
fit_sample_file="${script_dir}/Configs_FitSampleSet/fitSamples_${var}.yaml"
xs_sample_file="${script_dir}/Configs_FitSampleSet/fitSamples_True${var}.yaml"
# Systematic-subset parameterSetLists: EleE's live in frac_uncertainty/, other
# variables' (same lists, different template set) in frac_uncertainty/<var>/.
frac_rel="frac_uncertainty"
if [[ "${var}" != "EleE" ]]; then
    frac_rel="frac_uncertainty/${var}"
fi

# --split mode: resolve the chosen set to its parameterSet file, dial definition
# (for the parameter names), and the existing group-variant used as the full-set
# reference curve.
split_params=()
split_is_norm=false
if [[ -n "${split_set}" ]]; then
    case "${split_set}" in
        detector)   set_parset="Configs_ParameterSet/Detector/parameterSet_Detector.yaml"
                    set_dialdef="Configs_ParameterSet/Detector/DetectorDialsDef.yaml"
                    set_ref="detector_only" ;;
        genie|xsec) set_parset="Configs_ParameterSet/GENIE/parameterSet_GENIE.yaml"
                    set_dialdef="Configs_ParameterSet/GENIE/GENIEDialsDef.yaml"
                    set_ref="xsec_only" ;;
        flux)       set_parset="Configs_ParameterSet/Flux/parameterSet_Beamline.yaml"
                    set_dialdef="Configs_ParameterSet/Flux/NuMIBeamlineDialsDef.yaml"
                    set_ref="flux_only" ;;
        norm|fluxnorm|normalisation)
                    # Isolate the normaliseParameters in the xsec sample def
                    # (flux normalisation, targets). Handled by a dedicated flow below.
                    split_is_norm=true ;;
        *) echo "ERROR: unknown --split set '${split_set}' (expected detector|genie|flux|norm)" >&2; exit 1 ;;
    esac
    if ! $split_is_norm; then
        # Parameter names, in dial-def order.
        mapfile -t split_params < <(grep -oE 'parameterName:[[:space:]]*"[^"]+"' \
            "${script_dir}/${set_dialdef}" | sed -E 's/.*"([^"]+)".*/\1/')
        if [[ ${#split_params[@]} -eq 0 ]]; then
            echo "ERROR: no parameterName entries found in ${set_dialdef}" >&2; exit 1
        fi
        echo "Split mode: isolating ${#split_params[@]} '${split_set}' parameters."
    fi
fi

if $diagonal; then
    # Use true energy as the "reco" variable — response matrix is diagonal by construction.
    binning_ele="binning_True${var}.txt"
    binning_nominal="binning_True${var}.txt"
    default_name="frac_uncertainty_asimov_diagonal"
else
    binning_ele="binning_${var}.txt"
    binning_nominal="binning_${var}Nominal.txt"
    default_name="frac_uncertainty_asimov"
fi

# --zero-flux-norm: set the "Constant flux" normalisation stdDev to 0 in the xsec
# sample def, removing the 9.13% flux-normalisation uncertainty baked into every
# cross-section bin by gundamCalcXsec. Only affects calcXsec output, not the fit.
if $zero_flux_norm; then
    default_name="${default_name}_noflux"
fi

# --no-mc-throw: set enableEventMcThrow:false in the calcXsec config, removing the
# MC-statistical noise of the efficiency estimate from the cross-section error.
# Lets you see the physical systematic separate from the (uncorrelated, binning-
# sensitive) finite-MC term. Only affects calcXsec output.
if $no_mc_throw; then
    default_name="${default_name}_nomcthrow"
fi

# Split runs get their own data/plots tree, e.g. split_detector_asimov[_diagonal].
if [[ -n "${split_set}" ]]; then
    default_name="split_${split_set}_${default_name#frac_uncertainty_}"
fi
if $decompose; then
    default_name="decompose_${default_name#frac_uncertainty_}"
fi

# Output location: a custom --output/OUTPUT dir is used verbatim (variant subdirs
# go directly under it, plots under <dir>/plots); otherwise the default
# data/<var>/<name> and plots/<var>/<name> layout keyed off the diagonal/noflux flags.
if [[ -n "${output}" ]]; then
    output_base="${output}"
    plot_dir="${output}/plots"
else
    output_base="${script_dir}/data/${var}/${default_name}"
    plot_dir="${script_dir}/plots/${var}/${default_name}"
fi
# PLOT_DIR overrides where the plots go (same as run_fake_study.sh).
plot_dir="${PLOT_DIR:-${plot_dir}}"

mkdir -p "${plot_dir}"

# fitSamples_<var>.yaml is modified in-place to switch binning.
# In normal mode: fitter uses the fine reco binning (binning_<var>.txt), toy uses
# binning_<var>Nominal.txt (truth edges, matches calcXsec).
# In diagonal mode: both use binning_True<var>.txt (true = reco → diagonal matrix).
set_binning() {
  python3 - "${fit_sample_file}" "$1" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
bfile = sys.argv[2]
# Only touch the signal sample: the sideband keeps its own binning file.
path.write_text(re.sub(r'^(\s+binning: ")(?!.*Sideband).*(")', r'\1${GUNDAM_CONFIG_DIR}/binnings/' + bfile + r'\2',
                path.read_text(), flags=re.M))
PY
}

# Sideband: fine reco binning (binning_SidebandVisE.txt) for the fit, coarse
# binning_SidebandVisENominal.txt for toys/xsec, same scheme as the signal.
set_sideband_binning() {
  python3 - "${fit_sample_file}" "$1" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
path.write_text(re.sub(r'binnings/binning_SidebandVisE(Nominal)?\.txt',
                       'binnings/binning_SidebandVisE' + sys.argv[2] + '.txt', path.read_text()))
PY
}

set_ele_binning()     { set_binning "${binning_ele}"; set_sideband_binning ""; }
set_nominal_binning() { set_binning "${binning_nominal}"; set_sideband_binning "Nominal"; }

# Individual-contribution variants always run. The N-1 variants are skipped in
# diagonal mode: with a diagonal response the free per-bin templates are exactly
# degenerate with any reweighting systematic, so those fits are ill-posed (Migrad
# diverges to NaN). The individual variants are all that's needed for the
# individual-contributions plots.
if $split_is_norm; then
    VARIANTS=()   # the normalisation split uses its own flow (below) and exits early
elif [[ -n "${split_set}" ]]; then
    # Baseline (stat), the full-set reference curve, then each isolated parameter.
    VARIANTS=(template_only "${set_ref}" "${split_params[@]}")
elif $decompose; then
    # Group individual variants only; the stat decomposition runs after the loop.
    VARIANTS=(all template_only xsec_only flux_only detector_only)
else
    VARIANTS=(all template_only xsec_only flux_only detector_only)
    # G4 group: only where its parameter-set lists exist (EleE for now).
    has_g4=false
    [[ -f "${script_dir}/${frac_rel}/parameterSetList_g4_only.yaml" ]] && has_g4=true
    $has_g4 && VARIANTS+=(g4_only)
    if ! $diagonal; then
        VARIANTS+=(no_xsec no_flux no_detector)
        $has_g4 && VARIANTS+=(no_g4)
    fi
fi

tmp_fitter=""
tmp_xsec=""
tmp_toy=""
xs_backup=""
xsec_cfg_backup=""
split_gen_files=()   # per-parameter parameterSet/parameterSetList files to clean up
cleanup() {
    [[ -n "${tmp_fitter}" ]] && rm -f "${tmp_fitter}"
    [[ -n "${tmp_xsec}" ]]   && rm -f "${tmp_xsec}"
    [[ -n "${tmp_toy}" ]]    && rm -f "${tmp_toy}"
    set_ele_binning  # restore fitter binning on exit
    # restore the original xsec sample def if we zeroed the flux norm; the -f file
    # check keeps this idempotent if cleanup runs more than once.
    [[ -n "${xs_backup}" && -f "${xs_backup}" ]] && mv -f "${xs_backup}" "${xs_sample_file}"
    # restore the original calcXsec config if we disabled the MC throw
    [[ -n "${xsec_cfg_backup}" && -f "${xsec_cfg_backup}" ]] && mv -f "${xsec_cfg_backup}" "${xsec_base}"
    if [[ ${#split_gen_files[@]} -gt 0 ]]; then rm -f "${split_gen_files[@]}"; fi
}
# EXIT alone also fires on a set -e abort, so this avoids the double cleanup that
# a combined EXIT+ERR trap causes (second pass can't find an already-moved backup).
trap cleanup EXIT

# In --split mode, isolate each parameter by generating a dial-definition that
# contains ONLY that parameter (the enableOnlyParameters config key is silently
# ignored by this GUNDAM build), a parameterSet pointing at it, and a
# parameterSetList (Template + that set). All live under frac_uncertainty/ so the
# existing sed path logic finds them, and are removed on exit.
if [[ -n "${split_set}" ]]; then
    for p in "${split_params[@]}"; do
        dd="${script_dir}/frac_uncertainty/_split_dialdef_${p}.yaml"
        ps="${script_dir}/frac_uncertainty/_split_parameterSet_${p}.yaml"
        pl="${script_dir}/${frac_rel}/parameterSetList_${p}.yaml"
        # dial-def with only parameter p (keeps GUNDAM from building the others)
        python3 - "${script_dir}/${set_dialdef}" "${p}" "${dd}" <<'PY'
import yaml, sys
data = yaml.safe_load(open(sys.argv[1]))
kept = [e for e in data if e.get("parameterName") == sys.argv[2]]
assert len(kept) == 1, f"expected exactly 1 dial-def entry for {sys.argv[2]}, got {len(kept)}"
yaml.safe_dump(kept, open(sys.argv[3], "w"), default_flow_style=False, sort_keys=False)
PY
        # parameterSet copy whose parameterDefinitions points at the single-dial def
        python3 - "${script_dir}/${set_parset}" "${p}" "${ps}" <<'PY'
import re, sys, pathlib
t = pathlib.Path(sys.argv[1]).read_text()
dd = '"${GUNDAM_CONFIG_DIR}/frac_uncertainty/_split_dialdef_' + sys.argv[2] + '.yaml"'
t, n = re.subn(r'(parameterDefinitions:\s*)"[^"]*"', r'\g<1>' + dd, t, count=1)
assert n == 1, "could not rewrite parameterDefinitions path"
pathlib.Path(sys.argv[3]).write_text(t)
PY
        printf -- '- "${GUNDAM_CONFIG_DIR}/Configs_ParameterSet/TemplateParameters/parameterSet_True%s.yaml"\n- "${GUNDAM_CONFIG_DIR}/frac_uncertainty/_split_parameterSet_%s.yaml"\n' "${var}" "${p}" > "${pl}"
        split_gen_files+=("${dd}" "${ps}" "${pl}")
    done
fi

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

if $no_mc_throw; then
    xsec_cfg_backup=$(mktemp --suffix=.yaml)
    cp "${xsec_base}" "${xsec_cfg_backup}"
    python3 - "${xsec_base}" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
text, n = re.subn(r'(enableEventMcThrow:\s*)true', r'\g<1>false', path.read_text(), count=1)
assert n == 1, "could not find 'enableEventMcThrow: true' to disable"
path.write_text(text)
print("Set enableEventMcThrow: false in", path.name)
PY
fi

# --split norm/fluxnorm: isolate each normaliseParameter of the xsec sample def
# (flux normalisation, number of targets). They affect only calcXsec, so a single
# template_only fit is reused and calcXsec is re-run with each normaliser toggled:
#   nonorm   -> all normalisation stdDevs set to 0 (baseline)
#   <name>   -> only that normaliser kept, the rest zeroed (isolated contribution)
#   allnorm  -> the full normalisation (reference)
if $split_is_norm; then
    $zero_flux_norm && { echo "ERROR: --zero-flux-norm conflicts with --split ${split_set}" >&2; exit 1; }
    mapfile -t norm_names < <(python3 - "${xs_sample_file}" <<'PY'
import re, sys, pathlib
t = pathlib.Path(sys.argv[1]).read_text()
m = re.search(r'normaliseParameterList:(.*)$', t, flags=re.S)
for nm in re.findall(r'name:\s*"([^"]+)"', m.group(1) if m else ""):
    print(nm)
PY
)
    [[ ${#norm_names[@]} -eq 0 ]] && { echo "ERROR: no normaliseParameterList entries in ${xs_sample_file}" >&2; exit 1; }
    echo "Normalisation split: isolating: ${norm_names[*]}"

    xs_backup=$(mktemp --suffix=.yaml)
    cp "${xs_sample_file}" "${xs_backup}"

    # One shared template_only fit (normaliseParameters do not enter the fit).
    fitroot="${output_base}/_fit/fit.root"
    mkdir -p "$(dirname "${fitroot}")"
    tmp_fitter=$(mktemp --suffix=.yaml)
    sed "s|Configs_ParameterSetList/parameterSetList_${var}|${frac_rel}/parameterSetList_template_only|" \
        "${fitter_base}" > "${tmp_fitter}"
    gundamFitter -c "${tmp_fitter}" -o "${fitroot}" --asimov
    rm -f "${tmp_fitter}"; tmp_fitter=""

    # Keep only the named normaliser (or __none__ / __all__); zero the rest.
    norm_keep() {
        cp "${xs_backup}" "${xs_sample_file}"
        python3 - "${xs_sample_file}" "$1" <<'PY'
import re, sys, pathlib
p = pathlib.Path(sys.argv[1]); keep = sys.argv[2]; t = p.read_text()
m = re.search(r'normaliseParameterList:(.*)$', t, flags=re.S)
for nm in re.findall(r'name:\s*"([^"]+)"', m.group(1) if m else ""):
    if keep == "__all__" or (keep != "__none__" and nm == keep):
        continue
    t = re.sub(r'(name:\s*"' + re.escape(nm) + r'".*?stdDev:\s*)[0-9.eE+\-]+',
               r'\g<1>0.0', t, count=1, flags=re.S)
p.write_text(t)
PY
    }
    norm_xsec() {  # $1 = output variant dir
        local tmpx; tmpx=$(mktemp --suffix=.yaml)
        sed "s|Configs_ParameterSetList/parameterSetList_${var}|${frac_rel}/parameterSetList_template_only|" \
            "${xsec_base}" > "${tmpx}"
        mkdir -p "${output_base}/$1"
        gundamCalcXsec -c "${tmpx}" -f "${fitroot}" -o "${output_base}/$1/xsec.root" \
            -n 1000 --use-bf-as-xsec || echo "WARNING: calcXsec failed for norm variant '$1'"
        rm -f "${tmpx}"
    }

    norm_keep "__none__"; norm_xsec "nonorm"      # stat baseline, no normalisation
    for nm in "${norm_names[@]}"; do
        norm_keep "${nm}"; norm_xsec "${nm// /_}" # isolated normaliser
    done
    norm_keep "__all__";  norm_xsec "allnorm"     # full normalisation reference

    cp "${xs_backup}" "${xs_sample_file}"

    python3 "${script_dir}/scripts/plot_split_systematic.py" \
        --data-dir "${output_base}" --output-dir "${plot_dir}" \
        --set "${split_set}" --baseline nonorm --reference allnorm --pot "${pot}"
    echo "Done. Plots saved to ${plot_dir}"
    exit 0
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

    sed "s|Configs_ParameterSetList/parameterSetList_${var}|${frac_rel}/parameterSetList_${variant}|" \
        "${fitter_base}" > "${tmp_fitter}"
    sed "s|Configs_ParameterSetList/parameterSetList_${var}|${frac_rel}/parameterSetList_${variant}|" \
        "${xsec_base}" > "${tmp_xsec}"
    sed "s|Configs_ParameterSetList/parameterSetList_${var}|${frac_rel}/parameterSetList_${variant}|" \
        "${toy_base}" > "${tmp_toy}"

    # Process each variant tolerantly: a failed fit (e.g. a degenerate variant
    # where Migrad diverges to NaN) is logged and skipped, never aborts the run.
    variant_ok=true

    # Asimov fits start at the minimum; the Simplex pre-step can leave a parameter
    # at NaN there (seen for RPA_CCQE), so a failed fit is retried once without it.
    if ! gundamFitter -c "${tmp_fitter}" -o "${out_dir}/fit.root" --asimov; then
        echo "WARNING: gundamFitter failed for '${variant}' -- retrying without Simplex"
        if ! gundamFitter -c "${tmp_fitter}" -o "${out_dir}/fit.root" --asimov \
               -O /fitterEngineConfig/minimizerConfig/enableSimplexBeforeMinimize=false; then
            echo "WARNING: gundamFitter failed for '${variant}' (also without Simplex) -- skipping this variant"
            variant_ok=false
        fi
    fi

    # Self-check: a per-dial split variant must end up with exactly ONE systematic
    # parameter active. If more are active the isolation silently failed (as the
    # ignored enableOnlyParameters did), so abort rather than produce bogus splits.
    if ${variant_ok} && [[ -n "${split_set}" ]] && [[ " ${split_params[*]} " == *" ${variant} "* ]]; then
        nsys=$(python3 -c "
import uproot
try:
    labs=uproot.open('${out_dir}/fit.root')['FitterEngine/postFit/Migrad/hessian/postfitCovariance_TH2D'].axis(0).labels()
    print(sum(1 for l in labs if 'template' not in l.lower()))
except Exception: print(-1)
" 2>/dev/null)
        if [[ "${nsys}" == "0" ]]; then
            # Parameter is in the dial-def but has no covariance prior (e.g. GENIE
            # CC1piTPi / LowQ2Suppression), so it can't be a fit parameter -> skip it.
            echo "  [isolation check] '${variant}': not in the covariance (0 fit parameters) -- skipping."
            rm -f "${tmp_fitter}" "${tmp_xsec}" "${tmp_toy}"; tmp_fitter=""; tmp_xsec=""; tmp_toy=""
            rm -rf "${out_dir}"
            set_ele_binning
            continue
        elif [[ "${nsys}" != "1" ]]; then
            echo "ERROR: split variant '${variant}' has ${nsys} systematic parameters active (expected 1)." >&2
            echo "       Isolation FAILED -- per-dial results would be invalid. Aborting." >&2
            exit 1
        fi
        echo "  [isolation check] '${variant}': 1 systematic parameter active (OK)"
    fi

    if ${variant_ok}; then
        if ! gundamCalcXsec -c "${tmp_xsec}" -f "${out_dir}/fit.root" \
                -o "${out_dir}/xsec.root" -n 1000 --use-bf-as-xsec; then
            echo "WARNING: gundamCalcXsec failed for '${variant}' -- skipping this variant"
            variant_ok=false
        fi
    fi

    # Toys give the reconstructed-distribution bands. Normal runs always make them;
    # --split skips them (xsec-only) unless --split-dist is set; --decompose never.
    if ${variant_ok} && ! $decompose && { [[ -z "${split_set}" ]] || $split_dist; }; then
        set_nominal_binning  # switch to the truth-edge nominal binning for the toys (matches calcXsec)
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

if $decompose; then
    echo "============================================"
    echo " Stat decomposition (calcXsec toggles on the template_only fit)"
    echo "============================================"
    tpl_fit="${output_base}/template_only/fit.root"
    # pristine backups; restored by the main cleanup trap on exit
    xsec_cfg_backup=$(mktemp --suffix=.yaml); cp "${xsec_base}" "${xsec_cfg_backup}"
    xs_backup=$(mktemp --suffix=.yaml);       cp "${xs_sample_file}" "${xs_backup}"

    # Set calcXsec throws ($1=statThrow, $2=mcThrow as true/false) and normaliser
    # keep ($3 = __none__ / __all__ / a name), always from the pristine configs.
    decomp_set() {
        cp "${xsec_cfg_backup}" "${xsec_base}"
        python3 - "${xsec_base}" "$1" "$2" <<'PY'
import re, sys, pathlib
p = pathlib.Path(sys.argv[1]); t = p.read_text()
t = re.sub(r'(enableStatThrowInToys:\s*)\w+', r'\g<1>' + sys.argv[2], t, count=1)
t = re.sub(r'(enableEventMcThrow:\s*)\w+',    r'\g<1>' + sys.argv[3], t, count=1)
p.write_text(t)
PY
        cp "${xs_backup}" "${xs_sample_file}"
        python3 - "${xs_sample_file}" "$3" <<'PY'
import re, sys, pathlib
p = pathlib.Path(sys.argv[1]); keep = sys.argv[2]; t = p.read_text()
m = re.search(r'normaliseParameterList:(.*)$', t, flags=re.S)
for nm in re.findall(r'name:\s*"([^"]+)"', m.group(1) if m else ""):
    if keep == "__all__" or (keep != "__none__" and nm == keep):
        continue
    t = re.sub(r'(name:\s*"' + re.escape(nm) + r'".*?stdDev:\s*)[0-9.eE+\-]+',
               r'\g<1>0.0', t, count=1, flags=re.S)
p.write_text(t)
PY
    }
    decomp_xsec() {  # $1 = output variant dir
        local tmpx; tmpx=$(mktemp --suffix=.yaml)
        sed "s|Configs_ParameterSetList/parameterSetList_${var}|${frac_rel}/parameterSetList_template_only|" \
            "${xsec_base}" > "${tmpx}"
        mkdir -p "${output_base}/$1"
        gundamCalcXsec -c "${tmpx}" -f "${tpl_fit}" -o "${output_base}/$1/xsec.root" \
            -n 1000 --use-bf-as-xsec || echo "WARNING: decomp calcXsec failed for '$1'"
        rm -f "${tmpx}"
    }

    #          statThrow mcThrow  normKeep             -> variant dir
    decomp_set false     false    "__none__";          decomp_xsec template_base
    decomp_set true      false    "__none__";          decomp_xsec statthrow
    decomp_set false     true     "__none__";          decomp_xsec mcthrow
    decomp_set false     false    "Constant flux";     decomp_xsec fluxnorm
    decomp_set false     false    "Number of targets"; decomp_xsec targets
fi

echo "============================================"
echo " Generating comparison plots"
echo "============================================"

if $decompose; then
    python3 "${script_dir}/scripts/plot_full_breakdown.py" \
        --data-dir "${output_base}" \
        --output-dir "${plot_dir}" \
        --pot "${pot}"
elif [[ -n "${split_set}" ]]; then
    dist_arg=(); $split_dist && dist_arg=(--with-dist)
    python3 "${script_dir}/scripts/plot_split_systematic.py" \
        --data-dir "${output_base}" \
        --output-dir "${plot_dir}" \
        --set "${split_set}" \
        --baseline template_only \
        --reference "${set_ref}" \
        --pot "${pot}" \
        "${dist_arg[@]}"
else
    python3 "${script_dir}/scripts/plot_frac_uncertainty.py" \
        --data-dir "${output_base}" \
        --output-dir "${plot_dir}" \
        --pot "${pot}"
fi

echo "Done. Plots saved to ${plot_dir}"
