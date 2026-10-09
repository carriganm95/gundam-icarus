#!/bin/bash
#
# Grid job wrapper for gundamFitter (GUNDAM 1.9). Submitted by submitJobs.sh.
#
# Every job runs the same fit:
#   - GUNDAM-1.9.X/Install and configsNue/ come from the tarball, unpacked at
#     INPUT_TAR_DIR_LOCAL (cvmfs, read-only).
#   - The fitter config is shipped with -f and lands in CONDOR_DIR_INPUT.
#   - The input ROOT file(s) live on /pnfs and are copied here with ifdh, or
#     read in place over xrootd if inputMode=xrootd. Any path in the configs
#     ending in the same file name is rewritten to point at them.
#
# Environment (passed with -e by submitJobs.sh):
#   outputDir   /pnfs directory for outputs/ and logs/
#   configFile  config file name (shipped with -f)
#   inputFile   /pnfs path of the input file; comma-separate several files
#   inputMode   "ifdh" (default) or "xrootd"
#   fitterArgs  optional extra gundamFitter arguments; @JOB@ becomes the job number

jobNum=${PROCESS:-${CONDOR_PROC_ID:-}}
# Calculate subdirectory based on job number (groups of 100)
subdir=$(printf "%02d" $((${jobNum:-0} / 100)))
log_dir="${outputDir}/logs/${subdir}"

echo "PROCESS=$PROCESS  CONDOR_PROC_ID=$CONDOR_PROC_ID  jobNum=$jobNum"

errFile="${CONDOR_DIR_INPUT}/job_${jobNum}.err"
outFile="${CONDOR_DIR_INPUT}/job_${jobNum}.out"

# Redirect stdout and stderr to log files in subdirectories
exec > >(tee -a "$outFile") 2> >(tee -a "$errFile" >&2)

# Cleanup function to transfer logs before exit
cleanup_and_exit() {
    local exit_code=$1
    echo "Cleanup: transferring log files before exit with code $exit_code"

    # Check if log directory exists, create if it doesn't
    if ! ifdh ls "$log_dir" >/dev/null 2>&1; then
        echo "Creating log directory: $log_dir"
        ifdh mkdir_p "$log_dir" 2>/dev/null || true
    fi

    echo "Transferring log files to $log_dir"

    # Attempt to copy logs, but don't fail if it doesn't work
    if [[ -f "$errFile" ]]; then
        echo "Transferring err file to $log_dir"
        ifdh cp "$errFile" "$log_dir/" 2>/dev/null || echo "Warning: Failed to transfer $errFile"
    fi

    if [[ -f "$outFile" ]]; then
        echo "Transferring out file to $log_dir"
        ifdh cp "$outFile" "$log_dir/" 2>/dev/null || echo "Warning: Failed to transfer $outFile"
    fi

    exit $exit_code
}


required_vars=(outputDir jobNum CONDOR_DIR_INPUT INPUT_TAR_DIR_LOCAL configFile inputFile)
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: environment variable $var is unset or empty" >&2
        cleanup_and_exit 10
    fi
done

echo "Starting runJob.sh on grid worker node $HOSTNAME"
echo "outputDir is set to: $outputDir"
echo "logDir is set to: $log_dir"
echo "echo-ing JOBSUBJOBID: $JOBSUBJOBID"
echo "the job number is: $jobNum"

echo "doing general software setup"

# Run setup commands in the current shell using curly braces.
# Same spack packages GUNDAM-1.9.X was built against (see setup_spack.sh).
if ! {

    # --- Spack ---
    . /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
    spack load /hljl7gy           # root@6.28.12%gcc@12.2.0
    spack load gcc@12.2.0
    spack load /b2vqwgc           # ifdhc@2.8.0
    spack load xrootd@5.6.9%gcc@12.2.0

}; then
    echo "Issue setting up packages"
    cleanup_and_exit 15
fi

# --- GUNDAM from the tarball ---
echo "Tarball contents in ${INPUT_TAR_DIR_LOCAL}:"
ls -alh "${INPUT_TAR_DIR_LOCAL}"
INSTALL_DIR="${INPUT_TAR_DIR_LOCAL}/Install"
if [[ ! -x "${INSTALL_DIR}/gundam/bin/gundamFitter" ]]; then
    echo "ERROR: gundamFitter not found in ${INSTALL_DIR}/gundam/bin" >&2
    cleanup_and_exit 16
fi
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:${INSTALL_DIR}/yaml-cpp/:${INSTALL_DIR}/json/
export LD_LIBRARY_PATH=${INSTALL_DIR}/gundam/lib:${INSTALL_DIR}/yaml-cpp/lib64:$(root-config --libdir):$LD_LIBRARY_PATH
export PATH=${INSTALL_DIR}/gundam/bin:$PATH

echo "software setup complete"
echo "Using gundamFitter: $(which gundamFitter)"

cd "$CONDOR_DIR_INPUT"
work_dir="${PWD}"
echo "Working in ${work_dir}"
ls -ltrha

# --- Configs ---
# Copy configsNue out of the read-only tarball so input paths in any included
# config can be rewritten, then use the shipped top-level config.
if [[ ! -d "${INPUT_TAR_DIR_LOCAL}/configsNue" ]]; then
    echo "ERROR: configsNue not found in ${INPUT_TAR_DIR_LOCAL}" >&2
    cleanup_and_exit 20
fi
cp -r "${INPUT_TAR_DIR_LOCAL}/configsNue" "${work_dir}/configsNue" || cleanup_and_exit 21
# The copy keeps the tarball's read-only permissions
chmod -R u+w "${work_dir}/configsNue"
export GUNDAM_CONFIG_DIR="${work_dir}/configsNue"

configPath="${CONDOR_DIR_INPUT}/$(basename "${configFile}")"
if [[ ! -f "${configPath}" ]]; then
    echo "ERROR: Config ${configPath} not found" >&2
    cleanup_and_exit 22
fi

# Rewrite any quoted path ending in /<file name> to $2, in the shipped config
# and in configsNue
rewrite_input_path() {
    local base=$1 newPath=$2
    echo "Rewriting config paths ending in ${base} -> ${newPath}"
    if ! { echo "${configPath}"; grep -rl --include='*.yaml' --include='*.yml' --include='*.json' "/${base}\"" "${GUNDAM_CONFIG_DIR}"; } \
        | xargs -r sed -i "s|\"[^\"]*/${base}\"|\"${newPath}\"|g"; then
        echo "ERROR: failed to rewrite config paths for ${base}" >&2
        cleanup_and_exit 31
    fi
}

# --- Input file(s) from /pnfs ---
input_dir="${work_dir}/inputs"
mkdir -p "${input_dir}"
IFS=',' read -r -a inputList <<< "${inputFile}"
for f in "${inputList[@]}"; do
    [[ -z "$f" ]] && continue
    base=$(basename "$f")
    if [[ "${inputMode:-ifdh}" == "xrootd" ]]; then
        if [[ "$f" == /pnfs/* ]]; then
            url="root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/${f#/pnfs/}"
        else
            url="$f"
        fi
        echo "Reading input ${f} over xrootd"
        rewrite_input_path "${base}" "${url}"
    else
        echo "Copying input ${f}"
        if ! ifdh cp "$f" "${input_dir}/${base}"; then
            echo "ERROR: ifdh cp failed for input ${f}" >&2
            cleanup_and_exit 30
        fi
        rewrite_input_path "${base}" "${input_dir}/${base}"
    fi
done
ls -alh "${input_dir}"

# Anything still pointing at /exp will not be readable on the grid
if grep -n "\"/exp/[^\"]*\.root\"" "${configPath}" | grep -v "^[0-9]*:\s*#"; then
    echo "WARNING: ${configFile} still references /exp input files (lines above); add them to inputFile"
fi

# --- Run ---
outputFile="$(basename "${configFile}" .yaml)_${jobNum}.root"
extraArgs=${fitterArgs//@JOB@/${jobNum}}

echo "Executing: gundamFitter -c ${configPath} -o ${work_dir}/${outputFile} ${extraArgs}"
# extraArgs is intentionally unquoted so it splits into separate arguments
if ! gundamFitter -c "${configPath}" -o "${work_dir}/${outputFile}" ${extraArgs}; then
    echo "ERROR: gundamFitter failed" >&2
    cleanup_and_exit 40
fi

echo "Job completed checking contents of working directory"
ls -ltrha "${work_dir}"

if [[ ! -f "${work_dir}/${outputFile}" ]]; then
    echo "ERROR: Expected output ${work_dir}/${outputFile} not produced" >&2
    cleanup_and_exit 45
fi

echo "now trying to ifdh cp output file over to ${outputDir}"

target_dir="${outputDir}/outputs/${subdir}"

# Check if target directory exists, create if it doesn't
if ! ifdh ls "$target_dir" >/dev/null 2>&1; then
    echo "Creating output directory: $target_dir"
    ifdh mkdir_p "$target_dir" || { echo "ERROR: could not create $target_dir" >&2; cleanup_and_exit 50; }
fi

if ! ifdh cp "${work_dir}/${outputFile}" "$target_dir/${outputFile}"; then
    echo "ERROR: ifdh cp failed for ${outputFile} -> $target_dir/" >&2
    cleanup_and_exit 60
fi

echo "Job completed successfully"
cleanup_and_exit 0
#End of script
