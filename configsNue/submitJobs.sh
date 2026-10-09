#! /bin/bash

# Location of this script (configsNue/) and the GUNDAM 1.9 build next to it
scriptDir="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
gundamDir="${scriptDir}/../GUNDAM-1.9.X"

sourceDir="/pnfs/icarus/scratch/users/micarrig/nueCCInclusive_gundam/asimov"
exe="${scriptDir}/runJob.sh"
export tarFile="gundam.tar.gz"
# The tar file must stay on local disk: jobsub reads it to upload it to cvmfs,
# and files on /pnfs cannot be read directly from the gpvms
tarDir="/exp/icarus/data/users/micarrig/gundam_grid"

# Input file(s) for the fit. If not already on /pnfs it is copied to sourceDir.
# Comma-separate several files.
export inputFile="/exp/icarus/data/users/micarrig/gundam_v30.root"

# "ifdh" copies the input to the worker; "xrootd" reads it in place (large files)
export inputMode="ifdh"

# Fitter config, relative to configsNue/. Shipped to every job.
export configFile="fakeStudyConfigs/config_Fitter_FakeData_EleE_asimovGrid.yaml"

# Optional extra gundamFitter arguments; @JOB@ becomes the job number,
# e.g. "--toy @JOB@" or "-s @JOB@"
export fitterArgs="--toy @JOB@ -s @JOB@"
export outputDir=$sourceDir
nJobs=1000
recopy=false

# Rebuild the tar file from GUNDAM-1.9.X/Install and configsNue/ and upload it
makeTar=true

if [ ! -d $outputDir ]; then
    mkdir -p $outputDir
fi
if [ ! -d $outputDir/logs ]; then
    mkdir -p $outputDir/logs
fi
if [ ! -d "$outputDir/outputs" ]; then
    mkdir -p "$outputDir/outputs"
fi

# Create subdirectories for outputs and logs based on number of jobs
for ((i=0; i<$((($nJobs + 99) / 100)); i++)); do
    subdir=$(printf "%02d" $i)
    if [ ! -d $outputDir/outputs/$subdir ]; then
        mkdir -p $outputDir/outputs/$subdir
    fi
    if [ ! -d $outputDir/logs/$subdir ]; then
        mkdir -p $outputDir/logs/$subdir
    fi
done

# Build the tar file: Install/ (GUNDAM 1.9) and configsNue/ at the top level
if [ -n "$tarFile" ] && [ "$makeTar" = true ]; then
    echo "Building tar file ${tarFile}"
    if [ ! -x "${gundamDir}/Install/gundam/bin/gundamFitter" ]; then
        echo "ERROR: no GUNDAM install found in ${gundamDir}/Install" >&2
        exit 1
    fi
    stageDir=$(mktemp -d)
    cp -r "${gundamDir}/Install" "${stageDir}/Install"
    rsync -a \
        --exclude data --exclude plots --exclude .venv --exclude notes \
        --exclude '.ipynb_checkpoints' --exclude '__pycache__' --exclude '*.ipynb' \
        "${scriptDir}/" "${stageDir}/configsNue/"
    mkdir -p "${tarDir}"
    tar czf "${tarDir}/${tarFile}" -C "${stageDir}" Install configsNue
    ls -lh "${tarDir}/${tarFile}"
    rm -rf "${stageDir}"
fi

# Keep a copy of the submitted tar file next to the outputs for reference
if [ -n "$tarFile" ] && [ -f "${tarDir}/${tarFile}" ] && \
   ( [ "$makeTar" = true ] || [ ! -f "${sourceDir}/${tarFile}" ] || [ "$recopy" = true ] ); then
    echo "Copying tar file to ${sourceDir}"
    cp "${tarDir}/${tarFile}" "${sourceDir}/${tarFile}"
fi

# Stage input file(s) on /pnfs; jobs copy or stream them from there
pnfsInputs=()
IFS=',' read -r -a inputList <<< "${inputFile}"
for f in "${inputList[@]}"; do
    [ -z "$f" ] && continue
    if [[ "$f" == /pnfs/* ]]; then
        pnfsInputs+=("$f")
        continue
    fi
    base=$(basename "$f")
    if [ ! -f "${sourceDir}/${base}" ] || [ "$recopy" = true ]; then
        echo "Copying input file ${f} to scratch area"
        cp "${f}" "${sourceDir}/${base}"
    fi
    pnfsInputs+=("${sourceDir}/${base}")
done
export inputFile=$(IFS=','; echo "${pnfsInputs[*]}")

# Always copy the config so jobs never run a stale version
configBase=$(basename "${configFile}")
if [ -n "$configFile" ]; then
    if [ ! -f "${scriptDir}/${configFile}" ]; then
        echo "ERROR: config ${scriptDir}/${configFile} not found" >&2
        exit 1
    fi
    echo "Copying config file to scratch area"
    cp "${scriptDir}/${configFile}" "${sourceDir}/${configBase}"
fi

#not clear to me whether the following two lines are needed
export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u) && htgettoken -a htvaultprod.fnal.gov -i icarus

echo "Submitting $nJobs jobs to process"
echo "Using executable $exe and tar file $tarFile"
echo "Config: $configFile, input: $inputFile ($inputMode)"
echo "Output will be stored in $outputDir"

#
#jobsub_cmd="jobsub_submit -G icarus -N ${nJobs} --expected-lifetime=18h --disk=25GB --memory=8000MB -e IFDH_CP_MAXRETRIES=4 -e IFDH_CP_UNLINK_ON_ERROR=2 --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=3600' -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/carriganm95/uproot\:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true)'"
# Note: --maxConcurrent (HTCondor late materialization) causes $PROCESS to be 0
# for every job. Omit it, or if needed re-add it and also pass
#   --lines 'environment = "PROCESS=$(Process)"'
# so the process number is stamped at materialization time.
jobsub_cmd="jobsub_submit -G icarus -N ${nJobs} --expected-lifetime=18h --disk=25GB --memory=8000MB -e IFDH_CP_MAXRETRIES=4 -e IFDH_CP_UNLINK_ON_ERROR=2 --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=3600' -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9\:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true)'"

if [ -n "$tarFile" ] && [ -f "${tarDir}/${tarFile}" ]; then
    jobsub_cmd="${jobsub_cmd} --tar_file_name dropbox://${tarDir}/${tarFile} --use-cvmfs-dropbox -e tarFile"
else
    echo "ERROR: tar file ${tarDir}/${tarFile} not found; set makeTar=true to build it" >&2
    exit 1
fi

jobsub_cmd="${jobsub_cmd} -e inputFile -e inputMode"

if [ -n "$configFile" ]; then
    jobsub_cmd="${jobsub_cmd} -e configFile -f ${sourceDir}/${configBase}"
fi

if [ -n "$fitterArgs" ]; then
    jobsub_cmd="${jobsub_cmd} -e fitterArgs"
fi

jobsub_cmd="${jobsub_cmd} -e outputDir file://${exe}"

echo "Executing jobsub command:"
echo "$jobsub_cmd"
# jobsub writes a repacked copy of the tar file to the current directory
cd "${tarDir}"
eval $jobsub_cmd
