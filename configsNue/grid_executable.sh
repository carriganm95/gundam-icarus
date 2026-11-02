#!/bin/bash
set -e
set -o pipefail 

OutputDest=""
gundamTar=""
configTar=""
FitterConfigPath=""
OFFSET=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dest )              OutputDest=$2; shift ;;
    --gundamTar )         gundamTar=$2; shift ;;
    --configTar )         configTar=$2; shift ;;
    --FitterConfigPath )  FitterConfigPath=$2; shift ;;
    --offset )            OFFSET=$2; shift ;;
    * ) echo "Unknown option $1"; exit 1 ;;
  esac
  shift
done

nProcess=$((PROCESS + OFFSET + 1))   # seeds 1..N
echo "@@ Seed: $nProcess"

gridWD=$(pwd)

# --- Spack ---
. /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load /hljl7gy           # ROOT  -- VERIFY this hash works for icarus
spack load gcc@12.2.0
spack load /b2vqwgc           # ifdhc -- VERIFY
spack load xrootd@5.6.9%gcc@12.2.0

# --- GUNDAM tarball ---
export GUNDAM_DIR=$gridWD
export INSTALL_DIR=$GUNDAM_DIR/Install
export GUNDAM_CONFIG_DIR=$GUNDAM_DIR/work

ifdh cp $gundamTar $GUNDAM_DIR/local_gundam.tar
mkdir -p $INSTALL_DIR
tar xf local_gundam.tar -C $INSTALL_DIR

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$INSTALL_DIR/yaml-cpp/:$INSTALL_DIR/json/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_DIR/yaml-cpp/lib64/
export LIBRARY_PATH=$LIBRARY_PATH:$INSTALL_DIR/yaml-cpp/lib64/
export PATH=$INSTALL_DIR/gundam/bin:$PATH
export LD_LIBRARY_PATH=$INSTALL_DIR/gundam/lib:$LD_LIBRARY_PATH

# --- Config tarball ---
cd $gridWD
ifdh cp $configTar $GUNDAM_DIR/local_config.tar
mkdir -p $GUNDAM_CONFIG_DIR
tar xf local_config.tar -C $GUNDAM_CONFIG_DIR

# --- Stage input ROOT file to worker local scratch ---
INPUT_DIR=$gridWD/inputs
mkdir -p $INPUT_DIR
INPUT_FILE_PNFS=/pnfs/icarus/scratch/users/rvizarr/toy_ensemble/inputs/icarus_numi_numu_mc_onbeam_offbeam_syst_gundam_withPCA.root
INPUT_FILE_LOCAL=$INPUT_DIR/icarus_numi_numu_mc_onbeam_offbeam_syst_gundam_withPCA.root

echo "@@ Copying input ROOT file"
ifdh cp $INPUT_FILE_PNFS $INPUT_FILE_LOCAL
ls -lh $INPUT_FILE_LOCAL

# --- Rewrite paths in dataSetList to local file ---
EXP_PATH="/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_gundam_withPCA.root"
sed -i "s|${EXP_PATH}|${INPUT_FILE_LOCAL}|g" \
    $GUNDAM_CONFIG_DIR/Configs_DataSetList/dataSetListConfig_Asimov.yaml

echo "@@ Rewrote input path in dataSetListConfig_Asimov.yaml:"
grep -n "filePathList\|\.root" $GUNDAM_CONFIG_DIR/Configs_DataSetList/dataSetListConfig_Asimov.yaml

# --- Run ---
FitterOutput=output_Fitter_Toy_${nProcess}.root
FitterConfigFile=$GUNDAM_CONFIG_DIR/$FitterConfigPath

echo "@@ Running gundamFitter on $FitterConfigFile"
gundamFitter -c $FitterConfigFile -o $gridWD/$FitterOutput \
    -a --toy $nProcess 2>&1 | tee log.log

# --- Copy output back with retries ---
ifdh mkdir_p $OutputDest
for trial in $(seq 1 10); do
  if ifdh cp $gridWD/$FitterOutput $OutputDest/$FitterOutput; then
    echo "@@ Copy succeeded"
    break
  fi
  echo "@@ Copy failed, trial $trial; retrying in 10s"
  sleep 10
done
