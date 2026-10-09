export GUNDAM_DIR=$PWD

export INSTALL_DIR="$GUNDAM_DIR/Install/"
export BUILD_DIR="$GUNDAM_DIR/Build/"
export REPO_DIR="$GUNDAM_DIR/Repositories/"

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:${INSTALL_DIR}/yaml-cpp/
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:${INSTALL_DIR}/json/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INSTALL_DIR}/yaml-cpp/lib64/
export LIBRARY_PATH=$LIBRARY_PATH:${INSTALL_DIR}/yaml-cpp/lib64/

# Finalize setup
export PATH="$INSTALL_DIR/gundam/bin:$PATH"
export LD_LIBRARY_PATH="$INSTALL_DIR/gundam/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH=$(root-config --libdir):$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:/exp/icarus/data/users/micarrig/.local/lib/python3.9/site-packages/

