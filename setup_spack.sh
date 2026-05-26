#!/bin/bash
# Load packages from spack

# spack env
. /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh

# root
# - hljl7gy root@6.28.12%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
spack load /hljl7gy

# gcc
spack load gcc@12.2.0

# ifdhc
# - b2vqwgc ifdhc@2.8.0%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
spack load /b2vqwgc

# xrootd
spack load xrootd@5.6.9%gcc@12.2.0

#Path to example configs for 2025 GUNDAM workshop
export GUNDAM_CONFIG_DIR="$PWD/configsNue"
