#!/usr/bin/env bash

set -euo pipefail

BUILD_DIR="workflow/scripts/grups-rs/target"
[[ -d ${BUILD_DIR} ]] && rm -rf ${BUILD_DIR}


# ---- Install GRUPS-rs
export CARGO_INCREMENTAL=0
export CC="x86_64-conda-linux-gnu-gcc"
CMAKE_C_FLAGS="-I./ -L./ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" \
cargo install --path $(pwd)/workflow/scripts/grups-rs --root $CONDA_PREFIX

# ---- Install grups.plots
R --slave -e 'devtools::install_github("MaelLefeuvre/grups.plots", upgrade="always")'
