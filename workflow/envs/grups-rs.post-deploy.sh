#!/usr/bin/env bash

BUILD_DIR="workflow/scripts/grups/target"
[[ -d ${BUILD_DIR} ]] && rm -rf ${BUILD_DIR}


# ---- Install GRUPS-rs
export CARGO_INCREMENTAL=0
export CC="x86_64-conda-linux-gnu-gcc"
CMAKE_C_FLAGS="-I./ -L./ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" \
cargo install --path $(pwd)/workflow/scripts/grups-rs --root $CONDA_PREFIX

# ---- Install grups.plots
R --slave -e 'devtools::install_github("git@github.com:MaelLefeuvre/grups.plots.git", upgrade="always")'
