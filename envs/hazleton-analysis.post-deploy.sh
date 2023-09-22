#!/usr/bin/env bash

# ---- Run this script once, right after the deployment of the homonymous conda environment
# ---- Usage: ./grups-hazelton-benchmark.post-deploy.sh



eval "$(conda shell.bash hook)"
conda activate $(basename "${0%.post-deploy.sh}")

R -e 'IRkernel::installspec()'

#R --slave -e 'devtools::install_github("git@github.com:MaelLefeuvre/grups.plots.git", upgrade="always")'
R --slave -e 'devtools::install("workflow/scripts/grups.plots", upgrade="always")'

