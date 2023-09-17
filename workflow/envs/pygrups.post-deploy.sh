#!/usr/bin/env bash

# ---- Clone repository
REPO="https://github.com/sameoldmike/grups.git"

cd "${CONDA_PREFIX}/bin"
git clone $REPO && cd grups

# ---- Convert scripts to executable format.
add_shebang(){
    sed -i '1i #!/usr/bin/env '"$1"'' $2
    chmod +x $2
}

add_shebang python2 pedigree_sims.py
add_shebang python2 PWD_from_stdin.py
add_shebang bash pedigree_sims_concat.sh
add_shebang Rscript plot_pedigree_sims.R

# ---- symlink
ln -srt ${CONDA_PREFIX}/bin *.py *.R *.sh
