## Usage

### Setup
1. Install [miniconda3](https://docs.conda.io/projects/miniconda/en/latest/)
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. (Optional) Install mamba
   ```bash
   conda install -n base -c conda-forge mamba
   ```

4. Create the base snakemake environment
   ```
   conda env create -f workflow/envs/snakemake-7.20.0.yml
   ```

### Running GRUPS-rs' benchmark

Optionally tweak the `config/config.yml` configuration file beforehand. 

Then, simply run the snakemake pipeline using the following command

```bash
conda activate snakemake-7.20.0
snakemake grups_rs_bench --cores 1 --use-conda --conda-frontend mamba --printshellcmds
```

This should:
1. Download and preprocess the 1000g-phase3 dataset in `data/1000g-phase3`
2. Download the HapMapII recombination map in `data/recombination-maps`
3. Download the AADR v52.2 1240K variant callser in `data/aadr-1240-dataset`
4. FSA-encode the 1000g-phase3 dataset (if `mode` param is set to `fst` or `fst-mmap`)
5. Apply GRUPS-rs 10 times on the following parameter space:
  - number of pairwise comparisons: (1, 3, 6, 10, 15, 21, 28, 36, 45, 55]
  - Average pairwise SNP overlap: [2247, 4494, 8988, 17975, 35950, 71899, 143797, 287593, 575185, 1150369]
6. Plot a 3D surface plot of the average runtime and maximum resident set size

#### Output:

- raw benchmark results should be located in `benchmarks`
- 3D-surface plots should be located in `results/bench-plots`


### Running PyGrups' benchmark
```bash
conda activate snakemake-7.20.0
snakemake pygrups_bench --cores 1 --use-conda --conda-frontend mamba --printshellcmds
```

#### raw benchmark results should be located in `benchmarks/pygrups`

### Compute the average depth
```Shell
BAM="original-data/1X-EUR-pedigree-1240K.qQ20.L70.pileup"
for i in {0..10}; do 
  echo $i $(awk -vstep=$i 'BEGIN{n=0}{n+=$(4+(3*step))}END{print n/NR}' $BAM)
done
```
