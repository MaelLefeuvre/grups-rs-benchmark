# GRUPS-rs benchmark

A set of snakemake workflows and jupyter notebooks to benchmark the runtime, memory consumption, and performance of GRUPS-rs

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
   conda env create -f ./envs/snakemake-7.20.0.yml
   ```

## Requirements:

- This pipeline is self-contained and will therefore download any data required to run to completion. about 85GB of disk is required to run the whole workflow.
- If you wish to apply GRUPS-rs on the Hazleton-North dataset with this pipeline, a computer with at least 32GB of RAM is recommended. 


## Running GRUPS-rs' runtime benchmark

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

###### Output:
- raw benchmark results should be located in `benchmarks`
- 3D-surface plots should be located in `results/bench-plots`

### Compute the average local sequencing depth of each input pileup file
```Shell
BAM="original-data/1X-EUR-pedigree-1240K.qQ20.L70.pileup"
for i in {0..10}; do 
  echo $i $(awk -vstep=$i 'BEGIN{n=0}{n+=$(4+(3*step))}END{print n/NR}' $BAM)
done
```

## Running PyGrups' runtime benchmark
```bash
conda activate snakemake-7.20.0
snakemake pygrups_bench --cores 1 --use-conda --conda-frontend mamba --printshellcmds
```

###### Output:
- Raw benchmark results for PyGrups should be located in `benchmarks/pygrups`

## Applying GRUPS-rs on the Hazleton-North dataset
   
**Note:** Due to the presence of 630 pairwise comparaisons and relatively high sequencing coverages for some individuals, Applying GRUPS-rs on the Hazleton-North dataset will require about ~32 GB of virtual memory. 

1. Run the snakemake pipeline using the following command
   ```bash
   conda activate snakemake-7.20.0
   snakemake grups_rs_hazleton --cores `nproc` --use -conda --conda-frontend mamba --printshellcmds
   ```

   This should:
   1. Download Fowler et al's supplementary data file within `data/'
   2. Download All the Hazelton North samples (n=35) from the European Nucleotide Archive [PRJEB46958](https://www.ebi.ac.uk/ena/browser/view/PRJEB46958), within the `data/hazleton-bams` subdirectory.
   3. Pileup these samples using samtools, while targeting the AADR v52.2 '1240K' snp callset.
   4. Apply GRUPS-rs on this pileup (`pedigree-pop: "EUR"`, `seq-error-rate: 0.0`, `reps: 1000`)

3. Create a conda environment for jupyter-notebook:
   ```bash
   conda env create -f envs/hazleton-analysis.yml
   bash envs/hazleton-analysis.post-deploy.sh
   ```
4. Open and run the jupyter-notebook environment
   ```bash
   conda activate hazleton-analysis
   jupyter-notebook ./notebooks
   ```
###### Output:
- Results and plots should be located in `results/hazleton`

