# GRUPS-rs benchmark
[![DOI](https://zenodo.org/badge/695099930.svg)](https://zenodo.org/doi/10.5281/zenodo.10389549)

A set of snakemake workflows and jupyter notebooks to benchmark the runtime, memory consumption, and performance of GRUPS-rs

## Installation

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

## Requirements

- This pipeline is self-contained and will therefore download any data required to run to completion. about 85GB of disk is required to run the whole workflow.
- If you wish to apply GRUPS-rs on the Hazleton-North dataset with this pipeline, a computer with at least 32GB of RAM is recommended. 

## Recommendations

Note that as with many other softwares, the runtime performance of GRUPS-rs can be heavily influenced by the type of disk in which its input data is stored. If possible, make sure to run this workflow from within an SSD drive:, Having as much free space as possible within this drive (e.g.: at least 25%), can also be critical, as the read/write speed of SSDs is notoriously known to plummet when close to being full. If runtime performance is critical to you, running the following UNIX commands before starting this workflow might also be of help:

```bash
sudo fstrim -av
sudo cpupower frequency-set --governor performance
```

---

## Usage 
### Running GRUPS-rs' runtime benchmark

Optionally tweak the `config/config.yml` configuration file beforehand. 

Then, simply run the snakemake pipeline using the following command

```bash
conda activate snakemake-7.20.0
snakemake runtime --cores 1 --use-conda --conda-frontend mamba --printshellcmds
```

This should:
1. Download and preprocess the 1000g-phase3 dataset in `data/1000g-phase3`
2. Download the HapMapII recombination map in `data/recombination-maps`
3. Download the AADR v52.2 1240K variant callser in `data/aadr-1240-dataset`
4. FSA-encode the 1000g-phase3 dataset (if `mode` param is set to `fst` or `fst-mmap`)
5. Apply GRUPS-rs 10 times on the following parameter space:
  - number of pairwise comparisons: (1, 3, 6, 10, 15, 21, 28, 36, 45, 55]
  - Average number of positions within the input pileup: [2246, 4493, 8987, 17974, 35948, 71897, 143794, 287588, 575177, 1150354]
6. Plot a 3D surface plot of the average runtime and maximum resident set size

###### Output:
- raw benchmark results should be located in `benchmarks`
- 3D-surface plots should be located in `results/bench-plots`

## Testing the accuracy of GRUPS-rs on the benchmark input dataset

A rudimentary benchmark of the classification accuracy of GRUPS-rs can be obtained by running the following workflow:

```bash
conda activate snakemake-7.20.0
snakemake accuracy --cores `nproc` --use-conda --conda-frontend mamba --printshellcmds
```

This should:
1. Perform cascading subsamples from an input pileup file containing 1150354 positions, thus generating independant input files in powers of two: 
   | avg. number of positions within the input pileup | number of independent files |
   | ------------------------------------------------ | --------------------------- |
   | 1150354                                          |   1                         |
   | 575177                                           |   2                         |
   | 287588                                           |   4                         |
   | 143794                                           |   8                         |
   | 71897                                            |  16                         |
   | 35948                                            |  32                         |
   | 17974                                            |  64                         |
   | 8987                                             | 128                         |
   | 4493                                             | 256                         |
   | 2246                                             | 512                         |

2. Apply GRUPS-rs once on each of these files.
3. Generate per-class classification performance plots for every expected relationship found within the pileup file. The 11 individuals within this pileup file were generated using [ped-sim](https://github.com/williamslab/ped-sim) [(Caballero et al. 2019)](https://doi.org/10.1371/journal.pgen.1007979) and harbors the following pairwise relationships:
   | relationship | n   |
   | ------------ | --- |
   | Unrelated    | 31  |
   | First        | 12  |
   | Second       | 08  |
   | Third        | 03  |
   | Self         | 01  |


Hence the total number of predictions applied by GRUPS-rs during this benchmark is as follows:
| avg. number of positions within the input pileup | number of independent files | total | n-Unrelated | n-First | n-Second | n-Third | n-Self |
| ------------------------------------------------ | --------------------------- | ----- | ----------- | ------- | -------- | ------- | ------ |
| 1150354                                          |   1                         | 55    | 31          | 12      | 8        | 3       | 1      |
| 575177                                           |   2                         | 110   | 62          | 24      | 16       | 6       | 2      |
| 287588                                           |   4                         | 220   | 124         | 48      | 32       | 12      | 4      |
| 143794                                           |   8                         | 440   | 248         | 96      | 64       | 24      | 8      |
| 71897                                            |  16                         | 880   | 496         | 192     | 128      | 48      | 16     |
| 35948                                            |  32                         | 1760  | 992         | 384     | 256      | 96      | 32     |
| 17974                                            |  64                         | 3520  | 1984        | 768     | 512      | 192     | 64     |
| 8987                                             | 128                         | 7040  | 3968        | 1536    | 1024     | 384     | 128    |
| 4493                                             | 256                         | 14080 | 7936        | 3072    | 2048     | 768     | 256    |
| 2246                                             | 512                         | 28160 | 15872       | 6144    | 4096     | 1536    | 512    |

##### Output
- Per-class classification performance plots should be located in `results/accuracy-plots`


---

## Running PyGrups' runtime benchmark

This workflow will fire up a benchmark of the previous version of GRUPS, written in python (See [grups (GitHub)](https://github.com/sameoldmike/grups), and [(Martin et al 2017](https://doi.org/10.1111/mec.14188))

```bash
conda activate snakemake-7.20.0
snakemake pygrups --cores 1 --use-conda --conda-frontend mamba --printshellcmds
```

###### Output:
- Raw benchmark results for PyGrups should be located in `benchmarks/pygrups`

---

## Applying GRUPS-rs on the Hazleton-North dataset
   
**Note:** Due to the presence of 630 pairwise comparaisons and relatively high sequencing coverages for some individuals, Applying GRUPS-rs on the Hazleton-North dataset will require about ~32 GB of virtual memory. 

1. Run the snakemake pipeline using the following command
   ```bash
   conda activate snakemake-7.20.0
   snakemake hazleton --cores `nproc` --use-conda --conda-frontend mamba --printshellcmds
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

---

## Applying GRUPS-rs on the Koszyce dataset

1. Run the snakemake pipeline using the following command
   ```bash
   conda activate snakemake-7.20.0
   snakemake koszyce --cores `nproc` --use-conda --conda-frontend mamba --printshellcmds
   ```

###### Output:
- Raw results and output tables should be located in `results/koszyce`
- Interactive plots can be visualized using [grups.plots](https://github.com/MaelLefeuvre/grups.plots). e.g.:
    ```r
    grups.plots::app("results/koszyce/02-run-grups/koszyce-first-degree/")
    ```
