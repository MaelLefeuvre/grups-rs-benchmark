configfile: "./config/config.yml"

from os.path import dirname
import sys

def get_bench_replicates():
    1 if "grups_rs_perf_bench" in sys.argv else config['bench']['bench-replicates']


rule subset_1000g_samples_panel:
    input:
        panel = rules.fetch_samples_panel.output.panel
    output:
        panel = expand("{directory}/integrated_call_samples_v3.20130502.{source}-{cont}.panel",
            directory = dirname(rules.fetch_samples_panel.output.panel),
            source = "{source}",
            cont = "{cont}"
        )
    threads: 1
    shell: """
        grep -P '{wildcards.source}|{wildcards.cont}' {input.panel} > {output.panel}
    """


rule filter_1000g:
    input:
        vcf   = rules.download_1000_genomes.output.vcf,
        panel = expand(rules.subset_1000g_samples_panel.output.panel,
            source = config['grups-rs']['pedigree-pop'],
            cont   = config['grups-rs']['contam-pop']
        )
    output:
        vcf = expand("data/1000g-phase3/01-filtered/{source}-{cont}.chr{{chr}}.phase3_shapeit2_mvncall_integrated_v5b.20130502.snps.vcf.gz",
            source = config['grups-rs']['pedigree-pop'],
            cont   = config['grups-rs']['contam-pop']
        )
    conda: "../envs/bcftools-1.15.yml"
    threads: 1
    shell: """
        bcftools view --threads {threads} -vsnps -m2 -M2 -S <(awk '{{print $1}}' {input.panel}) -Oz -o {output.vcf} {input.vcf}
    """


rule tabix:
    input:
        vcf = "{directory}/{file}.{ext}"
    output:
        tbi = "{directory}/{file}.{ext}.tbi"
    wildcard_constraints:
        ext="vcf|vcf.gz"
    conda: "../envs/bcftools-1.15.yml"
    threads: 1
    shell: """
        tabix {input.vcf}
    """


rule GRUPS_generate_fst_set:
    """
    Generate a .fst and .fst.frq dataset from a set of VCF files.
    """
    input:
        data    = expand(rules.download_1000_genomes.output.vcf, chr=range(1,23)),
        panel   = rules.fetch_samples_panel.output.panel
    output:
        fst     = protected(expand(
            "data/fst/g1k-phase3-v5/{ped_pop}-{cont_pop}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-{ped_pop}-{cont_pop}.{ext}",
            ped_pop  = "{ped_pop}",
            cont_pop = "{cont_pop}",
            chrom    = range(1,23),
            ext      = ["fst", "fst.frq"]
        ))
    log:       "logs/grups-rs/grups_generate_fst_set/{ped_pop}-{cont_pop}-GRUPS_generate_fst_set.log"
    benchmark: "benchmarks/grups-rs/grups_generate_fst_set/{ped_pop}-{cont_pop}-GRUPS_generate_fst_set.tsv"
    conda:     "../envs/grups-rs.yml"
    threads:   22
    shell: """
        grups fst \
        --vcf-dir $(dirname {input.data} | uniq) \
        --output-dir $(dirname {output.fst} | uniq) \
        --pop-subset {wildcards.ped_pop} {wildcards.cont_pop} \
        --panel {input.panel} \
        --threads {threads} \
        --verbose > {log} 2>&1
    """


rule run_GRUPS:
    input:
        pileup     = "results/input-pileups/{overlap}/{depth}X-{description}-{overlap}.{id}.pileup",
        data       = expand(rules.filter_1000g.output.vcf, chr=range(1,23)) if config['grups-rs']['mode'] == "vcf" else expand(
            rules.GRUPS_generate_fst_set.output.fst,
            ped_pop  = config['grups-rs']['pedigree-pop'],
            cont_pop = config['grups-rs']['contam-pop']
        ),
        panel      = expand(rules.subset_1000g_samples_panel.output.panel,
            source = config['grups-rs']['pedigree-pop'],
            cont   = config['grups-rs']['contam-pop']
        ),
        recomb_map = rules.download_hapmap.output.map,
        targets    = rules.download_reich_1240K.output.eigenstrat[0],
        pedigree   = config['grups-rs']['pedigree']
    output:
        output_dir   = directory("results/{depth}X/{overlap}/{nsamples}/{description}-{id}/"),
        results      = multiext("results/{depth}X/{overlap}/{nsamples}/{description}-{id}/{depth}X-{description}-{overlap}.{id}", ".pwd", ".result")
    params:
        sample_names = config['grups-rs']['sample-names'],
        data_dir     = lambda wildcards, input: dirname(input.data[0]),
        recomb_dir   = lambda wildcards, input: dirname(input.recomb_map[0]),
        pedigree_pop = config['grups-rs']['pedigree-pop'],
        contam_pop   = config['grups-rs']['contam-pop'],
        reps         = config['grups-rs']['reps'],
        mode         = config['grups-rs']['mode'],
        min_depth    = config['grups-rs']['min-depth'],
        min_quality  = config['grups-rs']['min-qual'],
        maf          = config['grups-rs']['maf'],
        seed         = config['grups-rs']['seed']
    resources:
        # Crude estimation, based on a previous benchmark
        mem_mb       = lambda w: (int(w.overlap) * 0.001898) + 1900
    log: "logs/grups-rs-bench/{depth}X/{overlap}/{description}-{nsamples}.{id}.log"
    benchmark: repeat("benchmarks/grups-rs-bench/{depth}X/{overlap}/{description}-{nsamples}.{id}.log", get_bench_replicates())
    conda: "../envs/grups-rs.yml"
    shell: """
        grups-rs pedigree-sims \
        --pileup {input.pileup} \
        --data-dir {params.data_dir} \
        --panel {input.panel} \
        --recomb-dir {params.recomb_dir} \
        --pedigree {input.pedigree} \
        --pedigree-pop {params.pedigree_pop} \
        --contam-pop {params.contam_pop} \
        --min-depth {params.min_depth} \
        --samples 0-{wildcards.nsamples} \
        --sample-names {params.sample_names} \
        --reps {params.reps} \
        --mode {params.mode} \
        --output-dir {output.output_dir} \
        --maf {params.maf} \
        --min-qual {params.min_quality} \
        --seed {params.seed} \
        --overwrite \
        --quiet > {log} 2>&1
    """
