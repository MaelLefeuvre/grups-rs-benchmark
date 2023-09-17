configfile: "./config/config.yml"

from os.path import dirname

rule subset_1000g_samples_panel:
    input:
        panel = rules.fetch_samples_panel.output.panel
    output:
        panel = expand("{directory}/integrated_call_samples_v3.20130502.{source}-{cont}.panel",
            directory = dirname(rules.fetch_samples_panel.output.panel),
            source    = config['kinship']['GRUPS']['pedigree-pop'],
            cont      = config['kinship']['GRUPS']['contam-pop'],
        )
    params:
        source = config['kinship']['GRUPS']['pedigree-pop'],
        cont   = config['kinship']['GRUPS']['contam-pop'],
    threads: 1
    shell: """
        grep -P '{params.source}|{params.cont}' {input.panel} > {output.panel}
    """


rule filter_1000g:
    input:
        vcf   = rules.download_1000_genomes.output.vcf,
        panel = rules.subset_1000g_samples_panel.output.panel,
    output:
        vcf = expand("data/1000g-phase3/01-filtered/{source}-{cont}.chr{{chr}}.phase3_shapeit2_mvncall_integrated_v5b.20130502.snps.vcf.gz",
            source = config['kinship']['GRUPS']['pedigree-pop'],
            cont   = config['kinship']['GRUPS']['contam-pop']
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
            ped_pop  = config["kinship"]["GRUPS"]["pedigree-pop"],
            cont_pop = config["kinship"]["GRUPS"]["contam-pop"],
            chrom    = range(1,23),
            ext      = ["fst", "fst.frq"]
        ))
    params:
        pedigree_pop = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop   = config["kinship"]["GRUPS"]["contam-pop"]
    log:       "logs/grups-rs/grups_generate_fst_set/{params.pedigree_pop}-{params.contam_pop}-GRUPS_generate_fst_set.log"
    benchmark: "benchmarks/grups-rs/grups_generate_fst_set/{params.pedigree_pop}-{params.contam_pop}-GRUPS_generate_fst_set.tsv"
    conda:     "../envs/grups-rs.yml"
    threads:   22
    shell: """
        grups fst \
        --vcf-dir $(dirname {input.data} | uniq) \
        --output-dir $(dirname {output.fst} | uniq) \
        --pop-subset {params.pedigree_pop} {params.contam_pop} \
        --panel {input.panel} \
        --threads {threads} \
        --verbose > {log} 2>&1
    """


rule run_GRUPS:
    input:
        pileup     = "results/input-pileups/{coverage}X-EUR-pedigree-1240K.qQ20.L70-{SNPs}.pileup",
        data       = expand(rules.filter_1000g.output.vcf, chr=range(1,23)) if config["kinship"]["GRUPS"]["mode"] == "vcf" else rules.GRUPS_generate_fst_set.output.fst,
        panel      = rules.subset_1000g_samples_panel.output.panel, 
        recomb_map = rules.download_hapmap.output.map,
        targets    = rules.download_reich_1240K.output.eigenstrat[0],
        pedigree   = config["kinship"]["GRUPS"]["pedigree"]
    output:
        output_dir   = directory("results/{coverage}X/{SNPs}/{rep}/"),
        results      = multiext("results/{coverage}X/{SNPs}/{rep}/{coverage}X-EUR-pedigree-1240K.qQ20.L70-{SNPs}", ".pwd", ".result")
    params:
        sample_names = config['kinship']['GRUPS']['sample-names'],
        data_dir     = lambda wildcards, input: dirname(input.data[0]),
        recomb_dir   = lambda wildcards, input: dirname(input.recomb_map[0]),
        pedigree_pop = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop   = config["kinship"]["GRUPS"]["contam-pop"],
        reps         = config["kinship"]["GRUPS"]["reps"],
        mode         = config["kinship"]["GRUPS"]["mode"],
        min_depth    = config["kinship"]["GRUPS"]["min-depth"],
        min_quality  = config["kinship"]["GRUPS"]["min-qual"],     
        maf          = config["kinship"]["GRUPS"]["maf"],
        grups        = config["kinship"]["GRUPS"]["path"]
    resources:
        # Crude estimation, based on a previous benchmark
        mem_mb       = lambda w: (int(w.SNPs) * 0.001898) + 1900
    log: "logs/run_GRUPS-{coverage}X-{SNPs}-{rep}.log"
    benchmark: repeat("benchmarks/run_GRUPS-{coverage}X-{SNPs}-{rep}.log", 10)
    conda: "../envs/grups-rs.yml"
    shell: """
        grups pedigree-sims \
        --pileup {input.pileup} \
        --data-dir {params.data_dir} \
        --panel {input.panel} \
        --recomb-dir {params.recomb_dir} \
        --pedigree {input.pedigree} \
        --pedigree-pop {params.pedigree_pop} \
        --contam-pop {params.contam_pop} \
        --min-depth {params.min_depth} \
        --samples 0-{wildcards.rep} \
        --sample-names {params.sample_names} \
        --reps {params.reps} \
        --mode {params.mode} \
        --output-dir {output.output_dir} \
        --maf {params.maf} \
        --min-qual {params.min_quality} \
        --print-blocks \
        --overwrite \
        --quiet > {log} 2>&1
    """


