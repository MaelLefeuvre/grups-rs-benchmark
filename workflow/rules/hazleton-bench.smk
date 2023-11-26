configfile: "./config/config.yml"
configfile: "./config/bam-urls.yml"

from os.path import basename


rule download_hazleton_bams:
    params:
        urls = config["hazleton-bams"]
    output:
        bams = expand("data/hazleton-bams/{bam}", bam=[basename(url) for url in config['hazleton-bams']])
    threads: 8
    shell: """
        for url in {params.urls}; do
            bam="$(basename $url)"
            echo "Downloading data/hazleton-bams/$bam" $url
            wget -qO "data/hazleton-bams/$bam" $url
        done
    """


rule create_hazleton_bamlist:
    input:
        bams    = rules.download_hazleton_bams.output.bams,
    output:
        bamlist = "results/hazleton/01-pileup/hazleton-bams.list"
    shell: """
        ls -1 {input.bams} > {output.bamlist}
    """


rule pileup_hazleton_bams:
    input:
        bams      = rules.download_hazleton_bams.output.bams,
        bamlist   = rules.create_hazleton_bamlist.output.bamlist,
        targets   = rules.download_reich_1240K.output.eigenstrat[0],
        reference = rules.download_reference_genome.output.reference
    output:
        pileup = "results/hazleton/01-pileup/hazleton.RBq{mq}Q{bq}.1240K.pileup"
    params:
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools mpileup -RB -q {wildcards.mq} -Q {wildcards.bq} \
        -l <(awk 'BEGIN{{OFS="\t"}}{{print $2, $4}}' {input.targets}) \
        -f {input.reference} \
        -b {input.bamlist} > {output.pileup}
    """


rule grups_rs_hazleton:
    input:
        pileup     = expand(rules.pileup_hazleton_bams.output.pileup, bq=25, mq=25),
        bamlist    = rules.create_hazleton_bamlist.output.bamlist,
        data       = expand(rules.GRUPS_generate_fst_set.output.fst, ped_pop = "EUR", cont_pop="AFR"),
        recomb_map = rules.download_hapmap.output.map,
        targets    = rules.download_reich_1240K.output.eigenstrat[0],
        pedigree   = "workflow/scripts/grups-rs/resources/pedigrees/extended_pedigree.txt",
    output:
        output_dir = directory("results/hazleton/02-run-grups/grups-rs-hazleton-north"),
        pwd        = "results/hazleton/02-run-grups/grups-rs-hazleton-north/hazleton.RBq25Q25.1240K.pwd",
        results    = "results/hazleton/02-run-grups/grups-rs-hazleton-north/hazleton.RBq25Q25.1240K.result",

    params:
        data_dir       = lambda w, input: dirname(input.data[0]),
        recomb_dir     = lambda w, input: dirname(input.recomb_map[0]),
        mode           = "fst-mmap",
        pedigree_pop   = config['hazleton']['pedigree-pop'],
        contam_pop     = config['hazleton']['contam-pop'],
        seq_error_rate = config['hazleton']['seq-error-rate'],
        reps           = config['hazleton']['reps'],
        seed           = config['hazleton']['seed']
    log:       "logs/grups-rs-hazleton/GRUPS_rs_hazleton.log"
    benchmark: "benchmarks/grups-rs-hazleton/GRUPS_rs_hazleton.tsv"
    conda:     "../envs/grups-rs.yml"
    shell: """
        grups pedigree-sims \
        --pileup {input.pileup} \
        --data-dir {params.data_dir} \
        --recomb-dir {params.recomb_dir} \
        --pedigree {input.pedigree} \
        --pedigree-pop {params.pedigree_pop} \
        --contam-pop {params.contam_pop} \
        --seq-error-rate {params.seq_error_rate} \
        --samples 0-$(($(cat {input.bamlist} | wc -l)-1)) \
        --sample-names $(grep -oP '[A-Z0-9]+[m|f](?=_hg19.bam$)' {input.bamlist} | tr '\n' ' ') \
        --reps {params.reps} \
        --mode {params.mode} \
        --output-dir {output.output_dir} \
        --seed {params.seed} \
        --self-comparison \
        --overwrite \
        --verbose > {log} 2>&1
    """
