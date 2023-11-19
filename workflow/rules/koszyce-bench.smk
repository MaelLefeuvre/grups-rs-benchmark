configfile: "./config/config.yml"
configfile: "./config/bam-urls.yml"

from os.path import basename
import re

def search_koszyce_url(wildcards):
    for url in config['koszyce-bams']:
        if wildcards.sample in url:
            return url
    raise RuntimeError(f"Cannot find wildcard {wildcards.sample} in urls")

rule download_koszyce_bams:
    params:
        url = search_koszyce_url
    output:
        bam = "data/koszyce-bams/{sample}.bam"
    threads: 1
    shell: """
        echo Downloading {output.bam} using {params.url}
        wget -qO {output.bam} {params.url}
    """

rule samtools_index:
    input: "{bam}"
    output: "{bam}.bai"
    threads: 16
    conda: "../envs/samtools-1.15.yml"
    shell: "samtools index -@ {threads} {input}"


rule mapdamage_koszyce:
    input:
        bam       = rules.download_koszyce_bams.output.bam,
        bai       = "data/koszyce-bams/{sample}.bam.bai",
        reference = rules.download_reference_genome.output.reference
    output:
        rescale_dir = directory("results/koszyce/01-mapdamage/{sample}.rescaled"),
        c_to_t_freq = "results/koszyce/01-mapdamage/{sample}.rescaled/5pCtoT_freq.txt",
        g_to_a_freq = "results/koszyce/01-mapdamage/{sample}.rescaled/3pGtoA_freq.txt"
    log:       "logs/grups-rs-koszyce/mapdamage_koszyce/{sample}.log"
    benchmark: "benchmarks/grups-rs-koszyce/mapdamage_koszyce/{sample}.tsv"
    conda:     "../envs/mapdamage-2.2.1.yml"
    shell: """
        mapDamage -i {input.bam} -r {input.reference} --folder {output.rescale_dir} --verbose > {log} 2>&1
    """


def match_position(file, threshold=0.01):
    with open(file, 'r') as freqs:
        for line in freqs:
            pos, freq = line.split()
            try:
                if float(freq) < threshold:
                    return pos
            except ValueError:
                continue # header
        return pos

rule trimbam_koszyce:
    input:
        bam         = rules.download_koszyce_bams.output.bam,
        c_to_t_freq = rules.mapdamage_koszyce.output.c_to_t_freq,
        g_to_a_freq = rules.mapdamage_koszyce.output.g_to_a_freq
    output:
         bam  = "results/koszyce/02-trimmed/{sample}.trimmed.bam"
    params:
        left      = lambda w, input: match_position(input.c_to_t_freq, threshold =  0.01),
        right     = lambda w, input: match_position(input.g_to_a_freq, threshold =  0.01),
    log:       "logs/grups-rs-koszyce/trimbam_koszyce/{sample}.log"
    benchmark: "benchmarks/grups-rs-koszyce/trimbam_koszyce/{sample}.tsv"
    conda:     "../envs/bamutil-1.0.15.yml"
    shell: """
        bam trimBam {input.bam} {output.bam} --left {params.left} --right {params.right} > {log} 2>&1
    """

rule create_koszyce_bamlist:
    input:
        bams = expand(rules.trimbam_koszyce.output.bam, sample = [basename(url).removesuffix('.bam') for url in config['koszyce-bams']])
    output:
        bamlist = "results/koszyce/01-pileup/koszyce-bams.list"
    shell: """
        ls -1 {input.bams} | LC_ALL=C sort -n -k2 -t_ > {output.bamlist}
    """

rule pileup_koszyce_bams:
    input:
        bams      = rules.create_koszyce_bamlist.input.bams,
        bamlist   = rules.create_koszyce_bamlist.output.bamlist,
        targets   = rules.download_reich_1240K.output.eigenstrat[0],
        reference = rules.download_reference_genome.output.reference
    output:
        pileup = "results/koszyce/01-pileup/koszyce.RBq{mq}Q{bq}.1240K.pileup"
    params:
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools mpileup -RB -q {wildcards.mq} -Q {wildcards.bq} \
        -l <(awk 'BEGIN{{OFS="\t"}}{{print $2, $4}}' {input.targets}) \
        -f {input.reference} \
        -b {input.bamlist} > {output.pileup}
    """

def parse_sample_names(wildcards, input):
    regex = re.compile(r'Koszyce-3_[0-9]+')
    remove = 'oszyce-3_'
    with open(input.bamlist, "r") as bams:
        return [regex.search(basename(bam.strip())).group().replace(remove, '') for bam in bams]

rule grups_rs_koszyce:
    input:
        pileup     = expand(rules.pileup_koszyce_bams.output.pileup, bq=25, mq=25),
        bamlist    = rules.create_koszyce_bamlist.output.bamlist,
        data       = expand(rules.GRUPS_generate_fst_set.output.fst, ped_pop = "EUR", cont_pop="EUR"),
        recomb_map = rules.download_hapmap.output.map,
        targets    = rules.download_reich_1240K.output.eigenstrat[0],
        pedigree   = "resources/koszyce-pedigree.txt",
    output:
        output_dir = directory("results/koszyce/02-run-grups/grups-rs-koszyce"),
        pwd        = "results/koszyce/02-run-grups/grups-rs-koszyce/koszyce.RBq25Q25.1240K.pwd",
        results    = "results/koszyce/02-run-grups/grups-rs-koszyce/koszyce.RBq25Q25.1240K.result",

    params:
        data_dir       = lambda w, input: dirname(input.data[0]),
        recomb_dir     = lambda w, input: dirname(input.recomb_map[0]),
        mode           = "fst-mmap",
        pedigree_pop   = "EUR",
        contam_pop     = "EUR",
        seq_error_rate = 0.0,
        reps           = 500,
        seed           = 42,
        sample_names   = lambda w, input : parse_sample_names(w, input)
    log:       "logs/grups-rs-koszyce/GRUPS_rs_koszyce.log"
    benchmark: "benchmarks/grups-rs-koszyce/GRUPS_rs_koszyce.tsv"
    conda:     "../envs/grups-rs.yml"
    shell: """
    /data/mlefeuvre/dev/grups-rs/target/release/grups-rs pedigree-sims \
        --pileup {input.pileup} \
        --data-dir {params.data_dir} \
        --recomb-dir {params.recomb_dir} \
        --pedigree {input.pedigree} \
        --pedigree-pop {params.pedigree_pop} \
        --contam-pop {params.contam_pop} \
        --seq-error-rate {params.seq_error_rate} \
        --samples 0-$(($(cat {input.bamlist} | wc -l)-1)) \
        --sample-names {params.sample_names} \
        --reps {params.reps} \
        --mode fst-mmap \
        --output-dir {output.output_dir} \
        --seed {params.seed} \
        --self-comparison \
        --overwrite \
        --verbose > {log} 2>&1
    """