configfile: "./config/config.yml"

from os.path import basename, splitext

rule download_1000_genomes:
    """
    Download 1000-genomes phase 3 SNP variant callset
    """
    params:
        url = expand("{url}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            url = config["ftp"]["1000g-url"],
            chr = "{chr}"
        )
    output:
        vcf = "data/1000g-phase3/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    log: "logs/00-netrules/download_1000_genomes/download_1000_genomes-chr{chr}.log"
    threads: 1
    shell: """
        wget -qO {output.vcf} {params.url} > {log} 2>&1
    """


rule fetch_samples_panel:
    """
    Download samples metadata from the 1000g FTP website
    """
    params:
        url   = f"{config['ftp']['1000g-url']}/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        panel = "data/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.ALL.panel"
    log: "logs/00-netrules/fetch_samples_panel.log"
    threads: 1
    shell: """
        wget -qO {output.panel} {params.url} > {log} 2>&1
    """


rule download_reich_1240K:
    """
    Download the 1240K dataset from Reich Lab's website.
    """
    params:
        tarball    = config['ftp']['1240k-url'],
        output_dir = lambda wildcards, output: dirname(output.eigenstrat[0])
    output:
        eigenstrat = multiext("data/aadr-1240k-dataset/v52.2_1240K_public", ".snp", ".ind", ".geno")
    log: "logs/00-netrules/download_reich_1240K.log"
    threads: 1
    shell:"""
        wget -qO- {params.tarball} | tar -xvf- -C {params.output_dir} > {log} 2>&1
    """


rule download_hapmap:
    """
    Download the 2011 HapMapII recombination map from ncbi.
    """
    params:
        tarball    = config['ftp']['hapmap-url'],
        output_dir = lambda wildcards, output: dirname(output.map[0])
    output:
        map     = expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=range(1, 23)),
        exclude = temp(expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=["X", "X_par1", "X_par2"])),
        readme  = temp("data/recombination-maps/HapMapII_GRCh37/README.txt")
    log: "logs/00-netrules/download_hapmap.log"
    threads: 1
    shell: """
        wget -qO- {params.tarball} | tar -xvzf- -C {params.output_dir} > {log} 2>&1
    """

rule download_reference_genome:
    output:
        reference = "data/reference/{file}".format(file = splitext(basename(config['ftp']['reference']))[0]),
        fai       = "data/reference/{file}.fai".format(file=splitext(basename(config['ftp']['reference']))[0])
    params:
        url = config['ftp']['reference'],
        fai = splitext(config['ftp']['reference'])[0] + ".fai"
    shell: """
        wget -qO {output.fai} {params.fai}
        wget -qO {output.reference}.gz {params.url}
        truncate -s 891946027 {output.reference}.gz
        gunzip {output.reference}.gz
    """

