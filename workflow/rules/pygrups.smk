
rule get_observed_pwd:
    input:
        pileup = "results/input-pileups/{coverage}X-EUR-pedigree-1240K.qQ20.L70-{SNPs}.pileup", 
    output:
        pwd    = "results/pygrups/{coverage}X/{SNPs}/pwd-distributions/pwd-from-stdin.out"
    params:
        reps = config['kinship']['pygrups']['pwd-from-stdin']['reps'],
        min_depth = config['kinship']['pygrups']['min-depth'],
        min_qual = config['kinship']['pygrups']['min-qual'],
        ignore_dels = 1 if config['kinship']['pygrups']['pwd-from-stdin']['ignore-dels'] == True else 0
    conda: "../envs/pygrups.yml"
    threads: 1
    shell: """
        for ((i=0; i<{params.reps}; i++)); do 
            cat {input.pileup} | PWD_from_stdin.py --chr 1-22 --min_depth {params.min_depth},{params.min_depth] --min_qual {params.min_qual} --ignore_dels {params.ignore_dels} --quiet 1 >> {output.pwd}
        done
    """

rule compute_summary_stats:
    input:
        pwd = rules.get_observed_pwd.output.pwd
    output:
        stats = "results/pygrups/{coverage}X/{SNPs}/pwd-distributions/pwd-from-stdin.stats"
    conda: "../envs/pygrups.yml"
    threads: 1
    shell: """
        cat {input.pwd} | awk 'NF>0' | awk 'BEGIN{{
            max=0; min=10000;
        }}{{
            sum+=$1;
            sumsq+=$1*$1;
            if ($1>max){max=$1};
            if ($1<min){min=$1};
        }}END{{
            print "n mean std max min";
            print NR " " sum/NR " " sqrt(sumsq/NR - 'sum/NR)^2) " " max " " min;
        }}' > {output.stats}
    """

rule filter_sites:
    input:
        pileup = rules.get_observed_pwd.input.pileup
    output:
        pileup = "results/pygrups/{coverage}X/{SNPs}/filtered-pileups/pwd-from-stdin-simspecific.pileup.gz"
    params:
        min_depth = config['kinship']['pygrups']['min-depth'],
        min_qual = config['kinship']['pygrups']['min-qual'],
        ignore_dels = 1 if config['kinship']['pygrups']['pwd-from-stdin']['ignore-dels'] == True else 0
    conda: "../envs/pygrups.yml"
    threads: 1
    shell: """
        cat {input.pileup} | PWD_from_stdin.py --chr 1-22 --min_depth {params.min_depth},{params.min_depth} --min_qual {params.min_qual} --ignore_dels {params.ignore_dels} --quiet 1 --filter_sites 1 | gzip > {output.pileup}
    """

rule split_sites:
    input:
        pileup  = rules.filter_sites.output.pileup
    output:
        pileups = expand("results/pygrups/{{coverage}}X/{{SNPs}}/filtered-pileups/pwd-from-stdin-simspecific.pileup.chr{chr}.gz", chr=range(1,23))
    params:
        out_dir = "results/pygrups/{coverage}X/{SNPs}/filtered-pileups/"
    threads: 1
    shell: """
        mkdir -p {params.out_dir};
        for ((i=1; i<23; i++)); do 
            zcat {input.pileup} | awk -v y="$i" '($1==y)' | gzip > {params.out_dir}/pwd-from-stdin-simspecific.pileup.chr${{i}}.gz ;
        done
    """

rule symlink_1000g_dataset:
    input:
        vcf = expand(rules.filter_1000g.output.vcf, chr=range(1,23)),
        tbi = expand([vcf + ".tbi" for vcf in rules.filter_1000g.output.vcf], chr=range(1, 23))
    output:
        symlinks = expand("data/1000g-phase3/02-symlinked/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz",
          chr    = range(1, 23)
        ),
        symdir   = directory("data/1000g-phase3/02-symlinked")
    params:
        poptag   = config['kinship']['pygrups']['pedigree-sims']['pedigree-pop'] + "-" + config['kinship']['pygrups']['pedigree-sims']['contam-pop']
    threads: 1
    shell: """
        echo {params.poptag}; \
        mkdir -p {output.symdir}; \
        ln -srt {output.symdir} {input.vcf}; \
        ln -srt {output.symdir} {input.tbi}; \
        for file in {output.symdir}/*; do
            new=$(echo $file | sed 's/_v5b/_v5/' | sed 's/{params.poptag}/ALL/' | sed 's/snps/genotypes/');
            mv $file $new;
        done
    """

rule pygrups_pedigree_sims:
    input:
        pileup     = rules.filter_sites.output.pileup,
        split      = rules.split_sites.output.pileups,
        vcfs       = rules.symlink_1000g_dataset.output.symlinks,
        recomb_map = rules.download_hapmap.output.map,
        data_dir   = dirname(rules.symlink_1000g_dataset.output.symlinks[0]),
        recomb_dir = dirname(rules.download_hapmap.output.map[0]),
    output:
        done = touch("results/pygrups/{coverage}X/{SNPs}/pedigree-sims/pedigree-sims.done")
    params:
        out_dir    = "results/pygrups/{coverage}X/{SNPs}/pedigree-sims/",
        ped_pop    = config['kinship']['pygrups']['pedigree-sims']['pedigree-pop'],
        contam_pop = config['kinship']['pygrups']['pedigree-sims']['contam-pop'],
        ds_rate    = config['kinship']['pygrups']['pedigree-sims']['downsample-rate'],
        c_rate     = config['kinship']['pygrups']['pedigree-sims']['contam-rate'],
        q_rate     = config['kinship']['pygrups']['pedigree-sims']['seq-error-rate'],
        reps       = config['kinship']['pygrups']['pedigree-sims']['reps'],
        label      = config['kinship']['pygrups']['pedigree-sims']['labels'],
        min_qual   = config['kinship']['pygrups']['min-qual'],
    log: "logs/pygrups/pygrups_pedigree_sims-{coverage}X-{SNPs}.log"
    benchmark: "benchmarks/pygrups/pygrups_pedigree_sims-{coverage}X-{SNPs}.tsv"
    conda: "../envs/pygrups.yml"
    threads: 1
    shell: """
        pedigree_sims.py \
        --label {params.label} \
        --ds_rate {params.ds_rate} \
        --c_rate {params.c_rate},{params.c_rate} \
        --q_rate {params.q_rate},{params.q_rate} \
        --reps {params.reps} \
        --chr 1-22 \
        --include all \
        --min_qual {params.min_qual} \
        --data_dir {input.data_dir}/ \
        --out {params.out_dir} \
        --recomb_dir {input.recomb_dir}/ \
        --pedigree_pop {params.ped_pop} \
        --contam_pop {params.contam_pop},{params.contam_pop} \
        > {log} 2>&1
    """
