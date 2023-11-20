configfile: "./config/config.yml"

import gzip, re
from os.path import basename, splitext


rule untar_original_pileup:
    input:  config['bench']['input-pileup']
    output: expand("results/input-pileups/{overlap}/{depth}X-{description}-{overlap}.0.pileup", overlap=get_base_overlap(), description=get_description(), depth=get_depth())
    shell:  "zcat {input} > {output}"


rule subdivide_pileups:
    input: 
        original = rules.untar_original_pileup.output
    output: 
        subsamples = expand(
            expand("results/input-pileups/{{overlap}}/{depth}X-{description}-{{overlap}}.{{subsample_id}}.pileup",
                depth    = get_depth(),
                description = get_description()
            ),
            zip,
            overlap      = compute_overlap_powerlist(),
            subsample_id = compute_subsample_ids_powerlist(),
        )
    params:
        depth = config['bench']['n-subsamples']
    shell: "workflow/scripts/subdivide-pileup.py {input.original} {params.depth}"
