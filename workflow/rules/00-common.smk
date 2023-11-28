configfile: "./config/config.yml"

from os.path import basename, splitext
import re, sys

def get_bench_replicates():
    return 1 if "accuracy" in sys.argv else config['bench']['bench-replicates']

def koszyce_output_results(wildcards):
    hypotheses = [basename(splitext(pedigree)[0]) for pedigree in config['koszyce']['hypotheses']]
    return expand(rules.grups_rs_koszyce.output, hypothesis = hypotheses)


def get_expected_rels():
    with open(config['bench']['expected']) as f:
        _ = f.readline()
        return set([rel for rel in map(lambda x: x.strip('\n').split('\t')[3], f.readlines())])

def get_base_overlap():
    return int(re.findall(r'(?<=[-])[0-9]+(?=[.]pileup.gz)', config['bench']['input-pileup'])[0])

def get_description():
    return re.findall(r'(?<=[/][0-9]X-).*(?=[-][0-9]+[.]pileup.gz$)', config['bench']['input-pileup'])[0]

def get_depth():
    return re.findall(r'(?<=[/])[0-9](?=X)', config['bench']['input-pileup'])[0]

def compute_overlap_powerlist():
    return sum([[int(get_base_overlap()/(2**i))] * 2**i for i in range(1,config['bench']['n-subsamples'])], [])

def compute_subsample_ids_powerlist():
    return sum([[n for n in range(2**i)] for i in range(1,config['bench']['n-subsamples'])], [])

def get_overlap_list():
    return [int(get_base_overlap()/(2**i)) for i in range(config['bench']['n-subsamples'])]

def get_subsamples_list():
    return [2**i for i in range(config['bench']['n-subsamples'])]
