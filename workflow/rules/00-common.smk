configfile: "./config/config.yml"

import re

def get_base_overlap():
    return int(re.findall(r'(?<=[-])[0-9]+(?=[.])', config['bench']['input-pileup'])[0])

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
