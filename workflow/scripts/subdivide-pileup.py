#!/usr/bin/env python

import re
import shutil
from os import path
from pathlib import Path

def parse_coverage(file):
    return int(re.findall(r'(?<=[/])[0-9]+(?=[/])', file)[0])


def parse_child_filenames(parent_pileup):
    old_depth = parse_coverage(parent_pileup)

    parent_id = re.findall(r'(?<=[.])[0-9]+(?=[.]pileup)', parent_pileup)
    if len(parent_id) != 1:
        raise RuntimeError("Found multiple parent-id candidates")
    parent_id = int(parent_id[0])
    return [re.sub(f".{parent_id}.pileup$", f".{parent_id*2 + i}.pileup.tmp", parent_pileup) for i in range(2)]

def create_parent_directory(file):
    Path(path.dirname(file)).mkdir(parents=True, exist_ok=True)

def relocate_child(file, old_coverage):
    new_coverage = int(old_coverage/2)
    new_path     = re.sub(f"-{old_coverage}.", f"-{new_coverage}.", file)
    new_path     = re.sub(f"/{old_coverage}/", f"/{new_coverage}/", new_path)
    new_path     = path.splitext(new_path)[0]
    create_parent_directory(new_path)
    shutil.move(file, new_path)
    return new_path


def subdivide_pileup(file, maxdepth=10, depth=0):
    if depth >= maxdepth:
        return
    childs = parse_child_filenames(file)
    
    print(f"(depth {depth}): Subdividing {file}\n-> {childs[0]} | {childs[1]}")

    with open(file) as parent, open(childs[0], 'w') as left_child, open(childs[1], 'w') as right_child:
        for count, line in enumerate(parent):
            if not count % 2:
                left_child.write(line)
            else:
                right_child.write(line)

    old_coverage = parse_coverage(file)

    childs[0] = relocate_child(childs[0], old_coverage)
    childs[1] = relocate_child(childs[1], old_coverage)
    subdivide_pileup(childs[0], maxdepth=maxdepth, depth=depth+1)
    subdivide_pileup(childs[1], maxdepth=maxdepth, depth=depth+1)


if __name__ == '__main__':
    import sys
    subdivide_pileup(sys.argv[1], int(sys.argv[2]), depth=1)