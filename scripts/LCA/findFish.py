'''
Fish taxonomy IDs:

    EITHER

The clade of vertebrata 7742 (will need filtering for human and other animal contamination like seasnakes)

    OR

all fish, which includes

1476529 (Cyclostomata, jawless vertebrates - we should rarely have hits for this, it's only hagfish and lampreys)
7777 (Chondrichthyes, cartilaginous fishes which includes sharks)
Bony fishes don't have a taxonomy ID, but their subgroups have:
7898 (Actinopterygii, ray-finned fishes)
8287 (Sarcopterygii, lobe-finned fishes but it also includes all animals! So need to exclude taxonomy ID 32523, Tetrapoda)
There's also this weirdo 'outlier' taxonomy ID: 1476750  (fish environmental sample) which should probably be included

So 1476529 + 7777 + 7898 + 8287 + 1476750 - 32523 to pull out only fish from taxonomy IDs

TODO: Currently this script just checks for family/clade names. It would make more sense to do it properly,
      i.e. checking whether the given taxonomy ID is 'below' the above taxonomy IDs. Surprisingly hard using taxonkit
      
Currently, the script prints out a line for each lineage that contains fish.
'''

import argparse
import logging

def main():
    parser = argparse.ArgumentParser(description='Check whether computeLCA.py output contains fish.')
    parser.add_argument('input', help='Path to input file from computeLCA.py')
    args = parser.parse_args()

    with open(args.input) as fh:
        for line in fh:
            ll = line.rstrip().split('\t')
            try:
                asv, taxid, lineage, tax = ll
            except ValueError:
                # some lines are too short
                logging.warning(f'{ll[0]} has no taxa?')
                continue

            if 'Cyclostomata' in lineage or 'Chondrichthyes' in lineage \
                    or 'Sarcopterygii' in lineage or 'Actinopterygii' in lineage or \
                    'fish environmental' in lineage:
                print(line.rstrip())

main()
