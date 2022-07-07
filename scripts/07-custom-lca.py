#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import pytaxonkit as ptk
from joblib import Parallel, delayed

'''
Lowest common ancestor script for blast results of 16S and MiFish

Filtering done by maximum percent identity and maximum length

INPUT:
     Results from the script 06-blast-16S-MiFish.py (blast results)

RETURNS:
    Lowest common ancestor table based on percent identity and sequencing length, combined with taxonit lca call

Usage:
    python 07-custom-lca.py --blast_results [path to output from 16S/MiFish blast results] \
                             --dada2_asv_table [path to dada2 ASV count table output] \
                             --cores [number of cores for parallel LCA] \
                             --output [path to output file]


'''

def main():

    OUT_PATH="./"
    parser = argparse.ArgumentParser(description='Lowest Common Ancestor retrieval for blast results')
    parser.add_argument('--blast_results', metavar="BLAST_RES", type=str,
                        help='''Filename and path of blast output'''),
    parser.add_argument('--dada2_asv_table', metavar="DADA2_RES", type=str,
                       help='''Filename and path to dada2 ASV count table'''),
    parser.add_argument('--cores', metavar="CORES", type=int, default=100,
                       help='''Number of cores for parallel LCA''')
    parser.add_argument('--output', metavar="OUTPUT", type=str, default=OUT_PATH,
                        help=f'Path to output file defaulting to "{OUT_PATH}".')

    
    args = parser.parse_args()
    blast_results = args.blast_results
    dada2_asv = args.dada2_asv_table
    cores = args.cores
    output = args.output + "_LCA.csv"
    

    # Read the blast output and the dada2 ASV table
    dat = pd.read_csv(blast_results, sep = '\t')
    dada2 = pd.read_csv(dada2_asv, sep = '\t')

    filtered_data = filter_pident_length(dat)
    lca = assign_lineage(filtered_data, cores, dat)

    # Merge the LCA and the DADA2 results
    # We perform a left join, as we are only interested in reteining those ASV entries that have an LCA assigned.
    print(pd.merge(lca, dada2, on='ASV', how='left').head())
    final_lca = pd.merge(lca, dada2, on='ASV', how='left')
    
    final_lca.drop(final_lca.index[final_lca['ASV'].isnull()], inplace=True)
    final_lca = final_lca.reset_index(drop=True)
    
    final_lca.to_csv(output, index = False)



def filter_pident_length(data):
    '''
    Filter for highest percent identity of query and database sequence as well as sequence length

    For the LCA, we need to filter the blast results; this is a rough default filter that only retains
    the entries per ASV that have the highest pident and length entries

    Parameters:
    data: data set loaded from the path provided to the main function

    '''
    print('**************************')
    print('Filtering blast results...')
    print('**************************')

    lca = pd.DataFrame()
    for i in data.ASV.unique():

        dat_slice = data[data.ASV == i].copy()
        dat_slice['filter'] = dat_slice['pident'] + dat_slice['length']
        output = dat_slice[dat_slice['filter'] == dat_slice['filter'].max()]
        lca = pd.concat([lca, output])


    return lca.reset_index(drop=True)


def lca_runner(i, filtered_data):
    '''
    Function to generalise the creation of the lca output, so we can parallelize it.
    Internal function, not to be called by itself
    '''
    
    dat_slice = filtered_data[filtered_data.ASV == i].copy()
            # Because some entries in the 16S and MiFish databases do not map to a taxid,
            # this results in some NaN entries. This if else statement will handle those and
            # skip it by adding NaN to the output
    if dat_slice.taxid.isna().all() == True:
        out = dat_slice.taxid
    else:
        out = ptk.lca(dat_slice.taxid.tolist())
        out = [out, i]

    return out




def assign_lineage(filtered_data, cores, dat):

    '''

    '''

    print('**************************')
    print('Taxonkit lca algorithm  ..')
    print('**************************')

    filtered_data['taxid'] = filtered_data['taxid'].fillna(-1)
    filtered_data['taxid'] = filtered_data['taxid'].astype(int)
    filtered_data['taxid'] = filtered_data['taxid'].astype(str)
    filtered_data['taxid'] = filtered_data['taxid'].replace('-1', np.nan)

    #lineage = ptk.lineage(filtered_data.taxid.tolist())
    #filtered_data['FullLineage'] = lineage[['FullLineage']]

    lca = Parallel(n_jobs=cores)(delayed(lca_runner)(i, filtered_data) for i in dat.ASV.unique())
    lca = pd.DataFrame(columns=['taxid', 'ASV'], data=lca)


    lca['taxid'] = lca['taxid'].fillna(-1)
    lca['taxid'] = lca['taxid'].astype(int)
    lca['taxid'] = lca['taxid'].astype(str)
    lca['taxid'] = lca['taxid'].replace('-1', np.nan)
    
    # Need to reset the index, otherwise the additon of the lineage column does not work, as all indices in lca ar '0'
    #lca = lca.reset_index(drop=True)
    lineage = ptk.lineage(lca['taxid'].values.tolist())
    lca['LCA'] = lineage['Name']
    lca = lca.rename(columns = {0 : 'taxid'})
    #lca.colunms = ['taxid', 'ASV', 'LCA']
    lca = lca[['ASV', 'taxid', 'LCA']]
    #lineage = ptk.lineage(lca[[]])
    #lca.to_csv('lca_results_test.tsv', sep='\t', index = False)
    return lca


main()
