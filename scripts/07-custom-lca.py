#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import pytaxonkit as ptk

'''
Lowest common ancestor script for blast results of 16S and MiFish

Filtering done by maximum percent identity and maximum length

INPUT:
     Results from the script 06-blast-16S-MiFish.py (blast results)

RETURNS:
    Lowest common ancestor table based on percent identity and sequencing length, combined with taxonit lca call

Usage:
    python lca_16S_MiFish.py --blast_results [path to output from 16S/MiFish blast results] \
                             --dada2_asv_table [path to dada2 ASV count table output] \
                             --output [path to output file]


'''

def main():

    OUT_PATH="./"
    parser = argparse.ArgumentParser(description='Lowest Common Ancestor retrieval for blast results')
    parser.add_argument('--blast_results', metavar="BLAST_RES", type=str,
                        help='''Filename and path of blast output'''),
    parser.add_argument('--dada2_asv_table', metavar="DADA2_RES", type=str,
                       help='''Filename and path to dada2 ASV count table''')
    parser.add_argument('--output', metavar="OUTPUT", type=str, default=OUT_PATH,
                        help=f'Path to output file defaulting to "{OUT_PATH}".')

    args = parser.parse_args()
    blast_results = args.blast_results
    dada2_asv = args.dada2_asv_table
    output = args.output + "_LCA.csv"

    # Read the blast output and the dada2 ASV table
    dat = pd.read_csv(blast_results, sep = '\t')
    dada2 = pd.read_csv(dada2_asv)

    filtered_data = filter_pident_length(dat)
    assign_lineage(filtered_data)

    # Merge the LCA and the DADA2 results
    print(pd.merge(lca, dada2.set_index('sample_id').transpose().rename_axis('ASV').reset_index(), on='ASV', how='right').head())
    final_lca = pd.merge(lca, dada2.set_index('sample_id').transpose().rename_axis('ASV').reset_index(), on='ASV', how='right')

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


def assign_lineage(filtered_data):

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

    lca = pd.DataFrame()

    for i in filtered_data.ASV.unique():
        dat_slice = filtered_data[filtered_data.ASV ==i].copy()
        # Because some entries in the 16S and MiFish databases do not map to a taxid,
        # this results in some NaN entries. This if else statement will handle those and
        # skip it by adding NaN to the output
        if dat_slice.taxid.isna().all() == True:
            out = pd.DataFrame(dat_slice.taxid)
        else:
            out = pd.DataFrame([ptk.lca(dat_slice.taxid.tolist())])
            out = out.assign(ASV=i)

            lca = pd.concat([lca, out])

    # Need to reset the index, otherwise the additon of the lineage column does not work, as all indices in lca ar '0'
    lca = lca.reset_index(drop=True)
    lineage = ptk.lineage(lca[0].values.tolist())
    lca['LCA'] = lineage['Name']
    lca = lca.rename(columns = {0 : 'taxid'})
    #lca.colunms = ['taxid', 'ASV', 'LCA']
    lca = lca[['ASV', 'taxid', 'LCA']]
    #lineage = ptk.lineage(lca[[]])
    #lca.to_csv('lca_results_test.tsv', sep='\t', index = False)
    return lca


main()
