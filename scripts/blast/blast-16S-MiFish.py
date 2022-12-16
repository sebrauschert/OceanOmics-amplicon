#!/usr/bin/env python3
import argparse
import subprocess
import pandas as pd
import pytaxonkit

'''
Script to run blast with custom 16S and MiFish databases

This script will run blastn on the custom databases and return taxonomic annotations.
The script requires pandas and pytaxonkit in either an activated conda environment or installed on the system via
e.g. "pip install pandas" and "conda install -c bioconda pytaxonkit".

This pipeline contains the conda environmnent required for this script, which needs to be activated before running it.

It further requires a conda environment with blast installed, named 'blast-2.12.0'

Lastly, the taxonkit taxdump file will need to be installed:
mkdir ~/.taxonkit
cd ~/.taxonkit
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz

INPUT:
      - fasta file with ASV sequences from DADA2, containing the ASV IDs and sequences
      - output file name/path
      - database name either 16S or MiFish


RETURNS:
*******************************************
Blast against custom 16S/MiFish database finished!
Output file : Output_file_name.tsv
Glimpse at it:
*******************************************

      ASV                  Species   taxid  pident  length  mismatch  gaps  qstart  qend  sstart  send        evalue  qcovs  qcovhsp
0  ASV_64       Parupeneus indicus  334904  88.298     188        20     1       1   186     201    14  2.920000e-61    100      100
1  ASV_64       Parupeneus indicus  334904  88.298     188        20     1       1   186     201    14  2.920000e-61    100      100
2  ASV_64  Parupeneus heptacanthus  334900  89.941     169        17     0       1   169     200    32  1.360000e-59     91       91
3  ASV_64      Priacanthus tayenus  443711  78.659     164        30     5       1   161     202    41  3.970000e-25     87       87
4  ASV_64      Priacanthus tayenus  443711  78.659     164        30     5       1   161     202    41  3.970000e-25     87       87


Usage:
     python 06-blast-16S-MiFish.py --dada2_file [DADA2 fasta file] \
                                   --out_path [filename or path for output] \
                                   --database [either '16S' or 'MiFish']
'''

#............................................................................
# Code to activate conda environments
CONDA_BLAST    = 'eval "$(conda shell.bash hook)" ; conda activate blast-2.12.0'

# Locations of 16S and MiFish databases
DB16S  = "/data/tools/databases/16SDB/16S_MB_database.fasta"
MIFISH = "/data/tools/databases/MiFishDB/MiFishDB.fasta"
#............................................................................

def main():

    OUT_PATH="./"
    parser = argparse.ArgumentParser(description='Run blastn on a query ASV sequence from a DADA2 run and return taxonomic annotation.')
    parser.add_argument('--dada2_file', metavar='FILE', type=str,
                        help='''Filename of dada2 fasta file with query sequences of ASVs''')
    parser.add_argument('--out_path', metavar='DIR', type=str, default=OUT_PATH,
                        help=f'Output filename, defaulting to "{OUT_PATH}".')
    parser.add_argument('--database', metavar='DB', type=str, choices=['16S', 'MiFish'],
                        help='''Select either the 16S or the MiFish database.''')
    args = parser.parse_args()
    dada2_file = args.dada2_file
    out_path = args.out_path
    database = args.database

    OUT_NAME = out_path + database + '_blast_results.tsv'


    # Blast either against the 16S or MiFish database and return formatted file
    if database == '16S':
        # Run blastn for the 16S database
        blastn(dada2_file, OUT_NAME, database)
        # Process the result file to retrieve the correct format
        process_blast_output(OUT_NAME, database)
    else:
        # Run blastn for the MiFish database
        blastn(dada2_file, OUT_NAME, database)
        # Process the result file to retrieve the correct format
        process_blast_output(OUT_NAME, database)



# Define function to run blast on 16S and MiFish database

def blastn(dada2_file, out_file, database):

    '''
    blastn function

    Function to activate the blast conda environment and blast a query sequence
    against a custom 16S and the MiFish database http://mitofish.aori.u-tokyo.ac.jp/download.html

    Parameters:
    dada2_fasta (path): query amplicon sequences from the DADA2 R script
    out_file (string) : name and path for output file
    '''
    if database == '16S':
        blast_command = 'blastn -db ' + DB16S + ' -query ' + dada2_file + ' -outfmt "6 qseqid sseqid staxids stitle scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -num_threads 200 > ' + out_file
        subprocess.call(CONDA_BLAST + ';' + blast_command, shell = True)

    else:
        blast_command = 'blastn -db ' + MIFISH + ' -query ' + dada2_file + ' -outfmt "6 qseqid sseqid staxids stitle scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -num_threads 200 > ' + out_file
        subprocess.call(CONDA_BLAST + ';' + blast_command, shell = True)



# Prepare the output from the 16S database blast for taxonkit

def process_blast_output(out_file, database):

    '''
    Internal function to prepare the 16S query file for taxonkit annotation

    Parameters:
    out_file (string): taken from the command line, the blastn output file name taken as input for this function
    '''

    if database == '16S':

        blast_16 = pd.read_csv(out_file, sep = '\t', header = None)
        blast_16[blast_16.columns[3]] = blast_16[blast_16.columns[3]].str.split(',').str[0].values
        # VERY few species have two taxonomy IDs: Guyu wujalwujalensis, Xiphocheilus typus
        # we just choose the first one here
        taxid_df = pytaxonkit.name2taxid(blast_16[blast_16.columns[3]])
        # we can have duplicate IDs, so taxid_df is longer than the blast_mifish df
        # first, remove the second/third ones
        unique_taxids = taxid_df.groupby('Name').first().reset_index()
        # now remove duplicate taxids
        taxid_df = taxid_df.loc[taxid_df['TaxID'].isin(unique_taxids['TaxID'])].reset_index()
        # NOT resetting the index leads to curious bugs when adding the column below
        assert taxid_df.shape[0] == blast_16.shape[0], \
                f'ERROR: the pytaxonkit results have different rows ({taxid_df.shape[0]}) than the blast results ({blast_16.shape[0]})'
        blast_16['taxid'] = taxid_df['TaxID']
 
        blast_16 = blast_16[[0, 1, 'taxid', 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]
        blast_16.columns = ['ASV', 'sseqid', 'taxid', 'Species',  'scomnames', 'sskingdoms', 'pident', 'length', 'qlen', 'slen', 'mismatch', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', '20']
        blast_16.to_csv(out_file, sep = '\t', index = False)
        print('')
        print('*******************************************')
        print('Blast against custom 16S database finished!')
        print('Output file : ' + out_file)
        print('Glimpse at it:')
        print('*******************************************')
        print('')

        print(blast_16.head())

    else:

        taxa_list = []
        blast_mifish = pd.read_csv(out_file, sep = '\t', header = None)

        for hit in blast_mifish[blast_mifish.columns[3]]:

            if hit.startswith('gi'):
                tax = hit.split('|')[6].strip()
                taxa_list.append(tax)

            else:
                tax = hit.split('|')[2].split('([')[0].strip()
                taxa_list.append(tax)

        blast_mifish['taxa'] = taxa_list
        taxid_df = pytaxonkit.name2taxid(blast_mifish[blast_mifish.columns[len(blast_mifish.columns)-1]])
        unique_taxids = taxid_df.groupby('Name').first().reset_index()
        taxid_df = taxid_df.loc[taxid_df['TaxID'].isin(unique_taxids['TaxID'])].reset_index()
        assert taxid_df.shape[0] == blast_mifish.shape[0], \
                f'ERROR: the pytaxonkit results have different rows ({taxid_df.shape[0]}) than the blast results ({blast_mifish.shape[0]})'
        blast_mifish['taxid'] = taxid_df['TaxID']
        
        blast_mifish = blast_mifish[[0, 1, 'taxid', 'taxa', 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]
        blast_mifish.columns = ['ASV', 'sseqid', 'taxid', 'Species',  'scomnames', 'sskingdoms', 'pident', 'length', 'qlen', 'slen', 'mismatch', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', '20']
        blast_mifish.to_csv(out_file, sep = '\t', index = False)

        print('')
        print('*******************************************')
        print('Blast against custom MiFish database finished!')
        print('Output file : ' + out_file)
        print('Glimpse at it:')
        print('*******************************************')
        print('')

        print(blast_mifish.head())
main()
