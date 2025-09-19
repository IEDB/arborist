#!/usr/bin/env python3

import argparse
import pandas as pd

import warnings
warnings.filterwarnings('ignore')


def map_structure_ids(structure_path, peptide_path):
    structure_cols = ['structure_id', 'description', 'disc_region', 'mc_region', 'source_antigen_accession']
    structure = pd.read_csv(structure_path, sep='\t', names=structure_cols)
    peptide = pd.read_csv(peptide_path, sep='\t')

    # get the full regions for certain structures that have it
    structure['description'] = structure['mc_region'].combine_first(structure['disc_region']).combine_first(structure['description'])

    # separate PTM from description
    structure['description'] = structure['description'].str.split(r'\s*\+\s*').str[0]

    structure_with_acc = structure.dropna(subset=['source_antigen_accession'])
    structure_without_acc = structure[structure['source_antigen_accession'].isnull()]

    # Try to map the structure IDs using both sequence and source accession
    peptide['map_key'] = peptide['Sequence'].astype(str) + '|' + peptide['Source Accession'].astype(str)
    composite_map_series = (
        structure_with_acc['description'].astype(str) + '|' + structure_with_acc['source_antigen_accession'].astype(str)
    )
    structure_map_acc = pd.Series(
        structure_with_acc['structure_id'].values, index=composite_map_series
    ).to_dict()

    peptide['Epitope ID'] = peptide['map_key'].map(structure_map_acc)

    # Create a second map using only sequences from structures that LACK an accession
    structure_map_seq_only = structure_without_acc.drop_duplicates(
        subset=['description']
    ).set_index('description')['structure_id'].to_dict()
    fallback_matches = peptide['Sequence'].map(structure_map_seq_only)
    peptide['Epitope ID'] = peptide['Epitope ID'].combine_first(fallback_matches)

    peptide = peptide.drop(columns=['map_key'])
    return peptide


def main():
    argparser = argparse.ArgumentParser(description='Map structure IDs to peptides in peptide.tsv.')
    argparser.add_argument('structure', type=str, help='The input structure TSV file path.')
    argparser.add_argument('peptide', type=str, help='The input peptide TSV file path.')

    args = argparser.parse_args()

    peptide = map_structure_ids(args.structure, args.peptide)
    peptide.to_csv(args.peptide, sep='\t', index=False)


if __name__ == '__main__':
    main()