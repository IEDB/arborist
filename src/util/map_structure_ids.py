#!/usr/bin/env python3

import argparse
import pandas as pd

import warnings
warnings.filterwarnings('ignore')


def map_structure_ids(structure_path, peptide_path):
    structure = pd.read_csv(structure_path, sep='\t')
    peptide = pd.read_csv(peptide_path, sep='\t')

    # get the full regions for certain structures that have it
    structure['description'] = structure['mc_region'].combine_first(structure['disc_region']).combine_first(structure['description'])

    # separate PTM from description
    structure['description'] = structure['description'].str.split(r'\s*\+\s*').str[0]

    structure_map = structure.set_index('description')['structure_id'].to_dict()
    peptide['Epitope ID'] = peptide['Sequence'].map(structure_map)
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