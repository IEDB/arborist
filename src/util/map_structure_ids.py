#!/usr/bin/env python3

import argparse
import polars as pl

def map_structure_ids(structure_path, peptide_path):
  """Maps structure IDs to peptides using a two-step fallback logic with Polars."""
  structure_cols = ['structure_id', 'description', 'disc_region', 'mc_region', 'source_antigen_accession']
  structure = pl.read_csv(structure_path, separator='\t', has_header=False, new_columns=structure_cols)
  structure = structure.with_columns(
    pl.coalesce(['mc_region', 'disc_region', 'description']).alias('description')
  ).with_columns(
    pl.col('description').str.split(r'\s*\+\s*').list.get(0)
  )

  peptide = pl.read_csv(peptide_path, separator='\t')
  original_peptide_cols = peptide.columns

  map_with_acc = structure.filter(pl.col('source_antigen_accession').is_not_null())
  map_seq_only = structure.select(['description', 'structure_id']).unique(subset=['description'], keep='first')

  merged_df = peptide.join(
    map_with_acc,
    left_on=['Sequence', 'Source Accession'],
    right_on=['description', 'source_antigen_accession'],
    how='left'
  ).rename({'structure_id': 'id_from_acc'})

  merged_df = merged_df.join(
    map_seq_only,
    left_on='Sequence',
    right_on='description',
    how='left'
  ).rename({'structure_id': 'id_from_seq'})

  final_df = merged_df.with_columns(
    pl.coalesce(['id_from_acc', 'id_from_seq']).alias('Epitope ID')
  )
  
  return final_df.select(original_peptide_cols)

def main():
  argparser = argparse.ArgumentParser(description='Map structure IDs to peptides in peptide.tsv using Polars.')
  argparser.add_argument('structure', type=str, help='The input structure TSV file path.')
  argparser.add_argument('peptide', type=str, help='The input peptide TSV file path.')

  args = argparser.parse_args()

  peptide_mapped = map_structure_ids(args.structure, args.peptide)
  peptide_mapped.write_csv(args.peptide, separator='\t')

if __name__ == '__main__':
  main()