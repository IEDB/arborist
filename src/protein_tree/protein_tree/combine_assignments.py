#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path


def main():
  """Combine all source antigen and epitope assignments into two TSV files."""

  print('Combining all epitope and source assignments...\n')

  all_epitopes_df, all_sources_df = pd.DataFrame(), pd.DataFrame()

  data_path = Path(__file__).parent.parent / 'data'
  for species in (data_path / 'species').iterdir():
    try:
      epitopes_df = pd.read_csv(species / 'epitope_assignments.tsv', sep='\t')
      sources_df = pd.read_csv(species / 'source_assignments.tsv', sep='\t')
    except FileNotFoundError:
      continue
    
    # add taxon ID and species name columns
    taxon_id = species.name.split('-')[0]
    species_name = species.name.split('-')[1]

    epitopes_df['Taxon ID'] = taxon_id
    epitopes_df['Species Name'] = species_name

    sources_df['Taxon ID'] = taxon_id
    sources_df['Species Name'] = species_name

    # rearrange columns to put taxon ID and species name first
    epitopes_df = epitopes_df[['Taxon ID', 'Species Name'] + epitopes_df.columns[:-2].tolist()]
    sources_df = sources_df[['Taxon ID', 'Species Name'] + sources_df.columns[:-2].tolist()]

    all_epitopes_df = pd.concat([all_epitopes_df, epitopes_df])
    all_sources_df = pd.concat([all_sources_df, sources_df])

  # write to files
  all_epitopes_df.to_csv(data_path / 'all_epitope_assignments.tsv', sep='\t', index=False)
  all_sources_df.to_csv(data_path / 'all_source_assignments.tsv', sep='\t', index=False)

  # report total successful assignments
  all_epitopes_df.drop_duplicates(subset=['Sequence', 'Source Accession'], inplace=True)
  all_sources_df.drop_duplicates(subset=['Accession'], inplace=True)

  source_count = len(all_sources_df)
  epitope_count = len(all_epitopes_df)

  source_success_count = len(all_sources_df[all_sources_df['Assigned Protein ID'].notna()])
  epitope_success_count = len(all_epitopes_df[all_epitopes_df['Assigned Protein ID'].notna()])

  print(f'Total epitopes: {epitope_count}')
  print(f'Total sources: {source_count}\n')

  print(f'Total successful epitope assignments: {epitope_success_count}')
  print(f'Total successful source assignments: {source_success_count}\n')

  print(f'Epitope success rate: {epitope_success_count / epitope_count * 100:.2f}%')
  print(f'Source success rate: {source_success_count / source_count * 100:.2f}%')


if __name__ == '__main__':
  main()