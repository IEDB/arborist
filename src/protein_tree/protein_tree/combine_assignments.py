#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path


def main():
  """Combine all peptide and peptide source assignments into two TSV files."""
  
  parser = argparse.ArgumentParser()

  parser.add_argument(
    'build_path',
    type=str,
    default='build/',
    help='Path to species directory.'
  )

  args = parser.parse_args()
  build_path = Path(args.build_path)

  print('Combining all peptide and source assignments...\n')

  active_species = pd.read_csv(build_path / 'arborist' / 'active-species.tsv', sep='\t')
  all_peptides_df, all_sources_df = pd.DataFrame(), pd.DataFrame()
  
  for species in active_species['Species ID']:
    species_path = build_path / 'species' / f'{species}'

    try:
      peptides_df = pd.read_csv(species_path / 'peptide-assignments.tsv', sep='\t')
      sources_df = pd.read_csv(species_path / 'source-assignments.tsv', sep='\t')
      species_data = pd.read_csv(species_path / 'species-data.tsv', sep='\t')
    except FileNotFoundError:
      continue

    peptides_df['Species Taxon ID'] = species
    peptides_df['Species Name'] = species_data['Species Name'].iloc[0]

    sources_df['Species Taxon ID'] = species
    sources_df['Species Name'] = species_data['Species Name'].iloc[0]

    # rearrange columns to put taxon ID and species name first
    peptides_df = peptides_df[['Species Taxon ID', 'Species Name'] + peptides_df.columns[:-2].tolist()]
    sources_df = sources_df[['Species Taxon ID', 'Species Name'] + sources_df.columns[:-2].tolist()]

    all_peptides_df = pd.concat([all_peptides_df, peptides_df])
    all_sources_df = pd.concat([all_sources_df, sources_df])

  # fill in missing gene values with ARC assignments
  all_peptides_df['Parent Antigen Gene'].fillna(all_peptides_df['ARC Assignment'], inplace=True)
  all_sources_df['Assigned Gene'].fillna(all_sources_df['ARC Assignment'], inplace=True)

  # write to files
  all_peptides_df.to_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', sep='\t', index=False)
  all_sources_df.to_csv(build_path / 'arborist' / 'all-source-assignments.tsv', sep='\t', index=False)

  # report total successful assignments
  all_peptides_df.drop_duplicates(subset=['Sequence', 'Source Accession'], inplace=True)
  all_sources_df.drop_duplicates(subset=['Accession'], inplace=True)

  source_count = len(all_sources_df)
  peptide_count = len(all_peptides_df)

  source_success_count = len(all_sources_df[all_sources_df['Assigned Protein ID'].notna()])
  peptide_success_count = len(all_peptides_df[all_peptides_df['Parent Antigen Gene Isoform ID'].notna()])

  print(f'Total peptides: {peptide_count}')
  print(f'Total sources: {source_count}\n')

  print(f'Total successful peptide assignments: {peptide_success_count}')
  print(f'Total successful source assignments: {source_success_count}\n')

  print(f'Peptide success rate: {peptide_success_count / peptide_count * 100:.2f}%')
  print(f'Source success rate: {source_success_count / source_count * 100:.2f}%')


if __name__ == '__main__':
  main()
