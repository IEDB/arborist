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
  species_path = build_path / 'species'

  print('Combining all peptide and source assignments...\n')

  all_peptides_df, all_sources_df = pd.DataFrame(), pd.DataFrame()
  
  for species in species_path.iterdir():
    try:
      peptides_df = pd.read_csv(species / 'peptide-assignments.tsv', sep='\t')
      sources_df = pd.read_csv(species / 'source-assignments.tsv', sep='\t')
      species_data = pd.read_csv(species / 'species-data.tsv', sep='\t')
    except FileNotFoundError:
      continue

    peptides_df['Species Taxon ID'] = species.name
    peptides_df['Species Name'] = species_data['Species Name'].iloc[0]

    sources_df['Species Taxon ID'] = species.name
    sources_df['Species Name'] = species_data['Species Name'].iloc[0]
    sources_df['Proteome ID'] = species_data['Proteome ID'].iloc[0]
    sources_df['Proteome Label'] = species_data['Proteome Taxon'].iloc[0]


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