#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

import warnings
warnings.filterwarnings('ignore')


def combine_assignments(build_path):
  all_peptides_df, all_sources_df = pd.DataFrame(), pd.DataFrame()
  
  species_path = build_path / 'species'
  for species in species_path.iterdir():
    try:
      peptides_df = pd.read_csv(species / 'peptide-assignments.tsv', sep='\t')
      sources_df = pd.read_csv(species / 'source-assignments.tsv', sep='\t')
      species_data = pd.read_csv(species / 'species-data.tsv', sep='\t')
      proteome = pd.read_csv(species / 'proteome.tsv', sep='\t')
    except FileNotFoundError:
      continue
    
    proteome['Sequence Length'] = proteome['Sequence'].apply(len)
    id_to_seq_map = proteome.set_index('Protein ID')['Sequence'].to_dict()
    id_to_len_map = proteome.set_index('Protein ID')['Sequence Length'].to_dict()

    peptides_df['Species Taxon ID'] = species.name
    peptides_df['Species Name'] = species_data['Species Name'].iloc[0]

    sources_df['Species Taxon ID'] = species.name
    sources_df['Species Name'] = species_data['Species Name'].iloc[0]
    sources_df['Proteome ID'] = species_data['Proteome ID'].iloc[0]
    sources_df['Proteome Label'] = species_data['Proteome Taxon'].iloc[0]
    sources_df['Parent Sequence'] = sources_df['Assigned Protein ID'].map(id_to_seq_map)
    sources_df['Parent Sequence Length'] = sources_df['Assigned Protein ID'].map(id_to_len_map)

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


def generate_leidos_tables():
  def protein_strategy(row):
    if pd.notnull(row['ARC Assignment']):
      return 'antigen-receptor-classifier'
    elif row['Assignment Score'] > 80:
      return 'strong-blast-match'
    elif row['Assignment Score'] <= 80:
      return 'weak-blast-match'
    else:
      return 'manual'

  all_sources = pd.read_csv('build/arborist/all-source-assignments.tsv', sep='\t')
  all_sources['Protein Strategy'] = all_sources.apply(protein_strategy, axis=1)
  all_sources.loc[all_sources['Assigned Protein ID'].notnull(), 'Parent IRI'] = all_sources['Assigned Protein ID'].apply(lambda x: f'http://www.uniprot.org/uniprot/{x}')
  all_sources.loc[all_sources['Assigned Protein ID'].notnull(), 'Parent Protein Database'] = 'UniProt'

  source_parents = all_sources[[
    'Source ID', 'Accession', 'Database', 'Name', 'Aliases', 'Synonyms', 'Organism ID', 
    'Species Taxon ID', 'Species Name', 'Proteome ID', 'Proteome Label', 'Protein Strategy', 
    'Parent IRI', 'Parent Protein Database', 'Assigned Protein ID', 'Parent Sequence Length'
  ]]
  source_parents.rename(columns={
    'Organism ID': 'Taxon ID',
    'Species Taxon ID': 'Species ID',
    'Species Name': 'Species Label',
    'Assigned Protein ID': 'Parent Protein Accession',
  }, inplace=True)

  parent_proteins = all_sources[[
    'Assigned Protein ID', 'Parent Protein Database', 'Assigned Protein Name',
    'Proteome ID', 'Proteome Label', 'Parent Sequence'
  ]]
  parent_proteins.rename(columns={
    'Assigned Protein ID': 'Accession',
    'Parent Protein Database': 'Database',
    'Assigned Protein Name': 'Name',
    'Parent Sequence': 'Sequence'
  }, inplace=True)
  parent_proteins.drop_duplicates(subset=['Accession'], inplace=True)

  source_parents.to_csv('build/arborist/source-parents.tsv', sep='\t', index=False)
  parent_proteins.to_csv('build/arborist/parent-proteins.tsv', sep='\t', index=False)

def main():  
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
  combine_assignments(build_path)
  generate_leidos_tables()


if __name__ == '__main__':
  main()
