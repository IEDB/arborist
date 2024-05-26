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

  # make sure no extra spaces are in assigned protein ID column
  all_peptides_df['Parent ID'] = all_peptides_df['Parent ID'].str.split(' ').str[0]
  all_sources_df['Assigned Protein ID'] = all_sources_df['Assigned Protein ID'].str.split(' ').str[0]

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


def generate_leidos_tables(build_path):
  all_peptides_df = pd.read_csv('build/arborist/all-peptide-assignments.tsv', sep='\t')
  all_sources_df = pd.read_csv('build/arborist/all-source-assignments.tsv', sep='\t')
  
  generate_protein_tables(build_path, all_peptides_df, all_sources_df)
  generate_epitope_tables(build_path, all_peptides_df, all_sources_df)


def generate_protein_tables(build_path, all_peptides_df, all_sources_df):
  def protein_strategy(row):
    if pd.notnull(row['ARC Assignment_y']):
      return 'antigen-receptor-classifier'
    elif row['Source Assignment Score'] > 90:
      return 'strong-blast-match'
    elif row['Source Assignment Score'] <= 90:
      return 'weak-blast-match'
    else:
      return 'manual'
  
  df = pd.merge(all_peptides_df, all_sources_df, how='left', left_on='Source Accession', right_on='Accession')
  df['Protein Strategy'] = df.apply(protein_strategy, axis=1)
  df.loc[df['Parent ID'].notnull(), 'Parent IRI'] = df['Parent ID'].apply(lambda x: f'http://www.uniprot.org/uniprot/{x}')
  df.loc[df['Parent ID'].notnull(), 'Parent Protein Database'] = 'UniProt'

  source_parents = df[[
    'Source ID', 'Source Accession', 'Source Database', 'Source Name', 'Aliases', 'Synonyms', 'Organism ID_x', 
    'Species Taxon ID_x', 'Species Name_x', 'Proteome ID', 'Proteome Label', 'Protein Strategy', 
    'Parent IRI', 'Parent Protein Database', 'Parent ID', 'Parent Sequence Length', 'Assigned Gene'
  ]]
  source_parents.rename(columns={
    'Species Taxon ID_x': 'Species ID',
    'Species Name_x': 'Species Label',
    'Parent ID': 'Parent Protein Accession',
    'Assigned Gene': 'Parent Protein Gene'
  }, inplace=True)

  parent_proteins = df[[
    'Parent ID', 'Parent Protein Database', 'Parent Name',
    'Proteome ID', 'Proteome Label', 'Parent Sequence'
  ]]
  parent_proteins.rename(columns={
    'Parent ID': 'Accession',
    'Parent Protein Database': 'Database',
    'Parent Name': 'Name',
    'Parent Sequence': 'Sequence'
  }, inplace=True)

  source_parents.drop_duplicates(subset=['Source Accession'], inplace=True)
  source_parents.dropna(subset=['Parent IRI'], inplace=True)
  parent_proteins.drop_duplicates(subset=['Accession'], inplace=True)
  parent_proteins.dropna(subset=['Accession'], inplace=True)

  source_parents.to_csv(build_path / 'arborist' / 'source-parents.tsv', sep='\t', index=False)
  parent_proteins.to_csv(build_path / 'arborist' / 'parent-proteins.tsv', sep='\t', index=False)


def generate_epitope_tables(build_path, all_peptides_df, all_sources_df):
  epitope_mappings = all_peptides_df[[
    'Epitope ID', 'Sequence', 'Starting Position', 'Ending Position', 'Source Accession',
    'Parent Antigen ID', 'Parent Start', 'Parent End'
  ]]
  
  epitope_mappings['parent_seq'] = epitope_mappings['Parent Antigen ID'].map(all_sources_df.set_index('Assigned Protein ID')['Parent Sequence'].to_dict())
  epitope_mappings['identity_alignment'] = 1.0
  epitope_mappings['similarity_alignment'] = 1.0
  epitope_mappings['gaps_source_alignment'] = 0.0
  epitope_mappings['gaps_parent_alignment'] = 0.0
  epitope_mappings['all_gaps'] = 0.0

  epitope_mappings['source_alignment'] = epitope_mappings['Sequence']
  epitope_mappings['parent_alignment'] = epitope_mappings['Sequence']
  epitope_mappings['parent_alignment_modified'] = epitope_mappings['Sequence']

  epitope_mappings.rename(columns={
    'Epitope ID': 'epitope_id',
    'Sequence': 'epitope_seq',
    'Starting Position': 'epitope_start',
    'Ending Position': 'epitope_end',
    'Source Accession': 'source_accession',
    'Parent Antigen ID': 'parent_accession',
    'Parent Start': 'parent_start',
    'Parent End': 'parent_end'
  }, inplace=True)

  epitope_mappings = epitope_mappings[[
    'epitope_id', 'epitope_seq', 'epitope_start', 'epitope_end', 'source_accession',
    'parent_accession', 'parent_start', 'parent_end', 'identity_alignment',
    'similarity_alignment', 'gaps_source_alignment', 'gaps_parent_alignment',
    'all_gaps', 'source_alignment', 'parent_alignment', 'parent_alignment_modified'
  ]]
  epitope_mappings.to_csv(build_path / 'arborist' / 'epitope-mappings.tsv', sep='\t', index=False)


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
  generate_leidos_tables(build_path)

if __name__ == '__main__':
  main()
