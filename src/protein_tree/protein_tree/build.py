#!/usr/bin/env python3

import pandas as pd
import sqlite3
import argparse

from pathlib import Path

import warnings
warnings.filterwarnings('ignore')


def triple(subject, predicate, object, datatype='_IRI'):
  """Given subject, predicate, object, and optional datatype
  return a dictionary for an LDTab row in the protein_tree graph."""
  return {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': subject,
    'predicate': predicate,
    'object': object,
    'datatype': datatype,
    'annotation': None
  }


def owl_class(subject, label, parent):
  """Given a subject, label, and parent,
  return triples defining an owl:Class in the protein tree."""
  return [
    triple(subject, 'rdf:type', 'owl:Class'),
    triple(subject, 'rdfs:label', label, 'xsd:string'),
    triple(subject, 'rdfs:subClassOf', parent)
  ]


def build_old_tree(tree_df, source_assignments):
  """Given dataframes for the tree and peptide_assignments,
  insert LDTab rows for each protein under its 'species protein' parent."""
  
  new_rows = []
  nan_proteins = source_assignments['Assigned Protein ID'].isna()
  species_seen = set()

  # add parent entries as triplicate rows
  for _, row in source_assignments[~nan_proteins].iterrows():
    
    if row['ARC Assignment'] in ['TCR', 'BCR']:
      new_rows.extend(create_antigen_receptor_node(row, new_rows, species_seen))
      species_seen.add(row['Species Taxon ID'])
    else:
      new_rows.extend(
        owl_class(
          f"UP:{row['Assigned Protein ID']}",
          f"{row['Assigned Protein Name']} (UniProt:{row['Assigned Protein ID']})",
          f"iedb-protein:{row['Species Taxon ID']}"
      ))

      new_rows.extend(add_reviewed_status(row))
      # new_rows.extend(add_synonyms(row))
      # new_rows.extend(add_accession(row))
      # new_rows.extend(add_source_database(row))

  new_rows.extend(create_other_nodes(source_assignments[nan_proteins]))

  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

  return tree_df


def create_antigen_receptor_node(source_assignment_row, new_rows, species_seen):
  """Given a row from the source_assignments dataframe that is an antigen receptor
  (TCR or BCR), return a list of triples for the antigen receptor node and create the node
  if it does not already exist."""

  antigen_receptor = source_assignment_row['ARC Assignment']
  antigen_receptor_name = 'T Cell Receptor' if antigen_receptor == 'TCR' else 'B Cell Receptor / Immunoglobulin'

  if source_assignment_row['Species Taxon ID'] not in species_seen:
    new_rows.extend(
      owl_class(
        f"iedb-protein:{source_assignment_row['Species Taxon ID']}-{antigen_receptor.lower()}",
        f"{antigen_receptor_name} chain",
        f"iedb-protein:{source_assignment_row['Species Taxon ID']}"
      )
    )

  assignment_node = owl_class(
    f"{antigen_receptor}:{source_assignment_row['Accession']}",
    f"{source_assignment_row['Name']} [{source_assignment_row['Accession']}]",
    f"iedb-protein:{source_assignment_row['Species Taxon ID']}-{antigen_receptor.lower()}"
  )

  return assignment_node


def add_reviewed_status(row):
  return [triple(
    f"UP:{row['Assigned Protein ID']}",
    "UC:reviewed",
    "true" if row["Assigned Protein Reviewed"] == "sp" else "false",
    datatype="xsd:boolean"
  )]


def add_synonyms(row):
  pass


def add_accession(row):
  pass


def add_source_database(row):
  pass


def create_other_nodes(not_assigned_sources):
  """Given a dataframe of sources without an assigned protein,
  return a list of triples for each 'Other' node for specific species."""
  
  new_rows = []
  species_with_nan = not_assigned_sources['Species Taxon ID'].unique()
  species_id_to_name = {id: name for id, name in zip(not_assigned_sources['Species Taxon ID'], not_assigned_sources['Species Name'])}
  
  # create "Other" node for each certain species
  for species_id in species_with_nan:
    species_name = species_id_to_name[species_id]
    new_rows.extend(
      owl_class(
        f"iedb-protein:{species_id}-other",
        f"Other {species_name} protein",
        f"iedb-protein:{species_id}"
      )
    )

  # add proteins without a parent to the "Other" node
  for _, row in not_assigned_sources.iterrows():
    new_rows.extend(
      owl_class(
        f"UP:{row['Accession']}" if row['Database'] == 'UniProt' else f"NCBI:{row['Accession']}",
        f"{row['Name']} [{row['Accession']}]",
        f"iedb-protein:{row['Species Taxon ID']}-other"
      )
    )
  
  return new_rows


def build_new_tree(tree_df, peptide_assignments):
  """Given dataframes for the tree and peptide_assignments,
  insert LDTab rows for each gene under its 'species protein' parent,
  and each protein under its gene."""
 
  new_rows = []
  
  for _, row in peptide_assignments.iterrows():
    # WARN: This produces many duplicates?
    new_rows.extend(
      owl_class(
        f"{row['Species Taxon ID']}:{row['Parent Antigen Gene']}",
        f"{row['Parent Antigen Gene']}",
        f"iedb-protein:{row['Species Taxon ID']}"
      ))
    
    new_rows.extend(
      owl_class(
        f"UP:{row['Parent Antigen Gene Isoform ID']}",
        f"{row['Parent Antigen Gene Isoform Name']} (UniProt:{row['Parent Antigen Gene Isoform ID']})",
        f"{row['Species Taxon ID']}:{row['Parent Antigen Gene']}"
      ))

  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

  return tree_df


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

  # drop duplicates in source_assignments but keep NaNs
  source_assignments = pd.read_csv(build_path / 'arborist' / 'all-source-assignments.tsv', sep='\t')
  isna = source_assignments['Assigned Protein ID'].isna()
  source_assignments = pd.concat([
    source_assignments[isna],
    source_assignments[~isna].drop_duplicates(subset=['Assigned Protein ID'])
  ]).reset_index(drop=True)

  peptide_assignments = pd.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', sep='\t')
  peptide_assignments.drop_duplicates(subset=['Parent Antigen ID'], inplace=True)

  with sqlite3.connect(build_path / 'arborist' / 'nanobot.db') as connection:
    # Copy the organism_tree, but replace each taxon with 'taxon protein'.
    tree_df = pd.read_sql_query('''
      SELECT
        assertion,
        retraction,
        'iedb-taxon:protein_tree' AS graph,
        REPLACE(
          REPLACE(subject, 'iedb-taxon:', 'iedb-protein:'),
          'NCBITaxon:',
          'iedb-protein:'
        ) AS subject,
        predicate,
        REPLACE(
          REPLACE(object, 'iedb-taxon:', 'iedb-protein:'),
          'NCBITaxon:',
          'iedb-protein:'
        ) AS object,
        datatype,
        annotation
      FROM organism_tree
      WHERE subject NOT IN ('NCBITaxon:1', 'NCBITaxon:28384', 'OBI:0100026')''',
      connection
    )

    # Filter out subjects with iedb-taxon:level "lower" or blank.
    lower_subjects = tree_df[tree_df['object'].isin(['lower', ''])]['subject']
    tree_df = tree_df[-tree_df['subject'].isin(lower_subjects)]

    # Relabel 'taxon' to 'taxon protein'.
    tree_df.loc[(tree_df['subject'].str.startswith('iedb-protein:')) & (tree_df['predicate'] == 'rdfs:label'), 'object'] = tree_df['object'] + ' protein'

    # Add top-level 'protein'
    new_rows = owl_class('PR:000000001', 'protein', 'BFO:0000040')
    tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

    # Re-parent children of 'Root' and 'organism' to 'protein'
    tree_df.loc[tree_df['object'] == 'iedb-protein:1', 'object'] = 'PR:000000001'
    tree_df.loc[tree_df['object'] == 'iedb-protein:28384', 'object'] = 'PR:000000001'
    tree_df.loc[tree_df['object'] == 'OBI:0100026', 'object'] = 'PR:000000001'

    old_df = build_old_tree(tree_df, source_assignments)
    new_df = build_new_tree(tree_df, peptide_assignments)
  
    old_df.to_sql('protein_tree_old', connection, if_exists='replace', index=False)
    new_df.to_sql('protein_tree_new', connection, if_exists='replace', index=False)

if __name__ == "__main__":
  main()
