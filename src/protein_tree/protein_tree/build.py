#!/usr/bin/env python3


import json
import sqlite3
import argparse
import pandas as pd

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

  source_assignments['Assigned Protein Name'].fillna(source_assignments['Name'], inplace=True)
    
  new_rows = []
  nan_proteins = source_assignments['Assigned Protein ID'].isna()
  species_seen = set()

  # add parent entries as triplicate rows
  for _, row in source_assignments[~nan_proteins].iterrows():
    if row['ARC Assignment'] in ['TCR', 'BCR', 'MHC-I', 'MHC-II']:
      new_rows.extend(create_antigen_receptor_node(row, new_rows, species_seen))
      species_seen.add(row['Species Taxon ID'])
    else:
      new_rows.extend(
        owl_class(
          f"UP:{row['Assigned Protein ID']}",
          f"{row['Assigned Protein Name']} (UniProt:{row['Assigned Protein ID']})",
          f"iedb-protein:{row['Species Taxon ID']}"
      ))

    # fragments
    new_rows.extend(add_fragments(row))
    
    # annotations
    new_rows.extend(add_reviewed_status(row))
    new_rows.extend(add_synonyms(row))
    new_rows.extend(add_accession(row))
    new_rows.extend(add_source_database(row))

  new_rows.extend(create_other_nodes(source_assignments[nan_proteins]))

  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

  return tree_df


def create_antigen_receptor_node(source_assignment_row, new_rows, species_seen):
  """Given a row from the source_assignments dataframe that is an antigen receptor
  (TCR or BCR), return a list of triples for the antigen receptor node and create the node
  if it does not already exist."""

  antigen_receptor = source_assignment_row['ARC Assignment']
  antigen_receptor = 'ab' if antigen_receptor == 'BCR' else antigen_receptor

  if antigen_receptor == 'TCR':
    antigen_receptor_name = 'T Cell Receptor'
  elif antigen_receptor == 'ab':
    antigen_receptor_name = 'B Cell Receptor / Immunoglobulin'
  elif antigen_receptor == 'MHC-I':
    antigen_receptor_name = 'Major Histocompatibility Complex I'
  elif antigen_receptor == 'MHC-II':
    antigen_receptor_name = 'Major Histocompatibility Complex II'

  if source_assignment_row['Species Taxon ID'] not in species_seen:
    new_rows.extend(
      owl_class(
        f"iedb-protein:{source_assignment_row['Species Taxon ID']}-{antigen_receptor.lower()}",
        f"{antigen_receptor_name} chain",
        f"iedb-protein:{source_assignment_row['Species Taxon ID']}"
      )
    )
  
  prefix = 'UP' if source_assignment_row['Database'] == 'UniProt' else 'NCBI'
  assignment_node = owl_class(
    f"{prefix}:{source_assignment_row['Accession']}",
    f"{source_assignment_row['Name']} [{source_assignment_row['Accession']}]",
    f"iedb-protein:{source_assignment_row['Species Taxon ID']}-{antigen_receptor.lower()}"
  )

  return assignment_node


def add_fragments(row):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the fragments of the protein."""

  if pd.isna(row['Fragments']): return []
  
  with open(Path(__file__).parent.parent / 'data' / 'fragment-type.json', 'r') as f:
    fragment_type_map = json.load(f)
  
  fragment_count = len(row['Fragments'].split(', '))
  if fragment_count < 2: return []

  fragment_count = 0
  fragment_rows = []

  for fragment in row['Fragments'].split(', '):

    fragment_type = fragment.split('-')[0]
    fragment_type = fragment_type_map[fragment_type]
    fragment_start = fragment.split('-')[1]
    fragment_end = fragment.split('-')[2]

    try:
      int(fragment_start)
      int(fragment_end)
    except ValueError:
      continue

    fragment_rows.extend(owl_class(
      f"UP:{row['Assigned Protein ID']}-fragment-{fragment_count+1}",
      f"{fragment_type} ({fragment_start}-{fragment_end})",
      f"UP:{row['Assigned Protein ID']}"
    ))

    fragment_rows.extend(
      [triple(
        f"UP:{row['Assigned Protein ID']}-fragment-{fragment_count+1}",
        "ONTIE:0003627",
        fragment_start,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{row['Assigned Protein ID']}-fragment-{fragment_count+1}",
        "ONTIE:0003628",
        fragment_end,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{row['Assigned Protein ID']}-fragment-{fragment_count+1}",
        "ONTIE:0003620",
        f"{fragment_count+1} {fragment_type} ({fragment_start}-{fragment_end})",
        datatype="xsd:string"
      )]
    )

    fragment_count += 1

  return fragment_rows

def add_reviewed_status(row):
  """Given a row from the source_assignments dataframe, return a triple
  of the reviewed status of the protein."""
  return [triple(
    f"UP:{row['Assigned Protein ID']}",
    "UC:reviewed",
    "true" if row["Assigned Protein Reviewed"] == "sp" else "false",
    datatype="xsd:boolean"
  )]


def add_synonyms(row):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the synonyms of the protein."""
  synonyms = []
  if pd.isna(row['Synonyms']): return synonyms
  for synonym in row['Synonyms'].split(', '):
    if '@' in synonym or '{' in synonym: continue
    synonyms.append(triple(
      f"UP:{row['Assigned Protein ID']}",
      "ONTIE:0003622",
      synonym.split(' ')[0],
      datatype="xsd:string"
    ))
  return synonyms


def add_accession(row):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the accession of the protein and the URL to the UniProt entry."""
  return [triple(
    f"UP:{row['Assigned Protein ID']}",
    "ONTIE:0003623",
    row['Assigned Protein ID'],
    datatype="xsd:string"
  ),
  triple(
    f"UP:{row['Assigned Protein ID']}",
    "ONTIE:0003624",
    f"http://www.uniprot.org/uniprot/{row['Assigned Protein ID']}",
    datatype="xsd:string"
  )]


def add_source_database(row):
  """Given a row from the source_assignments dataframe, return a triple
  for the source database of the protein (always UniProt)."""
  return [triple(
    f"UP:{row['Assigned Protein ID']}",
    "ONTIE:0003625",
    "UniProt",
    datatype="xsd:string"
  )]


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
