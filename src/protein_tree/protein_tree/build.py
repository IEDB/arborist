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


def build_old_tree(tree_df, peptide_assignments):
  """Given dataframes for the tree and peptide_assignments,
  insert LDTab rows for each protein under its 'species protein' parent."""
  new_rows = []
  for _, row in peptide_assignments.iterrows():
    new_rows.extend(
      owl_class(
        f"UP:{row['Parent Antigen ID']}",
        f"{row['Parent Antigen Name']} (UniProt:{row['Parent Antigen ID']})",
        f"iedb-protein:{row['Species Taxon ID']}"
    ))

  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

  return tree_df


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

  peptide_assignments = pd.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', sep='\t')
  peptide_assignments.drop_duplicates(subset=['Parent Antigen Gene Isoform ID'], inplace=True)

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
      WHERE subject NOT IN ('NCBITaxon:1', 'OBI:0100026')''',
      connection
    )

    # Filter out subjects with iedb-taxon:level "lower"
    lower_subjects = tree_df[tree_df['object'] == 'lower']['subject']
    tree_df = tree_df[-tree_df['subject'].isin(lower_subjects)]

    # Relabel 'taxon' to 'taxon protein'.
    tree_df.loc[(tree_df['subject'].str.startswith('iedb-protein:')) & (tree_df['predicate'] == 'rdfs:label'), 'object'] = tree_df['object'] + ' protein'

    # Add top-level 'protein'
    new_rows = owl_class('PR:000000001', 'protein', 'owl:Thing')
    tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

    # Re-parent children of 'organism' to 'protein'
    tree_df.loc[tree_df['object'] == 'OBI:0100026', 'object'] = 'PR:000000001'

    old_df = build_old_tree(tree_df, peptide_assignments)
    new_df = build_new_tree(tree_df, peptide_assignments)
  
    old_df.to_sql('protein_tree_old', connection, if_exists='replace', index=False)
    new_df.to_sql('protein_tree_new', connection, if_exists='replace', index=False)

if __name__ == "__main__":
  main()
