#!/usr/bin/env python3

import pandas as pd
import sqlite3
import argparse

from pathlib import Path

import warnings
warnings.filterwarnings('ignore')



def build_old_tree(tree_df, peptide_assignments):
  new_rows = []
  for _, row in peptide_assignments.iterrows():
    new_rows.append(old_protein_class(row))
    new_rows.append(old_protein_label(row))
    # new_rows.append(old_gene_label(row))
  
  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)
  
  return tree_df

def old_protein_class(row):
  protein_class_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"UP:{row['Parent Antigen Gene Isoform ID']}",
    'predicate': 'rdfs:subClassOf',
    'object': f"NCBITaxon:{row['Species Taxon ID']}",
    'datatype': '_IRI',
    'annotation': None
  }
  return protein_class_row

def old_protein_label(row):
  protein_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"UP:{row['Parent Antigen Gene Isoform ID']}",
    'predicate': 'rdfs:label',
    'object': f"{row['Parent Antigen Gene Isoform Name']} (UniProt:{row['Parent Antigen Gene Isoform ID']})",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return protein_label_row

def old_gene_label(row):
  gene_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"UP:{row['Parent Antigen Gene Isoform ID']}",
    'predicate': 'from_gene',
    'object': f"{row['Parent Antigen Gene']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return gene_label_row

def build_new_tree(tree_df, peptide_assignments):
  new_rows = []
  for _, row in peptide_assignments.iterrows():
    new_rows.extend(new_gene_label(row))
    new_rows.append(new_protein_class(row))
    new_rows.append(new_protein_label(row))
    
  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)
  
  return tree_df

def new_gene_label(row):
  gene_class_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Species Taxon ID']}:{row['Parent Antigen Gene']}",
    'predicate': 'rdf:type',
    'object': 'owl:Class',
    'datatype': '_IRI',
    'annotation': None
  }
  gene_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Species Taxon ID']}:{row['Parent Antigen Gene']}",
    'predicate': 'rdfs:label',
    'object': f"{row['Parent Antigen Gene']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  gene_subclass_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Species Taxon ID']}:{row['Parent Antigen Gene']}",
    'predicate': 'rdfs:subClassOf',
    'object': f"NCBITaxon:{row['Species Taxon ID']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return gene_class_row, gene_label_row, gene_subclass_row

def new_protein_class(row):
  protein_class_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"UP:{row['Parent Antigen Gene Isoform ID']}",
    'predicate': 'rdfs:subClassOf',
    'object': f"{row['Species Taxon ID']}:{row['Parent Antigen Gene']}",
    'datatype': '_IRI',
    'annotation': None
  }
  return protein_class_row

def new_protein_label(row):
  protein_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"UP:{row['Parent Antigen Gene Isoform ID']}",
    'predicate': 'rdfs:label',
    'object': f"{row['Parent Antigen Gene Isoform Name']} (UniProt:{row['Parent Antigen Gene Isoform ID']})",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return protein_label_row

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
    tree_df = pd.read_sql_query("SELECT * FROM organism_tree", connection)
  
    predicates = ['rdfs:subClassOf', 'rdf:type', 'rdfs:label']
    lower_subjects = tree_df[tree_df['object'].isin(['upper', 'species'])]['subject']

    tree_df = tree_df[(tree_df['predicate'].isin(predicates))]
    tree_df = tree_df[tree_df['subject'].isin(lower_subjects)]
    tree_df['graph'] = 'iedb-taxon:protein_tree'
    tree_df.loc[tree_df['predicate'] == 'rdfs:label', 'object'] = tree_df['object'] + ' protein'
    tree_df.loc[tree_df['object'] == 'Organism protein', 'object'] = 'protein'
    
    old_df = build_old_tree(tree_df, peptide_assignments)
    new_df = build_new_tree(tree_df, peptide_assignments)
  
    old_df.to_sql('protein_tree_old', connection, if_exists='replace', index=False)
    new_df.to_sql('protein_tree_new', connection, if_exists='replace', index=False)

if __name__ == "__main__":
  main()
