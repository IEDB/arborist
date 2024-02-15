#!/usr/bin/env python3

import os
import pandas as pd
from anytree import Node, RenderTree
from pathlib import Path


def build_tree_for_species(taxon_id: str) -> dict:
  """Create a tree of assigned genes and proteins for a given species.

  Args:
    taxon_id: The taxon ID of the species to build the tree for."""

  root = Node(f'Taxon ID: {taxon_id}')
  tree = {root.name: root}
  protein_tracker = {}

  try:
    peptide_assignments_df = pd.read_csv(build_path / 'species' / taxon_id / 'peptide-assignments.tsv', sep='\t')
    peptide_assignments_df['Parent Antigen Gene'].fillna(peptide_assignments_df['ARC Assignment'], inplace=True)
  except FileNotFoundError:
    return tree

  for _, row in peptide_assignments_df.iterrows():
    gene = row['Parent Antigen Gene'] if not pd.isna(row['Parent Antigen Gene']) else 'Unknown Gene'
    protein_name = row['Parent Antigen Gene Isoform Name']
    protein_id = row['Parent Antigen Gene Isoform ID']

    if pd.isna(protein_id):
      continue

    protein_info = f"{protein_name} (UniProt: {protein_id})" if not pd.isna(protein_name) else f"UniProt: {protein_id}"

    if gene not in tree:
      tree[gene] = Node(gene, parent=root)
      protein_tracker[gene] = set()

    if protein_info not in protein_tracker[gene]:
      Node(protein_info, parent=tree[gene])
      protein_tracker[gene].add(protein_info)

  return tree

def write_tree_to_file(tree:dict, outfile: str) -> None:
  """Write the tree to a file.

  Args:
    tree: The tree to write to file.
    outfile: The file to write the tree to."""

  with open(outfile, 'w') as f:
    for taxon_id, tree in tree.items():
      root_node = tree[f'Taxon ID: {taxon_id}']
      for pre, _, node in RenderTree(root_node):
        f.write("%s%s\n" % (pre, node.name))

if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument(
    '-b', '--build_path', 
    type=str,
    default=Path(__file__).parent.parent / 'build',
    help='Path to data directory.'
  )
  parser.add_argument(
    '-a', '--all_species', 
    action='store_true', 
    help='Build protein tree for all IEDB species.'
  )
  parser.add_argument(
    '-t', '--taxon_id', 
    type=str,
    help='Taxon ID for the species to pull data for.'
  )
  parser.add_argument(
    '-o', '--outfile', 
    type=str,
    help='Output file for tree.'
  )

  args = parser.parse_args()

  build_path = Path(args.build_path)
  all_species = args.all_species
  taxon_id = args.taxon_id
  outfile = args.outfile

  tree = {}
  if all_species:
    for taxon_id in os.listdir(build_path / 'species'):
      tree[taxon_id] = build_tree_for_species(taxon_id)

  else:
    assert taxon_id is not None, 'Must provide a taxon ID if not building all species.'
    tree[taxon_id] = build_tree_for_species(taxon_id)

  if outfile:
    write_tree_to_file(tree, outfile)
  else:
    for taxon_id, tree in tree.items():
      root_node = tree[f'Taxon ID: {taxon_id}']
      for pre, _, node in RenderTree(root_node):
        print("%s%s" % (pre, node.name))