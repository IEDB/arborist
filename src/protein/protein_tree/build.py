import json
import sqlite3
import polars as pl
from pathlib import Path


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

def build_tree(tree_df, peptide_assignments):
  """Given dataframes for the tree and peptide_assignments,
  insert LDTab rows for each protein under its 'species protein' parent."""

  peptide_assignments['Parent Name'].fillna(peptide_assignments['Source Name'], inplace=True)

  new_rows = []
  species_seen = set()
  for parent_id, group in peptide_assignments.dropna(subset=['Parent ID']).groupby('Parent ID'):
    if group['ARC Assignment'].iloc[0] in ['TCR', 'BCR', 'MHC-I', 'MHC-II']:
      new_rows.extend(create_antigen_receptor_node(group, new_rows, species_seen))
      species_seen.add(group['Species Taxon ID'].iloc[0])
    else:
      new_rows.extend(
        owl_class(
          f"UP:{parent_id}",
          f"{group['Parent Name'].iloc[0]} (UniProt:{parent_id})",
          f"iedb-protein:{group['Species Taxon ID'].iloc[0]}"
      ))

    # fragments
    new_rows.extend(add_fragments(parent_id, group))
    
    # annotations
    new_rows.extend(add_reviewed_status(parent_id, group))
    new_rows.extend(add_synonyms(parent_id, group))
    new_rows.extend(add_accession(parent_id))
    new_rows.extend(add_source_database(parent_id))

  new_rows.extend(create_other_nodes(peptide_assignments[peptide_assignments['Parent ID'].isna()]))

  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

  return tree_df

def create_antigen_receptor_node(group, new_rows, species_seen):
  """Given a row from the source_assignments dataframe that is an antigen receptor
  (TCR or BCR), return a list of triples for the antigen receptor node and create the node
  if it does not already exist."""

  antigen_receptor = group['ARC Assignment'].iloc[0]
  antigen_receptor = 'ab' if antigen_receptor == 'BCR' else antigen_receptor

  if antigen_receptor == 'TCR':
    antigen_receptor_name = 'T Cell Receptor'
  elif antigen_receptor == 'ab':
    antigen_receptor_name = 'B Cell Receptor / Immunoglobulin'
  elif antigen_receptor == 'MHC-I':
    antigen_receptor_name = 'Major Histocompatibility Complex I'
  elif antigen_receptor == 'MHC-II':
    antigen_receptor_name = 'Major Histocompatibility Complex II'

  species_id = group['Species Taxon ID'].iloc[0]
  if species_id not in species_seen:
    new_rows.extend(
      owl_class(
        f"iedb-protein:{species_id}-{antigen_receptor.lower()}",
        f"{antigen_receptor_name} chain",
        f"iedb-protein:{species_id}"
      )
    )
  
  prefix = 'UP' if group['Source Database'].iloc[0] == 'UniProt' else 'NCBI'
  assignment_node = owl_class(
    f"{prefix}:{group['Source Accession'].iloc[0]}",
    f"{group['Source Name'].iloc[0]} [{group['Source Accession'].iloc[0]}]",
    f"iedb-protein:{species_id}-{antigen_receptor.lower()}"
  )

  return assignment_node

def add_fragments(parent_id, group):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the fragments of the protein."""

  if pd.isna(group['Fragments'].iloc[0]): return []
  
  with open(Path(__file__).parent.parent / 'data' / 'fragment-type.json', 'r') as f:
    fragment_type_map = json.load(f)
  
  fragment_count = len(group['Fragments'].iloc[0].split(', '))
  if fragment_count < 2: return []

  fragment_count = 0
  fragment_rows = []

  for fragment in group['Fragments'].iloc[0].split(', '):

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
      f"UP:{parent_id}-fragment-{fragment_count+1}",
      f"{fragment_type} ({fragment_start}-{fragment_end})",
      f"UP:{parent_id}"
    ))

    fragment_rows.extend(
      [triple(
        f"UP:{parent_id}-fragment-{fragment_count+1}",
        "ONTIE:0003627",
        fragment_start,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{parent_id}-fragment-{fragment_count+1}",
        "ONTIE:0003628",
        fragment_end,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{parent_id}-fragment-{fragment_count+1}",
        "ONTIE:0003620",
        f"{fragment_count+1} {fragment_type} ({fragment_start}-{fragment_end})",
        datatype="xsd:string"
      )]
    )

    fragment_count += 1

  return fragment_rows

def add_reviewed_status(parent_id, group):
  """Given a row from the source_assignments dataframe, return a triple
  of the reviewed status of the protein."""
  return [triple(
    f"UP:{parent_id}",
    "UC:reviewed",
    "true" if group["Assigned Protein Reviewed"].iloc[0] == "sp" else "false",
    datatype="xsd:boolean"
  )]

def add_synonyms(parent_id, group):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the synonyms of the protein."""
  synonyms = []
  if pd.isna(group['Source Synonyms'].iloc[0]): return synonyms
  for synonym in group['Source Synonyms'].iloc[0].split(', '):
    if '@' in synonym or '{' in synonym: continue
    synonyms.append(triple(
      f"UP:{parent_id}",
      "ONTIE:0003622",
      synonym.split(' ')[0],
      datatype="xsd:string"
    ))
  return synonyms

def add_accession(parent_id):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the accession of the protein and the URL to the UniProt entry."""
  return [triple(
    f"UP:{parent_id}",
    "ONTIE:0003623",
    parent_id,
    datatype="xsd:string"
  ),
  triple(
    f"UP:{parent_id}",
    "ONTIE:0003624",
    f"http://www.uniprot.org/uniprot/{parent_id}",
    datatype="xsd:string"
  )]

def add_source_database(parent_id):
  """Given a row from the source_assignments dataframe, return a triple
  for the source database of the protein (always UniProt)."""
  return [triple(
    f"UP:{parent_id}",
    "ONTIE:0003625",
    "UniProt",
    datatype="xsd:string"
  )]

def create_other_nodes(not_assigned):
  """Given a dataframe of sources without an assigned protein,
  return a list of triples for each 'Other' node for specific species."""
  
  new_rows = []
  species_with_nan = not_assigned['Species Taxon ID'].unique()
  species_id_to_name = {id: name for id, name in zip(not_assigned['Species Taxon ID'], not_assigned['Species Name'])}
  
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
  for source_accession, group in not_assigned.groupby('Source Accession'):
    new_rows.extend(
      owl_class(
        f"UP:{source_accession}" if group['Source Database'].iloc[0] == 'UniProt' else f"NCBI:{source_accession}",
        f"{group['Source Name'].iloc[0]} [{source_accession}]",
        f"iedb-protein:{group['Species Taxon ID'].iloc[0]}-other"
      )
    )
  
  return new_rows

if __name__ == "__main__":
  build_path = Path(__file__).parents[3] / 'build'

  # peptide_assignments = pl.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')

  with sqlite3.connect(build_path / 'arborist' / 'nanobot.db') as connection:
    # Copy the organism_tree, but replace each taxon with 'taxon protein'.
    tree_df = pl.read_database('''
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
      connection=connection, infer_schema_length=None
    )

    # Filter out subjects with iedb-taxon:level "lower" or blank.
    lower_subjects = tree_df.filter( pl.col('object').is_in(['lower', '']))['subject'].to_list()
    tree_df = tree_df.filter(~pl.col('subject').is_in(lower_subjects))


    exit(1)



    # Relabel 'taxon' to 'taxon protein'.
    tree_df.loc[(tree_df['subject'].str.startswith('iedb-protein:')) & (tree_df['predicate'] == 'rdfs:label'), 'object'] = tree_df['object'] + ' protein'

    # Add top-level 'protein'
    new_rows = owl_class('PR:000000001', 'protein', 'BFO:0000040')
    tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)

    # Re-parent children of 'Root' and 'organism' to 'protein'
    tree_df.loc[tree_df['object'] == 'iedb-protein:1', 'object'] = 'PR:000000001'
    tree_df.loc[tree_df['object'] == 'iedb-protein:28384', 'object'] = 'PR:000000001'
    tree_df.loc[tree_df['object'] == 'OBI:0100026', 'object'] = 'PR:000000001'

    old_df = build_old_tree(tree_df, peptide_assignments)
    new_df = build_new_tree(tree_df, peptide_assignments)
  
    old_df.to_sql('protein_tree_old', connection, if_exists='replace', index=False)
    new_df.to_sql('protein_tree_new', connection, if_exists='replace', index=False)