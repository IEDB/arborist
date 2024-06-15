import json
import sqlite3
import polars as pl
from pathlib import Path


def build_tree(tree_base, assignments):
  """Given the tree_base from the organism tree dataframes for the tree and 
  the assignments, insert LDTab rows for each protein under its 
  'species protein' parent."""

  new_rows = []

  arc_parents = assignments.filter(
    (pl.col('ARC Assignment').str.contains('TCR')) | (pl.col('ARC Assignment').str.contains('BCR')) | (pl.col('ARC Assignment').str.contains('MHC')),
  ).unique(subset=['Source Accession'])

  normal_parents = assignments.join(
    assignments.group_by('Source Accession').agg(
      pl.col('Assigned Protein ID').mode().first()
    ).rename({'Assigned Protein ID': 'Parent Protein ID'}), 
    how='left', on='Source Accession', coalesce=True
  ).filter(
    (pl.col('Assigned Protein ID') == pl.col('Parent Protein ID')),
    ~(pl.col('Source Accession').is_in(arc_parents['Source Accession']))
  ).unique(subset=['Parent Protein ID'])

  others = assignments.filter(
    (pl.col('Assigned Protein ID').is_null()),
  ).unique(subset=['Source Accession'])
  
  new_rows.extend(add_normal_parents(normal_parents))
  new_rows.extend(add_arc_parents(arc_parents))
  new_rows.extend(add_others(others))

  protein_tree = pl.concat([tree_base, pl.DataFrame(new_rows)])

  return protein_tree

def add_normal_parents(normal_parents):
  """Proteins that have been assigned a parent."""
  rows = []
  for parent in normal_parents.iter_rows(named=True):
    rows.extend(
      owl_class(
        f"UP:{parent['Parent Protein ID']}",
        f"{parent['Assigned Protein Name']} (UniProt:{parent['Parent Protein ID']})",
        f"iedb-protein:{parent['Species Taxon ID']}"
      )
    )
    rows.extend(add_metadata(parent))
  return rows

def add_arc_parents(arc_parents):
  """Proteins that have been assigned an antigen receptor via ARC."""
  rows = []
  nodes_seen = set()
  arc_parents = arc_parents.rename({'Assigned Protein ID': 'Parent Protein ID'})
  for parent in arc_parents.iter_rows(named=True):
    arc_assignment = parent['ARC Assignment']
    antigen_receptor_class = arc_assignment.split('_')[0]
    node_id = f"{str(parent['Species Taxon ID'])}-{antigen_receptor_class.lower()}"
    node_name = antigen_receptor_name(antigen_receptor_class)
    if node_id not in nodes_seen:
      rows.extend(
        owl_class(
          f"iedb-protein:{node_id}",
          f"{node_name} chain",
          f"iedb-protein:{parent['Species Taxon ID']}"
        )
      )
      nodes_seen.add(node_id)

    rows.extend(
      owl_class(
        f"UP:{parent['Parent Protein ID']}",
        f"{parent['Assigned Protein Name']} (UniProt:{parent['Parent Protein ID']})",
        f"iedb-protein:{node_id}"
      )
    )
    rows.extend(add_metadata(parent))
  return rows

def antigen_receptor_name(antigen_receptor_class):
  if antigen_receptor_class == 'TCR':
    return 'T Cell Receptor'
  elif antigen_receptor_class == 'BCR':
    return 'B Cell Receptor / Immunoglobulin'
  elif antigen_receptor_class == 'MHC-I':
    return 'Major Histocompatibility Complex I'
  elif antigen_receptor_class == 'MHC-II':
    return 'Major Histocompatibility Complex II'

def add_others(others):
  """Proteins that have not been assigned a parent nor as an antigen receptor."""
  rows = []
  species_seen = set()
  for other in others.iter_rows(named=True):
    if other['Species Taxon ID'] not in species_seen:
      rows.extend(
        owl_class(
          f"iedb-protein:{other['Species Taxon ID']}-other",
          f"Other {other['Species Name']} protein",
          f"iedb-protein:{other['Species Taxon ID']}"
        )
      )
      species_seen.add(other['Species Taxon ID'])

    rows.extend(
      owl_class(
        f"UP:{other['Source Accession']}" if other['Database'] == 'UniProt' else f"NCBI:{other['Source Accession']}",
        f"{other['Name']} [{other['Source Accession']}]",
        f"iedb-protein:{other['Species Taxon ID']}-other"
      )
    )
  return rows

def add_metadata(parent):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the metadata of the protein."""
  metadata_rows = []
  metadata_rows.extend(add_fragments(parent))
  metadata_rows.extend(add_reviewed_status(parent))
  metadata_rows.extend(add_synonyms(parent))
  metadata_rows.extend(add_accession(parent))
  metadata_rows.extend(add_source_database(parent))
  return metadata_rows

def add_fragments(parent):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the fragments of the protein."""

  if not parent['Assigned Protein Fragments']: return []

  with open(Path(__file__).parents[1] / 'data' / 'fragment-type.json', 'r') as f:
    fragment_type_map = json.load(f)

  fragments = parent['Assigned Protein Fragments'].split(', ')
  fragment_count = len(fragments)

  if fragment_count < 2: return []

  fragment_count = 0
  fragment_rows = []

  for fragment in fragments:

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
      f"UP:{parent['Parent Protein ID']}-fragment-{fragment_count+1}",
      f"{fragment_type} ({fragment_start}-{fragment_end})",
      f"UP:{parent['Parent Protein ID']}"
    ))

    fragment_rows.extend(
      [triple(
        f"UP:{parent['Parent Protein ID']}-fragment-{fragment_count+1}",
        "ONTIE:0003627",
        fragment_start,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{parent['Parent Protein ID']}-fragment-{fragment_count+1}",
        "ONTIE:0003628",
        fragment_end,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{parent['Parent Protein ID']}-fragment-{fragment_count+1}",
        "ONTIE:0003620",
        f"{fragment_count+1} {fragment_type} ({fragment_start}-{fragment_end})",
        datatype="xsd:string"
      )]
    )

    fragment_count += 1

  return fragment_rows

def add_reviewed_status(parent):
  """Given a row from the source_assignments dataframe, return a triple
  of the reviewed status of the protein."""
  return [triple(
    f"UP:{parent['Parent Protein ID']}",
    "UC:reviewed",
    "true" if parent['Assigned Protein Review Status'] else "false",
    datatype="xsd:boolean"
  )]

def add_synonyms(parent):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the synonyms of the protein."""
  synonyms = []
  for synonym in parent['Assigned Protein Synonyms'].split(', '):
    if '@' in synonym or '{' in synonym: continue
    synonyms.append(triple(
      f"UP:{parent['Parent Protein ID']}",
      "ONTIE:0003622",
      synonym,
      datatype="xsd:string"
    ))
  return synonyms

def add_accession(parent):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the accession of the protein and the URL to the UniProt entry."""
  return [triple(
    f"UP:{parent['Parent Protein ID']}",
    "ONTIE:0003623",
    parent['Parent Protein ID'],
    datatype="xsd:string"
  ),
  triple(
    f"UP:{parent['Parent Protein ID']}",
    "ONTIE:0003624",
    f"https://www.uniprot.org/uniprotkb/{parent['Parent Protein ID']}",
    datatype="xsd:string"
  )]

def add_source_database(parent):
  """Given a row from the source_assignments dataframe, return a triple
  for the source database of the protein (always UniProt)."""
  return [triple(
    f"UP:{parent['Parent Protein ID']}",
    "ONTIE:0003625",
    "UniProt",
    datatype="xsd:string"
  )]

def owl_class(subject, label, parent):
  """Given a subject, label, and parent,
  return triples defining an owl:Class in the protein tree."""
  return [
    triple(subject, 'rdf:type', 'owl:Class'),
    triple(subject, 'rdfs:label', label, 'xsd:string'),
    triple(subject, 'rdfs:subClassOf', parent)
  ]

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

if __name__ == "__main__":
  build_path = Path(__file__).parents[3] / 'build'
  
  assignments = pl.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')
  source_data = pl.read_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')
  assignments = assignments.join(source_data, how='left', on='Source Accession', coalesce=True)

  with sqlite3.connect(build_path / 'arborist' / 'nanobot.db') as connection:
    # Copy the organism_tree, but replace each taxon with 'taxon protein'.
    tree_base = pl.read_database('''
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
  lower_subjects = tree_base.filter( pl.col('object').is_in(['lower', '']))['subject'].to_list()
  tree_base = tree_base.filter(~pl.col('subject').is_in(lower_subjects))

  # Relabel 'taxon' to 'taxon protein'.
  tree_base = tree_base.with_columns(
    pl.when(
      pl.col('subject').str.starts_with('iedb-protein:'),
      pl.col('predicate') == 'rdfs:label'
    )
    .then(pl.col('object') + pl.lit(' protein'))
    .otherwise(pl.col('object'))
    .alias('object')
  )

  # Add top-level 'protein'
  new_rows = owl_class('PR:000000001', 'protein', 'BFO:0000040')
  tree_base = pl.concat([tree_base, pl.DataFrame(new_rows)])

  # Re-parent children of 'Root' and 'organism' to 'protein'
  tree_base = tree_base.with_columns(
    pl.when(
      (pl.col('object') == 'iedb-protein:1') | 
      (pl.col('object') == 'iedb-protein:28384') | 
      (pl.col('object') == 'OBI:0100026')
    )
    .then(pl.lit('PR:000000001'))
    .otherwise(pl.col('object'))
    .alias('object')
  )

  protein_tree = build_tree(tree_base, assignments)

  db = 'sqlite:///' + str(build_path / 'arborist' / 'nanobot.db')
  protein_tree.write_database('protein_tree_old', connection=db, if_table_exists='replace')
  protein_tree.write_database('protein_tree_new', connection=db, if_table_exists='replace')