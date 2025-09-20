import ast
import sqlite3
import polars as pl
from pathlib import Path

def build_tree(tree_base, assignments, gene_layer=False):
  """Given the tree_base from the organism tree dataframes for the tree and
  the assignments, insert LDTab rows for each protein under its
  'species protein' parent."""
  new_rows = []

  arc_parents = assignments.filter(
    (pl.col('ARC Assignment').str.contains('TCR')) |
    (pl.col('ARC Assignment').str.contains('BCR')) |
    (pl.col('ARC Assignment').str.contains('MHC'))
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

  new_rows.extend(add_normal_parents(normal_parents, gene_layer=gene_layer))
  new_rows.extend(add_arc_parents(arc_parents))
  new_rows.extend(add_others(others))

  protein_tree = pl.concat([tree_base, pl.DataFrame(new_rows)])

  return protein_tree

def add_normal_parents(normal_parents, gene_layer):
  """Proteins that have been assigned a parent."""
  rows = []
  unknown_gene_nodes_seen = set()
  gene_node_map = {}
  case_variant_counts = {}

  for parent in normal_parents.iter_rows(named=True):
    if gene_layer:
      gene = parent['Source Assigned Gene'] if parent['Source Assigned Gene'] else 'Unknown'
      gene = gene.replace('\\', '-')
      species_id = parent['Species Taxon ID']

      if gene == 'Unknown':
        gene_node_id = f"iedb-protein:{species_id}-Unknown"
        if species_id not in unknown_gene_nodes_seen:
          rows.extend(
            owl_class(
              gene_node_id,
              "Gene: Unknown",
              f"iedb-protein:{species_id}"
            )
          )
          rows.extend(
            [triple(
              gene_node_id,
              "ONTIE:0003674",
              "Unknown",
              "xsd:string"
            )]
          )
          unknown_gene_nodes_seen.add(species_id)
      else:  # Known gene
        if (species_id, gene) in gene_node_map:
          gene_node_id = gene_node_map[(species_id, gene)]
        else:
          # Create a new node for this unique gene symbol
          tracking_key = (species_id, gene.lower())
          count = case_variant_counts.get(tracking_key, 0)
          gene_for_id = f"{gene}_dup{count}" if count > 0 else gene
          gene_node_id = f"iedb-protein:{species_id}-{gene_for_id}"

          gene_node_map[(species_id, gene)] = gene_node_id
          case_variant_counts[tracking_key] = count + 1

          rows.extend(
            owl_class(
              gene_node_id,
              f"Gene: {gene}",
              f"iedb-protein:{species_id}"
            )
          )
          rows.extend(
            [triple(
              gene_node_id,
              "ONTIE:0003674",
              f"{gene}",
              "xsd:string"
            )]
          )
      
      # Attach the protein to its gene node
      rows.extend(
        owl_class(
          f"UP:{parent['Parent Protein ID']}",
          f"{parent['Assigned Protein Name']} (UniProt:{parent['Parent Protein ID']})",
          gene_node_id
        )
      )
    else:
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
  receptor_name_dict = {
    'TCR': 'T Cell Receptor', 'BCR': 'B Cell Receptor / Immunoglobulin',
    'MHC-I': 'Major Histocompatibility Complex I', 'MHC-II': 'Major Histocompatibility Complex II'
  }
  rows = []
  nodes_seen = set()
  genes_seen_in_receptor = set()  # group genes together
  arc_parents = arc_parents.rename({'Assigned Protein ID': 'Parent Protein ID'})
  for parent in arc_parents.iter_rows(named=True):
    arc_assignment = parent['ARC Assignment']
    antigen_receptor_class = arc_assignment.split('_')[0]
    node_name = receptor_name_dict[antigen_receptor_class]
    node_id = f"{str(parent['Species Taxon ID'])}-{antigen_receptor_class.lower()}"
    if antigen_receptor_class == 'BCR':  # add to end since there are "bcr" genes in some species
      node_id += '-ig'
    if node_id not in nodes_seen:
      rows.extend(
        owl_class(
          f"iedb-protein:{node_id}",
          f"{node_name} chain",
          f"iedb-protein:{parent['Species Taxon ID']}"
        )
      )
      nodes_seen.add(node_id)
    
    # group genes within receptor nodes
    gene = parent['Source Assigned Gene'] if parent['Source Assigned Gene'] else 'Unknown'
    gene = gene.replace('\\', '-').upper()
    gene_node_id = f"{node_id}-{gene}"
    if (node_id, gene) not in genes_seen_in_receptor:
      rows.extend(
        owl_class(
          f"iedb-protein:{gene_node_id}",
          f"Gene: {gene}",
          f"iedb-protein:{node_id}"
        )
      )
      genes_seen_in_receptor.add((node_id, gene))

    rows.extend(
      owl_class(
        f"UP:{parent['Parent Protein ID']}",
        f"{parent['Assigned Protein Name']} (UniProt:{parent['Parent Protein ID']})",
        f"iedb-protein:{gene_node_id}"
      )
    )
    rows.extend(add_metadata(parent))
  return rows

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

    prefix = 'UP' if other['Database'] == 'UniProt' else 'NCBI'
    rows.extend(
      owl_class(
        f"{prefix}:{other['Source Accession']}",
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
  metadata_rows.extend(add_canonical_status(parent))
  metadata_rows.extend(add_synonyms(parent))
  metadata_rows.extend(add_accession(parent))
  metadata_rows.extend(add_source_database(parent))
  return metadata_rows

def add_fragments(parent):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the fragments of the protein."""
  fragment_str = parent['Assigned Protein Fragments']
  if not fragment_str: return []

  fragments = ast.literal_eval(fragment_str)
  if len(fragments) < 2: return []  # don't fragment the protein if there is only one

  fragment_type_map = {
    "Chain": "mature protein", "Peptide": "peptide", "Propeptide": "propeptide",
    "Signal": "signal peptide", "Transit peptide": "transit peptide"
  }

  fragment_counter = 1
  fragment_rows = []

  for fragment in fragments:
    fragment_type  = fragment['type']
    fragment_type  = fragment_type_map[fragment_type]
    fragment_start = fragment['start']
    fragment_end   = fragment['end']
    fragment_desc  = fragment['description']
    fragment_id    = fragment['feature_id']

    try:
      int(fragment_start)
      int(fragment_end)
    except TypeError:  # these are nulls
      continue

    if fragment_start == 1 and fragment_end == parent['Assigned Protein Length']:
      continue  # do not fragment whole length chains

    if fragment_type == 'mature protein' and fragment_desc and fragment_desc != parent['Assigned Protein Name']:
      fragment_type = f'{fragment_type} ({fragment_desc})'

    if fragment_id == 'N/A':  # use counter for fragment id
      fragment_id = f'fragment-{fragment_counter}'
      fragment_counter += 1

    fragment_rows.extend(owl_class(
      f"UP:{parent['Parent Protein ID']}-{fragment_id}",
      f"{fragment_type} ({fragment_start}-{fragment_end})",
      f"UP:{parent['Parent Protein ID']}"
    ))

    fragment_rows.extend(
      [triple(
        f"UP:{parent['Parent Protein ID']}-{fragment_id}",
        "ONTIE:0003627",
        fragment_start,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{parent['Parent Protein ID']}-{fragment_id}",
        "ONTIE:0003628",
        fragment_end,
        datatype="xsd:integer"
      ),
      triple(
        f"UP:{parent['Parent Protein ID']}-{fragment_id}",
        "ONTIE:0003620",
        f"{fragment_id} {fragment_type} ({fragment_start}-{fragment_end})",
        datatype="xsd:string"
      )]
    )
  return fragment_rows

def add_canonical_status(parent):
  """Given a row from the source_assignments dataframe, return a triple
  of the canonical status protein."""
  return [triple(
    f"UP:{parent['Parent Protein ID']}",
    "ONTIE:0003673",
    "true" if parent['Assigned Protein Review Status'] and '-' not in parent['Parent Protein ID'] else "false",
    datatype="xsd:boolean"
  )]

def add_synonyms(parent):
  """Given a row from the source_assignments dataframe, return a list of triples
  for the synonyms of the protein."""
  synonyms = []
  protein_synonyms = parent['Assigned Protein Synonyms']
  if not protein_synonyms: return []
  
  for synonym in protein_synonyms.split(', '):
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
    f"http://www.uniprot.org/uniprot/{parent['Parent Protein ID']}",
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
  numeric_cols = [
    'Epitope ID',
    'Species Taxon ID',
    'Organism ID',
    'Source Starting Position',
    'Source Ending Position',
    'Assigned Protein Starting Position',
    'Assigned Protein Ending Position'
  ]
  assignments = pl.read_csv(
    build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t', schema_overrides={col: pl.Float64 for col in numeric_cols}
  ).with_columns(pl.col(numeric_cols).cast(pl.Int64))

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
  lower_subjects = tree_base.filter(pl.col('object').is_in(['lower', '']))['subject'].to_list()
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
  protein_tree_w_gene = build_tree(tree_base, assignments, gene_layer=True)

  db = 'sqlite:///' + str(build_path / 'arborist' / 'nanobot.db')
  protein_tree.write_database('protein_tree_old', connection=db, if_table_exists='replace')
  protein_tree_w_gene.write_database('protein_tree_new', connection=db, if_table_exists='replace')
