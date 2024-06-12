import polars as pl
from pathlib import Path


def get_all_parent_data(assignments, source_data, species_data):
  top_parents = assignments.group_by('Source Accession').agg(
    pl.col('Assigned Protein ID').mode().first()
  ).rename({'Assigned Protein ID': 'Parent Protein ID'})

  assignments = assignments.join(
    top_parents, how='left', on='Source Accession', coalesce=True
  ).filter(
    (pl.col('Assigned Protein ID') == pl.col('Parent Protein ID')) | 
    (pl.col('Assigned Protein ID').is_null())
  ).unique(subset=['Source Accession', 'Parent Protein ID'])

  all_parent_data = assignments.join(
    source_data, how='left', on='Source Accession', coalesce=True
  ).join(
    species_data, how='left', left_on='Species Taxon ID', right_on='Species ID', coalesce=True
  )
  return all_parent_data

def make_source_parents(all_parent_data):
  all_parent_data = all_parent_data.with_columns(
    pl.when(pl.col('Source Alignment Score') >= 90)
      .then(pl.lit('strong-blast-match'))
      .when(pl.col('Source Alignment Score') < 90)
      .then(pl.lit('weak-blast-match'))
      .otherwise(pl.lit('manual')).alias('Protein Strategy'),
    pl.when(pl.col('Parent Protein ID').is_not_null())
      .then(pl.lit('https://www.uniprot.org/uniprotkb/') + pl.col('Parent Protein ID'))
      .otherwise(
        pl.lit('https://ontology.iedb.org/taxon-protein/') + 
        pl.col('Species Taxon ID').cast(pl.String) + pl.lit('-other')
      )
      .alias('Parent IRI'),
    pl.when(pl.col('Parent Protein ID').is_not_null())
      .then(pl.lit('UniProt'))
      .otherwise(pl.lit('IEDB'))
      .alias('Parent Protein Database')
  )
  source_parents = all_parent_data.select(
    'Source ID', 'Source Accession', 'Database', 'Name', 'Aliases', 'Assigned Protein Synonyms',
    'Organism ID', 'Organism Name', 'Species Taxon ID', 'Species Name', 'Proteome ID',
    'Proteome Label', 'Protein Strategy', 'Parent IRI', 'Parent Protein Database',
    'Parent Protein ID', 'Assigned Protein Length', 'Assigned Protein Sequence'
  ).rename({
    'Source Accession': 'Accession', 'Assigned Protein Synonyms': 'Synonyms',
    'Organism ID': 'Taxon ID', 'Organism Name': 'Taxon Name', 'Species Taxon ID': 'Species ID',
    'Species Name': 'Species Label', 'Parent Protein ID': 'Parent Protein Accession', 
    'Assigned Protein Length': 'Parent Sequence Length', 'Assigned Protein Sequence': 'Sequence'
  })
  source_parents.write_csv(build_path / 'arborist' / 'source-parents.tsv', separator='\t')

def make_parent_proteins(all_parent_data):
  pass


class EpitopeMapping:
  def __init__(self, assignments):
    self.assignments = assignments

  def make_epitope_mapping(self):
    exact_match_df, non_exact_match_df = self.split_peptide_assignments()
    self.blast_non_exact_matches(non_exact_match_df)

  def split_peptide_assignments(self):
    exact_match_df = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_not_null()) & 
      (pl.col('Assigned Protein ID').is_not_null())
    )
    non_exact_match_df = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_null()) & 
      (pl.col('Assigned Protein ID').is_not_null())
    )
    return exact_match_df, non_exact_match_df
  
  def blast_non_exact_matches(self, non_exact_match_df):
    pass

if __name__ == '__main__':
  build_path = Path(__file__).parents[3] / 'build'
  assignments = pl.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')
  source_data = pl.read_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')
  species_data = pl.read_csv(build_path / 'arborist' / 'all-species-data.tsv', separator='\t')

  all_parent_data = get_all_parent_data(assignments, source_data, species_data)
  make_source_parents(all_parent_data)
  make_parent_proteins(all_parent_data)
  epitope_mapping = EpitopeMapping(assignments)