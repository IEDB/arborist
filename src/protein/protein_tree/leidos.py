import polars as pl
from pathlib import Path


def get_parents(assignments):
  top_parents = assignments.group_by('Source Accession').agg(
    pl.col('Assigned Protein ID').mode().first()
  ).rename({'Assigned Protein ID': 'Parent Protein ID'})

  assignments = assignments.join(
    top_parents, how='left', on='Source Accession', coalesce=True
  )

  parents = assignments.filter(
    (pl.col('Assigned Protein ID') == pl.col('Parent Protein ID')) | 
    (pl.col('Assigned Protein ID').is_null())
  ).unique(subset=['Source Accession', 'Parent Protein ID'])

  return parents

def make_source_parents(parents):
  all_data = parents.join(
    source_data, how='left', on='Source Accession', coalesce=True
  )
  all_data = all_data.join(
    species_data, how='left', left_on='Species Taxon ID', right_on='Species ID', coalesce=True
  )

def make_parent_proteins():
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

  parents = get_parents(assignments)
  make_source_parents(parents)
  make_parent_proteins()
  epitope_mapping = EpitopeMapping(assignments)