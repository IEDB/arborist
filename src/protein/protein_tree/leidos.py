import argparse
import subprocess
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
  parent_protein_name_tail = (
    pl.lit('|') + pl.col('Parent Protein ID') + pl.lit('|') + pl.col('Assigned Protein Entry Name')
  )
  unique_parents = all_parent_data.unique(subset=['Parent Protein ID']).with_columns(
    pl.lit('UniProt').alias('Parent Protein Database'),
    pl.when(pl.col('Assigned Protein Review Status'))
    .then(pl.lit('sp') + parent_protein_name_tail)
    .otherwise(pl.lit('tr') + parent_protein_name_tail)
    .alias('Parent Protein Name')
    ).filter(
      pl.col('Parent Protein ID').is_not_null()
    )

  parent_proteins = unique_parents.select(
    'Parent Protein ID', 'Parent Protein Database', 'Parent Protein Name', 'Assigned Protein Name',
    'Proteome ID', 'Proteome Label', 'Assigned Protein Sequence'
  ).rename({
    'Parent Protein ID': 'Accession', 'Parent Protein Database': 'Database',
    'Assigned Protein Name': 'Title', 'Assigned Protein Sequence': 'Sequence'
  })
  parent_proteins.write_csv(build_path / 'arborist' / 'parent-proteins.tsv', separator='\t')


class EpitopeMapping:
  def __init__(self, assignments, num_threads):
    self.assignments = assignments
    self.num_threads = num_threads
    self.bin_path = Path(__file__).parents[3] / 'bin'

  def make_epitope_mapping(self):
    exact_match_df, non_exact_match_df, discontinuous_df = self.split_peptide_assignments()
    exact_mappings = self.calculate_exact_mappings(exact_match_df)
    non_exact_mappings = self.blast_non_exact_matches(non_exact_match_df)
    
  def split_peptide_assignments(self):
    exact_match_df = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_not_null()) & 
      (pl.col('Assigned Protein ID').is_not_null())
    )
    non_exact_match_df = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_null()) & 
      (pl.col('Assigned Protein ID').is_not_null()) &
      (~pl.col('Epitope Sequence').str.contains('0|1|2|3|4|5|6|7|8|9')) &
      (pl.col('Epitope ID').is_not_null())
    )
    discontinuous_df = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_null()) & 
      (pl.col('Assigned Protein ID').is_not_null()) &
      (pl.col('Epitope Sequence').str.contains('0|1|2|3|4|5|6|7|8|9')) &
      (pl.col('Epitope ID').is_not_null())
    )
    return exact_match_df, non_exact_match_df, discontinuous_df
  
  def calculate_exact_mappings(self, exact_match_df):
    pass

  def blast_non_exact_matches(self, non_exact_match_df):
    blast_temp = build_path / 'blast_temp'
    blast_temp.mkdir(exist_ok=True)
    self.make_blast_db(non_exact_match_df, blast_temp)
    self.make_epitope_file(non_exact_match_df, blast_temp)
    cmd = [
      str(self.bin_path / 'blastp'),
      '-query', str(blast_temp / 'epitopes.fasta'),
      '-db', str(blast_temp / 'proteins.fasta'),
      '-outfmt', '10',
      '-num_threads', str(self.num_threads),
      '-out', str(blast_temp / 'blast_results.csv')
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return self.calculate_blast_mappings(blast_temp)

  def make_blast_db(self, non_exact_match_df, blast_temp):
    with open(str(blast_temp / 'proteins.fasta'), 'w') as f:
      for row in non_exact_match_df.iter_rows(named=True):
        f.write(f'>{row["Assigned Protein ID"]}\n{row["Assigned Protein Sequence"]}\n')
    cmd = [
      str(self.bin_path / 'makeblastdb'),
      '-in', str(blast_temp / 'proteins.fasta'),
      '-dbtype', 'prot',
      '-out', str(blast_temp / 'proteins.fasta')
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  def make_epitope_file(self, non_exact_match_df, blast_temp):
    with open(str(blast_temp / 'epitopes.fasta'), 'w') as f:
      for row in non_exact_match_df.iter_rows(named=True):
        f.write(f'>{int(row["Epitope ID"])}\n{row["Epitope Sequence"]}\n')

  def calculate_blast_mappings(self, blast_temp):
    blast_cols = [
      'Query', 'Subject', '% Identity', 'Alignment Length', 'Mismatches', 'Gap Openings', 
      'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-value', 'Bit Score'
    ]
    blast_results = pl.read_csv(
      blast_temp / 'blast_results.csv', separator=',', has_header=False, new_columns=blast_cols
    )


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-n', '--num_threads', type=int, default=1, help='Number of threads to use.')
  args = parser.parse_args()

  build_path = Path(__file__).parents[3] / 'build'
  assignments = pl.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')
  source_data = pl.read_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')
  species_data = pl.read_csv(build_path / 'arborist' / 'all-species-data.tsv', separator='\t')

  # all_parent_data = get_all_parent_data(assignments, source_data, species_data)
  # make_source_parents(all_parent_data)
  # make_parent_proteins(all_parent_data)

  epitope_mapping = EpitopeMapping(assignments, args.num_threads)
  epitope_mapping.make_epitope_mapping()