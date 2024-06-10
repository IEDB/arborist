import json
import argparse
import subprocess
import polars as pl
from pathlib import Path
from pepmatch import Preprocessor, Matcher
from ARC.classifier import SeqClassifier
from protein_tree.data_fetch import DataFetcher


class AssignmentHandler:
  def __init__(self, taxon_id, species_name, group, peptides, sources, num_threads):
    self.taxon_id = taxon_id
    self.species_name = species_name
    self.group = group
    self.peptides = peptides
    self.sources = sources
    self.num_threads = num_threads
    self.build_path = Path(__file__).parents[3] / 'build'
    self.species_path = self.build_path / 'species' / str(self.taxon_id)

  def process_species(self):
    source_processor = SourceProcessor(
      self.taxon_id, self.group, self.sources, self.num_threads, self.species_path
    )
    source_assignments = source_processor.process()
    self.peptides = self.peptides.join(source_assignments, how='left', on='Source Accession', coalesce=True)

    peptide_processor = PeptideProcessor(
      self.taxon_id, self.species_name, self.peptides, source_assignments, self.species_path
    )
    peptide_processor.process()

  def cleanup_files(self):
    pass


class SourceProcessor:
  def __init__(self, taxon_id, group, sources, num_threads, species_path):
    self.taxon_id = taxon_id
    self.group = group
    self.sources = sources
    self.species_path = species_path
    self.num_threads = num_threads
    self.bin_path = Path(__file__).parents[3] / 'bin'
    self.fasta_path = self.species_path / 'sources.fasta'
    self.proteome_path = self.species_path / 'proteome.fasta'
    self.run_mmseqs2 = self.taxon_id in [9606, 10090, 10116] # human, mouse, rat

  def process(self):
    self.write_source_data()
    self.write_sources_to_fasta(self.sources)

    if self.run_mmseqs2:
      self.mmseqs2_sources()
    else:
      self.blast_sources()

    top_proteins = self.select_top_proteins()
    top_protein_data = self.get_protein_data(top_proteins)

    if self.group == 'vertebrate':
      arc_df = self.run_arc()
      source_assignments = self.combine_arc_data(top_protein_data, arc_df)
    else:
      source_assignments = top_protein_data.with_columns(
        pl.lit(None).alias('ARC Assignment')
      )
    return source_assignments

  def write_source_data(self):
    self.sources.write_csv(self.species_path / 'source-data.tsv', separator='\t')

  def write_sources_to_fasta(self, df):
    with open(self.fasta_path, 'w') as fasta_file:
      for row in df.iter_rows(named=True):
        fasta_file.write(f">{row['Source Accession']}\n{row['Sequence']}\n")

  def blast_sources(self):
    makeblastdb_path = self.bin_path / 'makeblastdb'
    blastp_path = self.bin_path / 'blastp'

    makeblastdb_cmd = [
      str(makeblastdb_path), 
      '-in', str(self.proteome_path), 
      '-dbtype', 'prot', 
      '-out', str(self.proteome_path)
    ]
    blastp_cmd = [
      str(blastp_path), 
      '-query', str(self.fasta_path), 
      '-db', str(self.proteome_path), 
      '-evalue', '1', 
      '-outfmt', '10', 
      '-num_threads', str(self.num_threads), 
      '-out', str(self.species_path / 'alignments.csv')
    ]
    subprocess.run(makeblastdb_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(blastp_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  def mmseqs2_sources(self):
    mmseqs2_path = self.bin_path / 'mmseqs'

    mmseqs2_cmd = [
      str(mmseqs2_path), 'easy-search',
      str(self.fasta_path), str(self.proteome_path),
      str(self.species_path / 'alignments.tsv'),
      str(self.species_path / 'tmp'),
      '--threads', str(self.num_threads),
      '-s', '7.0'
    ]
    subprocess.run(mmseqs2_cmd, check=True)

  def run_arc(self):
    arc_results_exist = (self.species_path / 'arc-results.tsv').exists()
    if arc_results_exist:
      old_arc_df = pl.read_csv(self.species_path / 'arc-results.tsv', separator='\t')
      new_arc_sources = self.sources.filter(
        ~pl.col('Source Accession').str.contains('SRC'),
        ~pl.col('Source Accession').is_in(old_arc_df['id'].to_list())
      )
      if new_arc_sources.shape[0] == 0:
        return old_arc_df
      self.write_sources_to_fasta(new_arc_sources)

    temp_results_path = self.species_path / 'arc-temp-results.tsv'
    SeqClassifier(
      outfile=str(temp_results_path),
      threads=self.num_threads,
      blast_path=str(self.bin_path)+'/'
    ).classify_seqfile(f'{self.species_path}/sources.fasta')

    temp_results = pl.read_csv(temp_results_path, separator='\t')
    if arc_results_exist:
      arc_df = pl.concat([old_arc_df, temp_results])
    else:
      arc_df = temp_results

    arc_df.write_csv(self.species_path / 'arc-results.tsv', separator='\t')

    return arc_df
  
  def combine_arc_data(self, top_protein_data, arc_df):
    arc_df = arc_df.with_columns(
      pl.concat_str(
        [pl.col('class'), pl.col('chain_type'), pl.col('calc_mhc_allele')], separator='_', ignore_nulls=True
      ).alias('ARC Assignment')
    )
    arc_df = arc_df.select('id', 'ARC Assignment')
    top_protein_data = top_protein_data.join(
      arc_df, how='left', left_on='Source Accession', right_on='id', coalesce=True
    )
    return top_protein_data

  def select_top_proteins(self):
    alignment_cols = [
      'Query', 'Subject', '% Identity', 'Alignment Length', 'Mismatches', 'Gap Openings', 
      'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-value', 'Bit Score'
    ]
    if self.run_mmseqs2:
      alignments = pl.read_csv(self.species_path / 'alignments.tsv', separator='\t', has_header=False, new_columns=alignment_cols)
      alignments = alignments.with_columns(pl.col('% Identity').mul(100).alias('% Identity'))
    else:
      alignments = pl.read_csv(self.species_path / 'alignments.csv', separator=',', has_header=False, new_columns=alignment_cols)
      alignments = alignments.with_columns(pl.col('Subject').str.split('|').list.get(1))

    query_length_map = dict(self.sources.select('Source Accession', 'Length').iter_rows())
    alignments = alignments.with_columns(
      pl.col('Query').replace(query_length_map).cast(pl.Int32).alias('Query Length')
    )
    alignments = alignments.with_columns(
      pl.col('Alignment Length').mul(pl.col('% Identity')).truediv(pl.col('Query Length')).alias('Score')
    )
    return alignments.group_by('Query').agg(pl.all().sort_by('Score').last())

  def get_protein_data(self, top_proteins):
    proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t')
    protein_data = top_proteins.join(
      proteome, how='left', left_on='Subject', right_on='Protein ID', coalesce=False
    )
    protein_data = protein_data.select(pl.col(
      'Query', 'Score', 'Gene', 'Protein ID', 'Protein Name'
    )).rename({
      'Query': 'Source Accession', 'Score': 'Source Alignment Score', 'Gene': 'Source Assigned Gene',
      'Protein ID': 'Source Assigned Protein ID', 'Protein Name': 'Source Assigned Protein Name'}
    )
    return protein_data


class PeptideProcessor:
  def __init__(self, taxon_id, species_name, peptides, source_assignments, species_path):
    self.taxon_id = taxon_id
    self.species_name = species_name
    self.peptides = peptides
    self.source_assignments = source_assignments
    self.species_path = species_path

  def process(self):
    self.preprocess_proteome()
    self.search_peptides()
    assignments = self.assign_parents()
    assignments = self.get_protein_data(assignments)
    assignments = self.handle_allergens(assignments)
    assignments = self.handle_manuals(assignments)
    self.write_assignments(assignments)

  def preprocess_proteome(self):
    if (self.species_path / 'proteome.db').exists():
      return
    gp_proteome = self.species_path / 'gp_proteome.fasta' if (self.species_path / 'gp_proteome.fasta').exists() else ''
    Preprocessor(
      proteome = self.species_path / 'proteome.fasta',
      preprocessed_files_path = self.species_path,
      gene_priority_proteome=gp_proteome
    ).sql_proteome(k = 5)

  def search_peptides(self):
    peptides = self.peptides['Sequence'].to_list()
    Matcher(
      query=peptides,
      proteome_file=self.species_path / 'proteome.fasta',
      max_mismatches=0,
      k=5,
      preprocessed_files_path=self.species_path,
      best_match=False, 
      output_format='dataframe',
      output_name=self.species_path / 'peptide-matches.tsv',
      sequence_version=False
    ).match()
  
  def assign_parents(self):
    matches = pl.read_csv(self.species_path / 'peptide-matches.tsv', separator='\t')
    peptides_with_genes = self.peptides.filter(pl.col('Source Assigned Gene').is_not_null())
    peptides_without_genes = self.peptides.filter(pl.col('Source Assigned Gene').is_null())

    matches_with_genes = peptides_with_genes.join(
      matches, how="left", coalesce=False,
      left_on=["Sequence", "Source Assigned Gene"], 
      right_on=["Query Sequence", "Gene"]
    )
    matches_without_genes = peptides_without_genes.join(
      matches, how="left", coalesce=False,
      left_on=["Sequence", "Source Assigned Protein ID"], 
      right_on=["Query Sequence", "Protein ID"]
    )

    top_matches_with_genes = matches_with_genes.sort(
      ["Sequence", "SwissProt Reviewed", "Gene Priority", "Protein Existence Level"],
      descending=[False, True, True, False]
    ).group_by("Sequence").first()

    top_matches_without_genes = matches_without_genes.sort(
      ["Sequence", "SwissProt Reviewed", "Gene Priority", "Protein Existence Level"],
      descending=[False, True, True, False]
    ).group_by("Sequence").first()

    match_cols = ['Sequence', 'Gene', 'Protein ID', 'Protein Name', 'Index start', 'Index end', 'SwissProt Reviewed']
    top_matches_with_genes = top_matches_with_genes.select(match_cols)
    top_matches_without_genes = top_matches_without_genes.select(match_cols)

    assignments_with_genes = peptides_with_genes.join(
      top_matches_with_genes, how="left", coalesce=False,
      left_on=["Sequence", "Source Assigned Gene"],
      right_on=["Sequence", "Gene"],
    )

    assignments_without_genes = peptides_without_genes.join(
      top_matches_without_genes, how="left", coalesce=False,
      left_on=["Sequence", "Source Assigned Protein ID"], 
      right_on=["Sequence", "Protein ID"]
    )
  
    assignments = pl.concat([assignments_with_genes, assignments_without_genes])
    assignments = assignments.drop('Sequence_right', 'Gene')
    assignments = assignments.unique(subset=['Sequence', 'Source Accession'])
    assignments.write_csv(self.species_path / 'peptide-assignments.tsv', separator='\t')
    return assignments

  def get_protein_data(self, assignments):
    proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t')
    proteome = proteome.select(pl.col('Protein ID', 'Sequence'))
    proteome = proteome.rename({'Sequence': 'Assigned Protein Sequence'})
    fragments = self.get_fragment_data()
    assignments = assignments.join(
      proteome, how='left', on='Protein ID', coalesce=True
    )
    assignments = assignments.with_columns(
      pl.col('Assigned Protein Sequence').str.len_chars().alias('Assigned Protein Length'),
      pl.col('Protein ID').replace(fragments, default="").alias('Assigned Protein Fragments'),
      pl.lit(str(self.taxon_id)).alias('Species Taxon ID'),
      pl.lit(self.species_name).alias('Species Name')
    )

    return assignments
  
  def get_fragment_data(self):
    if (self.species_path / 'fragment-data.json').exists():
      with open(self.species_path / 'fragment-data.json', 'r') as f:
        return json.load(f)
      
  def handle_allergens(self, assignments):
    pass

  def handle_manuals(self, assignments):
    pass

  def write_assignments(self, assignments):
    assignments = assignments.rename({
      'Sequence': 'Epitope Sequence',
      'Starting Position': 'Source Starting Position',
      'Ending Position': 'Source Ending Position',
      'Protein ID': 'Assigned Protein ID',
      'Protein Name': 'Assigned Protein Name',
      'Index start': 'Assigned Protein Starting Position',
      'Index end': 'Assigned Protein Ending Position',
      'SwissProt Reviewed': 'Assigned Protein Review Status'
    })
    col_order = [
      'Species Taxon ID', 'Species Name', 'Organism ID', 'Source Accession', 
      'Source Alignment Score', 'Source Assigned Gene', 'Source Assigned Protein ID', 
      'Source Assigned Protein Name', 'ARC Assignment', 'Epitope ID', 'Epitope Sequence', 
      'Source Starting Position', 'Source Ending Position', 'Assigned Protein ID', 
      'Assigned Protein Name', 'Assigned Protein Review Status', 
      'Assigned Protein Starting Position', 'Assigned Protein Ending Position', 
      'Assigned Protein Sequence', 'Assigned Protein Length', 'Assigned Protein Fragments'
    ]
    assignments = assignments.select(col_order)
    assignments.write_csv(self.species_path / 'peptide-assignments.tsv', separator='\t')


def do_assignments(taxon_id):
  species_row = active_species.row(by_predicate=pl.col('Species ID') == taxon_id)
  species_name = species_row[2]
  group = species_row[4]
  active_taxa = [int(taxon_id) for taxon_id in species_row[3].split(', ')]
  peptides = data_fetcher.get_peptides_for_species(all_peptides, active_taxa)
  sources = data_fetcher.get_sources_for_species(all_sources, peptides['Source Accession'].to_list())
  sources = sources.with_columns(pl.col('Sequence').str.len_chars().alias('Length'))
  config = {
    'taxon_id': taxon_id,
    'species_name': species_name,
    'group': group,
    'peptides': peptides,
    'sources': sources,
    'num_threads': args.num_threads
  }
  assignment_handler = AssignmentHandler(**config)
  assignment_handler.process_species()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-n', '--num_threads', type=int, default=1, help='Number of threads to use.')
  parser.add_argument('-t', '--taxon_id', type=int, help='Taxon ID of the species to process.')
  args = parser.parse_args()

  taxon_id = args.taxon_id
  all_species = not bool(taxon_id)

  build_path = Path(__file__).parents[3] / 'build'
  active_species = pl.read_csv(build_path / 'arborist' / 'active-species.tsv', separator='\t')
  data_fetcher = DataFetcher(build_path)
  all_peptides = data_fetcher.get_all_peptides()
  all_sources = data_fetcher.get_all_sources()

  if all_species:
    for row in active_species.iter_rows(named=True):
      do_assignments(row['Species ID'])
  else:
    do_assignments(taxon_id)