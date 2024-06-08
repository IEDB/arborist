import argparse
import subprocess
import polars as pl
from pathlib import Path

from protein_tree.data_fetch import DataFetcher


class AssignmentHandler:
  def __init__(self, taxon_id, group, peptides, sources, num_threads):
    self.taxon_id = taxon_id
    self.group = group
    self.peptides = peptides
    self.sources = sources
    self.num_threads = num_threads
    self.build_path = Path(__file__).parents[3] / 'build' 
    self.species_path = self.build_path / 'species' / str(self.taxon_id)

  def process_species(self):
    source_processor = SourceProcessor(self.taxon_id, self.group, self.sources, self.num_threads, self.species_path)
    source_processor.process()

    peptide_processor = PeptideProcessor()
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
    self.write_isoforms_to_fasta()

    if self.group == 'vertebrate':
      self.run_arc()
      # arc_results = pd.read_csv(self.species_path / 'arc-results.tsv', sep='\t')

    if self.run_mmseqs2:
      self.mmseqs2_isoforms()
    else:
      self.blast_isoforms()

    top_proteins = self.select_top_protein()
    top_protein_data = self.get_protein_data(top_proteins)

  def write_source_data(self):
    self.sources.write_csv(self.species_path / 'source-data.tsv', separator='\t')

  def write_isoforms_to_fasta(self):
    with open(self.fasta_path, 'w') as fasta_file:
      for row in self.sources.iter_rows(named=True):
        fasta_file.write(f">{row['Source Accession']}\n{row['Sequence']}\n")

  def run_arc(self):
    pass

  def blast_isoforms(self):
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
    subprocess.run(makeblastdb_cmd, check=True)
    subprocess.run(blastp_cmd, check=True)

  def mmseqs2_isoforms(self):
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

  def select_top_protein(self):
    alignment_cols = [
      'Query', 'Subject', '% Identity', 'Alignment Length', 'Mismatches', 'Gap Openings', 
      'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-value', 'Bit Score'
    ]
    if self.run_mmseqs2:
      alignments = pl.read_csv(self.species_path / 'alignments.tsv', separator='\t', has_header=False, new_columns=alignment_cols)
      alignments['% Identity'] = alignments['% Identity'] * 100
    else:
      alignments = pl.read_csv(self.species_path / 'alignments.csv', separator=',', has_header=False, new_columns=alignment_cols)
      alignments = alignments.with_columns(pl.col('Subject').str.split('|').list.get(1))

    query_length_map = dict(self.sources.select('Source Accession', 'Length').iter_rows())
    alignments = alignments.with_columns(
      pl.col('Query').map_dict(query_length_map).alias('Query Length')
    )
    alignments = alignments.with_columns(
      pl.col('Alignment Length').mul(pl.col('% Identity')).truediv(pl.col('Query Length')).alias('Score')
    )
    return alignments.groupby('Query').agg(pl.all().sort_by('Score').last())

  def get_protein_data(self, top_proteins):
    proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t')
    protein_data = top_proteins.join(proteome, how='left', left_on='Subject', right_on='Protein ID')
    print(protein_data)



class PeptideProcessor:
  def __init__(self):
    pass

  def process(self):
    pass


def do_assignments(taxon_id):
  species_row = active_species.row(by_predicate=pl.col('Species ID') == taxon_id)
  group = species_row[4]
  active_taxa = [int(taxon_id) for taxon_id in species_row[3].split(', ')]
  peptides = data_fetcher.get_peptides_for_species(all_peptides, active_taxa)
  sources = data_fetcher.get_sources_for_species(all_sources, peptides['Source Accession'].to_list())
  sources = sources.with_columns(pl.col('Sequence').str.len_chars().alias('Length'))
  config = {
    'taxon_id': taxon_id,
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
    for _, species_row in active_species.iterrows():
      do_assignments(species_row['Species ID'])
  else:
    do_assignments(taxon_id)