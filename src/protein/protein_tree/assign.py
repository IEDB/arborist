import json
import argparse
import subprocess
import polars as pl
from pathlib import Path
from pepmatch import Preprocessor, Matcher
from ARC.classifier import SeqClassifier
from protein_tree.data_fetch import DataFetcher
from protein_tree.select_proteome import ProteomeSelector


class AssignmentHandler:
  def __init__(self, taxon_id, species_name, group, peptides, sources, num_threads):
    self.taxon_id = taxon_id
    self.species_name = species_name
    self.group = group
    self.peptides = peptides
    self.sources = sources
    self.num_threads = num_threads
    self.species_path = build_path / 'species' / str(self.taxon_id)

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
    files_to_remove = [
      'alignments.csv', 'alignments.tsv', 'arc-temp-results.tsv', 'proteome.fasta.pdb', 
      'proteome.fasta.phr', 'proteome.fasta.pin', 'proteome.fasta.pjs', 'proteome.fasta.pot',
      'proteome.fasta.psq', 'proteome.fasta.ptf', 'proteome.fasta.pto', 'sources.fasta'
    ]
    for file in files_to_remove:
      file_path = self.species_path / file
      if file_path.exists():
        file_path.unlink()
    subprocess.run(['rm', '-rf', str(self.species_path / 'tmp')])

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
    top_proteins = self.assign_manuals(top_proteins)
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
    subprocess.run(mmseqs2_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

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
      hmmer_path=str(self.bin_path)+'/',
      blast_path=str(self.bin_path)+'/'
    ).classify_seqfile(f'{self.species_path}/sources.fasta')

    temp_results = pl.read_csv(temp_results_path, separator='\t')
    if temp_results.shape[0] == 0:
      return old_arc_df

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
    try:
      if self.run_mmseqs2:
        alignments = pl.read_csv(self.species_path / 'alignments.tsv', separator='\t', has_header=False, new_columns=alignment_cols, infer_schema_length=10000)
        alignments = alignments.with_columns(pl.col('% Identity').mul(100).alias('% Identity'))
      else:
        alignments = pl.read_csv(self.species_path / 'alignments.csv', separator=',', has_header=False, new_columns=alignment_cols, infer_schema_length=10000)
        alignments = alignments.with_columns(pl.col('Subject').str.split('|').list.get(1))

      query_length_map = dict(self.sources.select('Source Accession', 'Length').iter_rows())
      alignments = alignments.with_columns(
        pl.col('Query').cast(pl.String).replace_strict(query_length_map).cast(pl.Int32).alias('Query Length')
      )
      alignments = alignments.with_columns(
        pl.col('Alignment Length').mul(pl.col('% Identity')).truediv(pl.col('Query Length')).alias('Score')
      )

    except pl.NoDataError: # no alignments were found
      alignments = pl.DataFrame({col: [] for col in alignment_cols})
      alignments = alignments.with_columns(
        pl.lit(0).alias('Query Length'),
        pl.lit(0).alias('Score')
      )
    alignments = alignments.with_columns(
      pl.col('Query').cast(pl.String).alias('Query'),
    )
    top_proteins = alignments.group_by('Query').agg(pl.all().sort_by('Score').last())
    missing_sources = self.sources.filter(
      ~pl.col('Source Accession').is_in(top_proteins['Query'].to_list())
    )
    if missing_sources.shape[0] != 0:
      missing_sources = pl.DataFrame({'Query': missing_sources['Source Accession'].to_list()})
      top_proteins = top_proteins.join(missing_sources, how='full', on='Query', coalesce=True)
    
    return top_proteins

  def assign_manuals(self, top_proteins):
    manual_parents = pl.read_csv(build_path / 'arborist' / 'manual-parents.tsv', separator='\t')
    top_proteins = top_proteins.join(
      manual_parents, how='left', left_on='Query', right_on='Accession', coalesce=True
    ).with_columns(
      pl.when(pl.col('Parent Accession').is_not_null())
      .then(pl.col('Parent Accession')).otherwise(pl.col('Subject')).alias('Subject')
    )
    return top_proteins

  def get_protein_data(self, top_proteins):
    proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t')
    protein_data = top_proteins.join(
      proteome, how='left', left_on='Subject', right_on='Protein ID', coalesce=False
    )
    allergy_data = pl.read_json(build_path / 'arborist' / 'allergens.json')
    protein_data = protein_data.join(
      allergy_data, left_on='Protein ID', right_on='uniprot_id', how='left', coalesce=True
    )
    protein_data = protein_data.with_columns(
      pl.when(pl.col('allergen_name').is_not_null())
      .then(pl.col('allergen_name')).otherwise(pl.col('Protein Name')).alias('Protein Name')
    )
    protein_data = protein_data.with_columns(
      pl.when(pl.col('mapped_id').is_not_null())
      .then(pl.col('mapped_id')).otherwise(pl.col('Protein ID')).alias('Protein ID')
    )
    protein_data = protein_data.select(pl.col(
      'Query', 'Score', 'Gene', 'Protein ID', 'Protein Name'
    )).rename({
      'Query': 'Source Accession', 'Score': 'Source Alignment Score', 'Gene': 'Source Assigned Gene',
      'Protein ID': 'Source Assigned Protein ID', 'Protein Name': 'Source Assigned Protein Name'}
    )
    protein_data = protein_data.with_columns(
      pl.col('Source Assigned Gene').cast(pl.String).alias('Source Assigned Gene')
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
    assignments = self.handle_allergens(assignments)
    assignments = self.get_protein_data(assignments)
    assignments = self.add_synonyms(assignments)
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
    peptides = [peptide for peptide in self.peptides['Sequence'].to_list() if peptide]
    Matcher(
      query=peptides,
      proteome_file=self.species_path / 'proteome.fasta',
      max_mismatches=0,
      k=5,
      preprocessed_files_path=self.species_path,
      best_match=False, 
      output_format='tsv',
      output_name=self.species_path / 'peptide-matches.tsv',
      sequence_version=False
    ).match()
  
  def assign_parents(self):
    matches = pl.read_csv(self.species_path / 'peptide-matches.tsv', separator='\t').with_columns(
      pl.col('Gene').cast(pl.String).alias('Gene'),
    )
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
    assignments = assignments.drop('Sequence_right', 'Gene', 'SwissProt Reviewed')
    assignments = assignments.unique(subset=['Sequence', 'Source Accession'])
    assignments = assignments.with_columns(
      pl.col('Protein ID').fill_null(pl.col('Source Assigned Protein ID')),
      pl.col('Protein Name').fill_null(pl.col('Source Assigned Protein Name'))
    )
    return assignments

  def get_protein_data(self, assignments):
    proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t')
    proteome = proteome.select(pl.col('Database', 'Entry Name', 'Protein ID', 'Sequence'))
    proteome = proteome.rename({'Sequence': 'Assigned Protein Sequence'})
    fragments = self.get_fragment_data()
    assignments = assignments.join(
      proteome, how='left', on='Protein ID', coalesce=True
    )
    assignments = assignments.with_columns(
      (pl.col('Database') == 'sp').alias('Assigned Protein Review Status'),
      pl.col('Assigned Protein Sequence').str.len_chars().alias('Assigned Protein Length'),
      pl.col('Protein ID').replace_strict(fragments, default='').alias('Assigned Protein Fragments'),
      pl.lit(str(self.taxon_id)).alias('Species Taxon ID'),
      pl.lit(self.species_name).alias('Species Name')
    )
    return assignments
  
  def get_fragment_data(self):
    if (self.species_path / 'fragment-data.json').exists():
      with open(self.species_path / 'fragment-data.json', 'r') as f:
        return json.load(f)
    else:
      return {}
      
  def handle_allergens(self, assignments):
    allergy_data = pl.read_json(build_path / 'arborist' / 'allergens.json')
    assignments = assignments.join(
      allergy_data, left_on='Protein ID', right_on='uniprot_id', how='left', coalesce=True
    )
    assignments = assignments.with_columns(
      pl.when(pl.col('allergen_name').is_not_null())
      .then(pl.col('allergen_name')).otherwise(pl.col('Protein Name')).alias('Protein Name')
    )
    assignments = assignments.with_columns(
      pl.when(pl.col('mapped_id').is_not_null())
      .then(pl.col('mapped_id')).otherwise(pl.col('Protein ID')).alias('Protein ID')
    )
    return assignments
  
  def add_synonyms(self, assignments):
    if (self.species_path / 'synonym-data.json').exists():
      with open(self.species_path / 'synonym-data.json', 'r') as f:
        synonym_data = json.load(f)
      synonym_data = {k: ', '.join(v) for k, v in synonym_data.items()}
    else:
      synonym_data = {}
    assignments = assignments.with_columns(
      pl.col('Protein ID').replace_strict(synonym_data, default='').alias('Assigned Protein Synonyms'),
    )
    assignments = assignments.with_columns(
      pl.when(pl.col('Assigned Protein Synonyms') != "")
      .then(pl.concat_str(
        [pl.col('Assigned Protein Synonyms'), pl.col('Entry Name')], separator=', ')
      )
      .otherwise(pl.col('Entry Name'))
      .alias('Assigned Protein Synonyms')
    )
    return assignments

  def write_assignments(self, assignments):
    assignments = assignments.rename({
      'Sequence': 'Epitope Sequence',
      'Starting Position': 'Source Starting Position',
      'Ending Position': 'Source Ending Position',
      'Protein ID': 'Assigned Protein ID',
      'Protein Name': 'Assigned Protein Name',
      'Entry Name': 'Assigned Protein Entry Name',
      'Index start': 'Assigned Protein Starting Position',
      'Index end': 'Assigned Protein Ending Position',
    })
    col_order = [
      'Species Taxon ID', 'Species Name', 'Organism ID', 'Organism Name','Source Accession', 
      'Source Alignment Score', 'Source Assigned Gene', 'Source Assigned Protein ID', 
      'Source Assigned Protein Name', 'ARC Assignment', 'Epitope ID', 'Epitope Sequence', 
      'Source Starting Position', 'Source Ending Position', 'Assigned Protein ID', 
      'Assigned Protein Name', 'Assigned Protein Entry Name', 'Assigned Protein Review Status',
      'Assigned Protein Starting Position', 'Assigned Protein Ending Position', 
      'Assigned Protein Sequence', 'Assigned Protein Length', 'Assigned Protein Fragments',
      'Assigned Protein Synonyms'
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
  check_for_proteome(taxon_id, active_taxa, species_name, group)
  print(f'Assigning peptides for {species_name} (ID: {taxon_id})')
  skip = check_for_skips(taxon_id)
  if skip:
    return
  assignment_handler = AssignmentHandler(**config)
  assignment_handler.process_species()
  assignment_handler.cleanup_files()

def check_for_proteome(taxon_id, active_taxa, species_name, group):
  species_path = build_path / 'species' / str(taxon_id)
  if not species_path.exists():
    Fetcher = DataFetcher(build_path)
    peptides_df = Fetcher.get_peptides_for_species(all_peptides, active_taxa).to_pandas()
    Selector = ProteomeSelector(
      taxon_id, species_name, group, build_path
    )
    Selector.select_best_proteome(peptides_df)
    Selector.proteome_to_tsv()

def check_for_skips(taxon_id):
  if (build_path / 'species' / str(taxon_id) / 'proteome.fasta').stat().st_size == 0:
    print(f'Proteome is empty, skipping.')
    return True

def combine_data():
  print('Combining assignment, source, and species data.')
  all_assignments = pl.DataFrame()
  all_source_data = pl.DataFrame()
  all_species_data = pl.DataFrame()

  for species_path in sorted((build_path / 'species').iterdir(), key=lambda x: int(x.name)):
    if (species_path / 'peptide-assignments.tsv').exists():
      assignments = pl.read_csv(species_path / 'peptide-assignments.tsv', separator='\t', infer_schema_length=0)
      all_assignments = pl.concat([all_assignments, assignments])
    if (species_path / 'source-data.tsv').exists():
      source_data = pl.read_csv(species_path / 'source-data.tsv', separator='\t', infer_schema_length=0)
      all_source_data = pl.concat([all_source_data, source_data])
    if (species_path / 'species-data.tsv').exists():
      species_data = pl.read_csv(species_path / 'species-data.tsv', separator='\t', infer_schema_length=0)
      all_species_data = pl.concat([all_species_data, species_data])

  all_assignments.write_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')
  all_source_data.write_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')
  all_species_data.write_csv(build_path / 'arborist' / 'all-species-data.tsv', separator='\t')

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
    for row in active_species.rows(named=True):
      do_assignments(row['Species ID'])
    combine_data()
  else:
    do_assignments(taxon_id)