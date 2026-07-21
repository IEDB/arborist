import re
import json
import argparse
import subprocess
import os
# Cap polars' thread pool: this pipeline is RAM-bound (31GB box, no swap) and the
# per-core memory arenas measurably inflated footprint during the human step.
os.environ.setdefault("POLARS_MAX_THREADS", "4")
import polars as pl
from polars.exceptions import NoDataError
from pathlib import Path
from Bio import SeqIO
from pepmatch import Preprocessor, Matcher
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

  def generate_proteome_tsv(self):
    proteome_fasta = self.species_path / 'proteome.fasta'
    if not proteome_fasta.exists() or proteome_fasta.stat().st_size == 0: return
    regexes = {
      'protein_id': re.compile(r"\|([^|]*)\|"),         # between | and |
      'entry_name': re.compile(r"\|([^\|]*?)\s"),       # between | and space
      'protein_name': re.compile(r"\s(.+?)\sOS"),       # between space and space before OS
      'gene': re.compile(r"GN=(.*?)(?=\s[A-Z]{2}=|$)"), # between GN= and PE= or eol
      'pe_level': re.compile(r"PE=(.+?)\s"),            # between PE= and space
    }

    try:
      proteins = list(SeqIO.parse(proteome_fasta, 'fasta'))
    except ValueError:  # non-FASTA payload (e.g. a UniProt stream-error body) slipped through
      return
    if not proteins: return  # non-empty file but zero parseable records

    proteome_data = []
    for protein in proteins:
      metadata = []
      for key in regexes:
        match = regexes[key].search(str(protein.description))

        if match:
          metadata.append(match.group(1))
        else:
          if key == 'protein_id':
            metadata.append(str(protein.id))
          elif key == 'pe_level':
            metadata.append('0')
          else:
            metadata.append('')

      metadata.append(str(protein.seq))
      metadata.append(protein.id.split('|')[0])
      proteome_data.append(metadata)

    schema = {
      'Protein ID': pl.String, 'Entry Name': pl.String, 'Protein Name': pl.String,
      'Gene': pl.String, 'Protein Existence Level': pl.String,
      'Sequence': pl.String, 'Database': pl.String,
    }
    proteome = pl.DataFrame(proteome_data, schema=schema, orient='row').with_columns(
      (pl.when(pl.col('Protein ID').str.contains('-'))
      .then(pl.col('Protein ID').str.split('-').list.last())
      .otherwise(pl.lit('1')).alias('Isoform Count')),
      (pl.when(pl.col('Protein ID').str.contains('-'))
      .then(
        pl.col('Protein Existence Level').filter(pl.col('Protein ID').str.split('-').list.first() == pl.col('Protein ID')).first()
      )
      .otherwise(pl.col('Protein Existence Level')).alias('Protein Existence Level')),
      pl.col('Gene').str.replace_all(r'[ \.\(\)]', '_').alias('Gene')
    ).select([
      'Database', 'Gene', 'Protein ID', 'Entry Name', 'Isoform Count', 'Protein Name', 
      'Protein Existence Level', 'Sequence'
    ])
    proteome.write_csv(self.species_path / 'proteome.tsv', separator='\t')

  def process_species(self):
    self.generate_proteome_tsv()

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
      'alignments.csv', 'alignments.tsv', 'arc-temp-results.tsv', 'peptide-matches.tsv', 'proteome.fasta.pdb',
      'proteome.fasta.phr', 'proteome.fasta.pin', 'proteome.fasta.pjs', 'proteome.fasta.pot',
      'proteome.fasta.psq', 'proteome.fasta.ptf', 'proteome.fasta.pto', 'sources.fasta', 'proteome.tsv',
      'allergens.fasta.pdb', 'allergens.fasta.phr', 'allergens.fasta.pin', 'allergens.fasta.pjs', 'allergens.fasta.pot',
      'allergens.fasta.psq', 'allergens.fasta.ptf', 'allergens.fasta.pto'
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
    self.proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t')
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
    from ARC.classifier import SeqClassifier  # heavy git dep; only vertebrates need it
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
    else:
      old_arc_df = pl.DataFrame({'id': [''], 'class': [''], 'chain_type': [''], 'calc_mhc_allele': ['']})

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
    alignment_dtypes = {
      'Query': pl.String,
      'Subject': pl.String,
      '% Identity': pl.Float64,
      'Alignment Length': pl.Int64,
      'Mismatches': pl.Int64,
      'Gap Openings': pl.Int64,
      'Query Start': pl.Int64,
      'Query End': pl.Int64,
      'Subject Start': pl.Int64,
      'Subject End': pl.Int64,
      'E-value': pl.Float64,
      'Bit Score': pl.Float64
    }

    try:
      if self.run_mmseqs2:
        alignments = pl.read_csv(
          self.species_path / 'alignments.tsv', 
          separator='\t', 
          has_header=False, 
          new_columns=alignment_cols, 
          schema_overrides=alignment_dtypes,
          infer_schema_length=0
        )
        alignments = alignments.with_columns(pl.col('% Identity').mul(100).alias('% Identity'))
      else:
        alignments = pl.read_csv(
          self.species_path / 'alignments.csv', 
          separator=',', 
          has_header=False, 
          new_columns=alignment_cols, 
          schema_overrides=alignment_dtypes,
          infer_schema_length=0
        )
        alignments = alignments.with_columns(pl.col('Subject').str.split('|').list.get(1))

      query_length_map = dict(self.sources.select('Source Accession', 'Length').iter_rows())
      alignments = alignments.with_columns(
        pl.col('Query').cast(pl.String).replace_strict(query_length_map).cast(pl.Int32).alias('Query Length')
      )
      alignments = alignments.with_columns(
        pl.col('Alignment Length').mul(pl.col('% Identity')).truediv(pl.col('Query Length')).alias('Score')
      )

    except NoDataError: # no alignments were found
      schema = {
          'Query': pl.String, 'Subject': pl.String, '% Identity': pl.Float64,
          'Alignment Length': pl.Int32, 'Mismatches': pl.Int32, 'Gap Openings': pl.Int32,
          'Query Start': pl.Int32, 'Query End': pl.Int32, 'Subject Start': pl.Int32,
          'Subject End': pl.Int32, 'E-value': pl.Float64, 'Bit Score': pl.Float64,
          'Query Length': pl.Int32, 'Score': pl.Float64
      }
      alignments = pl.DataFrame(schema=schema)

    proteome = pl.read_csv(self.species_path / 'proteome.tsv', separator='\t').select(['Protein ID', 'Database'])
    alignments = alignments.join(
      proteome, how='left', left_on='Subject', right_on='Protein ID'
    )
    alignments = alignments.with_columns(
      (pl.col('Database') == 'sp').alias('is_reviewed')
    )

    alignments = alignments.with_columns(
      pl.col('Query').cast(pl.String).alias('Query'),
    )
    top_proteins = alignments.group_by('Query').agg(
      pl.all().sort_by(['Score', 'is_reviewed', 'Subject'], descending=[True, True, False]).first()
    )
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
    protein_data = top_proteins.join(
      self.proteome, how='left', left_on='Subject', right_on='Protein ID', coalesce=False
    )
    protein_data = protein_data.select(pl.col(
      'Query', 'Score', 'Gene_right', 'Protein ID', 'Protein Name'
    )).rename({
      'Query': 'Source Accession', 'Score': 'Source Alignment Score', 'Gene_right': 'Source Assigned Gene',
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
    self.bin_path = Path(__file__).parents[3] / 'bin'

  def process(self):
    has_allergens = self.create_allergen_fasta()
    self.preprocess_proteome()
    self.search_peptides()
    assignments = self.assign_parents()
    if has_allergens:
      assignments = self.align_to_allergens(assignments)
    assignments = self.get_protein_data(assignments)
    assignments = self.add_synonyms(assignments)
    self.write_assignments(assignments)

  def create_allergen_fasta(self):
    allergens = pl.read_csv(build_path / 'arborist' / 'allergens.tsv', separator='\t')
    species_allergens = allergens.filter(pl.col('SpeciesID') == self.taxon_id)

    if species_allergens.height == 0:
      return False

    fasta_path = self.species_path / 'allergens.fasta'
    with open(fasta_path, 'w') as f:
      for row in species_allergens.iter_rows(named=True):
        allergen_name = row['Name']
        sequence = row['Sequence']

        if not sequence:
          continue

        allergen_name_safe = allergen_name.replace(' ', '_')
        if '>' in sequence:
          fragments = re.split(r'>', sequence)
          for i, fragment in enumerate(fragments[1:], start=1):
            parts = fragment.split(':', 1)
            if len(parts) == 2:
              seq_part = parts[1]
            else:
              seq_part = parts[0]
            clean_seq = re.sub(r'[^A-Z]', '', seq_part.upper())
            if clean_seq and len(clean_seq) >= 5:
              header = f"{allergen_name_safe}_frag_{i}"
              f.write(f">{header}\n{clean_seq}\n")
        else:
          clean_seq = re.sub(r'[^A-Z]', '', sequence.upper())
          if clean_seq and len(clean_seq) >= 5:
            header = f"{allergen_name_safe}"
            f.write(f">{header}\n{clean_seq}\n")

    return True

  def preprocess_proteome(self):
    if (self.species_path / 'proteome_5mers.pepidx').exists():
      return
    Preprocessor(
      proteome = self.species_path / 'proteome.fasta',
      preprocessed_files_path = self.species_path,
    ).preprocess(k = 5)

  def search_peptides(self):
    out = self.species_path / 'peptide-matches.tsv'
    # Resume/idempotency guard, mirroring preprocess_proteome: the matcher is the
    # most expensive step and its output can be multi-GB. If a prior run already
    # produced peptide-matches.tsv (e.g. one that died downstream in
    # assign_parents), reuse it instead of recomputing. cleanup_files removes this
    # file after a successful species, so the normal weekly run is unaffected.
    if out.exists() and out.stat().st_size > 0:
      return
    peptides = [peptide for peptide in self.peptides['Sequence'].to_list() if peptide]
    matcher_args = dict(
      proteome_file=self.species_path / 'proteome.fasta',
      max_mismatches=0,
      k=5,
      preprocessed_files_path=self.species_path,
      best_match=False,
      sequence_version=False,
    )
    CHUNK = 100_000
    if len(peptides) <= CHUNK:
      Matcher(query=peptides, output_format='tsv', output_name=out, **matcher_args).match()
      return
    # Large query sets (e.g. human ~1M peptides) OOM if every match is held in
    # memory before writing. Exact matching is per-peptide independent, so we
    # match in chunks and stream each batch to disk: the concatenated rows are
    # identical to a single-shot run, with the header written once.
    first = True
    for i in range(0, len(peptides), CHUNK):
      df = Matcher(query=peptides[i:i + CHUNK], output_format='dataframe', **matcher_args).match()
      with open(out, 'wb' if first else 'ab') as f:
        df.write_csv(f, separator='\t', include_header=first)
      first = False

  def assign_parents(self):
    matches_path = self.species_path / 'peptide-matches.tsv'
    peptides_with_genes = self.peptides.filter(pl.col('Source Assigned Gene').is_not_null())
    peptides_without_genes = self.peptides.filter(pl.col('Source Assigned Gene').is_null())

    match_cols = ['Sequence', 'Gene', 'Protein ID', 'Protein Name', 'Index start', 'Index end', 'SwissProt Reviewed']
    used_cols = ['Query Sequence', 'Gene', 'Protein ID', 'Protein Name', 'Index start',
                 'Index end', 'SwissProt Reviewed', 'Protein Existence Level', 'Gene Priority']

    # The human peptide-matches.tsv is ~34M rows (~4GB). Joining it whole against
    # the ~2.9M peptide occurrences and then sorting materializes an N x M cartesian
    # on promiscuous ("hot") sequences and OOMs the 31GB box; hash-partitioning the
    # skewed key cannot bound it (one hot sequence is a single key). Instead we
    # REDUCE before joining: collapse the peptide side to compact per-key aggregates,
    # join the matches many-to-ONE, and pick the per-Sequence top-1. Peak is then
    # linear in the match count (no cartesian). The match read is still streamed in
    # hash buckets to keep peak low; the reduction is per-Sequence so bucketing is
    # exact.
    #
    # The pick is source-dependent, so the peptide side is carried, not dropped:
    #   with_genes:    per (Sequence, Source Assigned Gene) -> the SET of Source
    #                  Assigned Protein IDs; matches join many-to-one on
    #                  (Query Sequence, Gene); is_source_match = Protein ID in set.
    #   without_genes: distinct (Sequence, Source Assigned Protein ID); matches join
    #                  on the exact protein.
    # Both use an inner join, so a peptide occurrence whose source protein has no
    # match no longer injects a null row that (via nulls-first sort) would poison
    # the whole sequence's pick -- real matches win. Assigned Protein ID is
    # unchanged; only previously-blanked coordinates/names recover.
    src_gene_sets = peptides_with_genes.group_by(
      ['Sequence', 'Source Assigned Gene']
    ).agg(
      pl.col('Source Assigned Protein ID').unique().alias('_src_pids')
    )
    src_prot_keys = peptides_without_genes.select(
      ['Sequence', 'Source Assigned Protein ID']
    ).unique()

    n_buckets = 8
    tops_with_genes = []
    tops_without_genes = []
    for b in range(n_buckets):
      matches_b = (
        pl.scan_csv(matches_path, separator='\t')
        .select(used_cols)
        .filter(pl.col('Query Sequence').hash().mod(n_buckets) == b)
        .with_columns(pl.col('Gene').cast(pl.String).alias('Gene'))
        .collect()
      )

      candidates_with_genes = matches_b.join(
        src_gene_sets, how="inner",
        left_on=["Query Sequence", "Gene"],
        right_on=["Sequence", "Source Assigned Gene"],
      ).with_columns(
        pl.col('_src_pids').list.contains(pl.col('Protein ID')).alias('is_source_match'),
        pl.col('Query Sequence').alias('Sequence'),
      )
      tops_with_genes.append(candidates_with_genes.sort(
        ["Sequence", "Gene Priority", "is_source_match", "SwissProt Reviewed", "Protein Existence Level", "Protein ID", "Index start", "Index end"],
        descending=[False, True, True, True, False, False, False, False]
      ).group_by("Sequence").first().select(match_cols))

      candidates_without_genes = matches_b.join(
        src_prot_keys, how="inner",
        left_on=["Query Sequence", "Protein ID"],
        right_on=["Sequence", "Source Assigned Protein ID"],
      ).with_columns(
        pl.col('Query Sequence').alias('Sequence')
      )
      tops_without_genes.append(candidates_without_genes.sort(
        ["Sequence", "Gene Priority", "SwissProt Reviewed", "Protein Existence Level", "Protein ID", "Index start", "Index end"],
        descending=[False, True, True, False, False, False, False]
      ).group_by("Sequence").first().select(match_cols))

    top_matches_with_genes = pl.concat(tops_with_genes)
    top_matches_without_genes = pl.concat(tops_without_genes)

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

  def align_to_allergens(self, assignments):
    allergen_fasta = self.species_path / 'allergens.fasta'
    if not allergen_fasta.exists():
      return assignments

    sources_fasta = self.species_path / 'sources.fasta'
    alignments_file = self.species_path / 'allergen-alignments.csv'

    makeblastdb_cmd = [
      str(self.bin_path / 'makeblastdb'),
      '-in', str(allergen_fasta),
      '-dbtype', 'prot',
      '-out', str(allergen_fasta)
    ]

    blastp_cmd = [
      str(self.bin_path / 'blastp'),
      '-query', str(sources_fasta),
      '-db', str(allergen_fasta),
      '-evalue', '1',
      '-outfmt', '10',
      '-num_threads', '1',
      '-out', str(alignments_file)
    ]
    subprocess.run(makeblastdb_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(blastp_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    alignment_cols = [
      'Query', 'Subject', '% Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
      'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-value', 'Bit Score'
    ]

    try:
      alignments = pl.read_csv(alignments_file, separator=',', has_header=False, new_columns=alignment_cols)
      alignments = alignments.with_columns(
        pl.col('Subject')
          .str.replace(r'_frag_\d+$', '')
          .str.replace_all('_', ' ')
          .alias('Allergen Name')
      )
      high_identity = alignments.filter(pl.col('% Identity') >= 95.0)
      allergen_map = dict(high_identity.select('Query', 'Allergen Name').iter_rows())
      assignments = assignments.with_columns(
        pl.col('Source Accession').replace_strict(allergen_map, default=None).alias('Allergen Match')
      )
      assignments = assignments.with_columns(
        pl.when(pl.col('Allergen Match').is_not_null())
        .then(pl.col('Allergen Match'))
        .otherwise(pl.col('Protein Name'))
        .alias('Protein Name')
      )
      assignments = assignments.with_columns(
        pl.when(pl.col('Allergen Match').is_not_null())
        .then(pl.col('Allergen Match'))
        .otherwise(pl.col('Source Assigned Protein Name'))
        .alias('Source Assigned Protein Name')
      )
      assignments = assignments.drop('Allergen Match')

    except:
      pass

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
        data = json.load(f)
      fragment_map = {}
      for uniprot_id, fragments in data.items():
        if not fragments:
          fragment_map[uniprot_id] = ""
        else:
          fragment_map[uniprot_id] = str(fragments)
      return fragment_map
    else:
      return {}

  def add_synonyms(self, assignments):
    if (self.species_path / 'synonym-data.json').exists():
      with open(self.species_path / 'synonym-data.json', 'r') as f:
        species_synonym_data = json.load(f)
      species_synonym_data = {k: ', '.join(v) for k, v in species_synonym_data.items()}
    else:
      species_synonym_data = {}

    manual_synonyms = pl.read_csv(build_path / 'arborist' / 'manual-synonyms.tsv', separator='\t')
    manual_synonym_data = dict(manual_synonyms.select(pl.col('Accession'), pl.col('Synonyms')).iter_rows())
    for accession, synonyms in manual_synonym_data.items():
      if accession in species_synonym_data:
        combined_synonyms = set(species_synonym_data[accession].split(', ') + synonyms.split(', '))
        species_synonym_data[accession] = ', '.join(combined_synonyms)
      else:
        species_synonym_data[accession] = synonyms

    assignments = assignments.with_columns(
      pl.col('Protein ID').replace_strict(species_synonym_data, default='').alias('Assigned Protein Synonyms'),
    )
    assignments = assignments.with_columns(
      pl.concat_str(
        [
          pl.when(pl.col('Assigned Protein Synonyms') != '').then(pl.col('Assigned Protein Synonyms')), 
          pl.when(pl.col('Source Assigned Gene') != '').then(pl.col('Source Assigned Gene')),
          pl.when(pl.col('Entry Name') != '').then(pl.col('Entry Name'))
        ], 
        separator=', ', 
        ignore_nulls=True
      ).alias('Assigned Protein Synonyms')
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
  proteome_fasta = species_path / 'proteome.fasta'
  # Trigger on-demand selection when proteome.fasta is MISSING or EMPTY, not
  # just when the species directory is absent. Upstream Makefile stages create
  # build/species/<id>/ (epitopes/taxa/sources) for a newly-active species
  # before its proteome exists; guarding on the directory alone skipped
  # selection and left a dir-but-no-fasta hole that crashed the whole run.
  if not proteome_fasta.exists() or proteome_fasta.stat().st_size == 0:
    print(f'No proteome detected for {species_name} (ID: {taxon_id}), selecting best one...')
    data_fetcher = DataFetcher(build_path)
    peptides = data_fetcher.get_peptides_for_species(all_peptides, active_taxa)
    proteome_selector = ProteomeSelector(
      taxon_id, species_name, group, peptides, build_path
    )
    proteome_selector.select()

def check_for_skips(taxon_id):
  # Belt-and-suspenders: treat a missing proteome.fasta the same as an empty
  # one -> skip this species cleanly, never raise. On-demand selection can
  # legitimately produce nothing (no UniProt candidates, no orphan proteins),
  # and one bad species must never kill the entire run again.
  proteome_fasta = build_path / 'species' / str(taxon_id) / 'proteome.fasta'
  if not proteome_fasta.exists() or proteome_fasta.stat().st_size == 0:
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
      try:
        all_species_data = pl.concat([all_species_data, species_data])
      except:
        species_data = species_data.with_columns(
          pl.lit(species_path.name).alias('Species ID'),
          pl.col('Species Name').cast(pl.String).alias('Proteome Label'),
        )
        all_species_data = pl.concat([all_species_data, species_data])

  all_assignments.write_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')
  all_source_data.write_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')
  all_species_data.write_csv(build_path / 'arborist' / 'all-species-data.tsv', separator='\t')

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-t', '--taxon_id', type=int, help='Taxon ID of the species to process.',
  )
  parser.add_argument(
    '-b', '--build_path', type=str, help='Path for all Arborist build files.',
    default=Path(__file__).parents[3] / 'build'
  )
  parser.add_argument(
    '-n', '--num_threads', type=int, help='Number of threads to use.',
    default=1
  )
  args = parser.parse_args()

  taxon_id = args.taxon_id
  build_path = Path(args.build_path)
  all_species = not bool(taxon_id)

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
