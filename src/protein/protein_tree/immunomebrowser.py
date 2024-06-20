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
      .then(pl.lit('http://www.uniprot.org/uniprot/') + pl.col('Parent Protein ID'))
      .otherwise(
        pl.lit('https://ontology.iedb.org/taxon-protein/') + 
        pl.col('Species Taxon ID').cast(pl.String) + pl.lit('-other')
      )
      .alias('Parent IRI'),
    pl.when(pl.col('Parent Protein ID').is_not_null())
      .then(pl.lit('UniProt'))
      .otherwise(pl.lit('IEDB'))
      .alias('Parent Protein Database')
  ).unique(subset=['Source Accession'])
  source_parents = all_parent_data.select(
    'Source ID', 'Source Accession', 'Database', 'Name', 'Aliases', 'Assigned Protein Synonyms',
    'Organism ID', 'Organism Name', 'Species Taxon ID', 'Species Name', 'Proteome ID',
    'Proteome Label', 'Protein Strategy', 'Parent IRI', 'Parent Protein Database',
    'Parent Protein ID', 'Assigned Protein Length', 'Assigned Protein Sequence',
    'Source Assigned Gene'
  ).rename({
    'Source Accession': 'Accession', 'Assigned Protein Synonyms': 'Synonyms',
    'Organism ID': 'Taxon ID', 'Organism Name': 'Taxon Name', 'Species Taxon ID': 'Species ID',
    'Species Name': 'Species Label', 'Parent Protein ID': 'Parent Protein Accession', 
    'Assigned Protein Length': 'Parent Sequence Length', 'Assigned Protein Sequence': 'Sequence',
    'Source Assigned Gene': 'Parent Protein Gene'
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
    'Parent Protein Name': 'Name', 'Assigned Protein Name': 'Title', 
    'Assigned Protein Sequence': 'Sequence'
  })
  parent_proteins.write_csv(build_path / 'arborist' / 'parent-proteins.tsv', separator='\t')


class EpitopeMapper:
  def __init__(self, assignments, num_threads):
    self.assignments = assignments
    self.num_threads = num_threads
    self.bin_path = Path(__file__).parents[3] / 'bin'
    self.blast_temp = build_path / 'blast_temp'

  def make_epitope_mappings(self):
    exact_matches, non_exact_linear, non_exact_discontinuous = self.split_assignments()
    exact_mappings = self.make_exact_mappings(exact_matches)
    linear_mappings = self.make_linear_mappings(non_exact_linear)
    discontinous_mappings = self.make_discontinuous_mappings(non_exact_discontinuous)
    self.combine_and_write_mappings(exact_mappings, linear_mappings, discontinous_mappings)

  def split_assignments(self):
    exact_matches = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_not_null()) & 
      (pl.col('Assigned Protein ID').is_not_null())
    )
    non_exact_linear = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_null()) & 
      (pl.col('Assigned Protein ID').is_not_null()) &
      (~pl.col('Epitope Sequence').str.contains('0|1|2|3|4|5|6|7|8|9')) &
      (pl.col('Epitope ID').is_not_null())
    )
    non_exact_discontinuous = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_null()) & 
      (pl.col('Assigned Protein ID').is_not_null()) &
      (pl.col('Epitope Sequence').str.contains('0|1|2|3|4|5|6|7|8|9')) &
      (pl.col('Epitope ID').is_not_null())
    )
    return exact_matches, non_exact_linear, non_exact_discontinuous
  
  def combine_mappings(self, exact_mappings, linear_mappings, discontinous_mappings):
    return pl.concat([exact_mappings, linear_mappings, discontinous_mappings])

  def make_exact_mappings(self, exact_matches):
    exact_mappings = exact_matches.with_columns(
      (pl.col('Epitope Sequence')).alias('parent_seq'),
      (pl.lit(1.0)).alias('identity_alignment'),
      (pl.lit(1.0)).alias('similarity_alignment'),
      (pl.lit(0.0)).alias('gaps_source_alignment'),
      (pl.lit(0.0)).alias('gaps_parent_alignment'),
      (pl.lit(0.0)).alias('all_gaps'),
      (pl.col('Epitope Sequence')).alias('source_alignment'),
      (pl.col('Epitope Sequence')).alias('parent_alignment'),
      (pl.col('Epitope Sequence')).alias('parent_alignment_modified')
    ).select(
      'Epitope ID', 'Epitope Sequence', 'Source Starting Position', 'Source Ending Position',
      'Source Accession', 'Assigned Protein ID', 'parent_seq', 'Assigned Protein Starting Position', 
      'Assigned Protein Ending Position', 'identity_alignment', 'similarity_alignment',
      'gaps_source_alignment', 'gaps_parent_alignment', 'all_gaps', 'source_alignment',
      'parent_alignment', 'parent_alignment_modified'
    ).rename({
      'Epitope ID': 'epitope_id', 'Epitope Sequence': 'epitope_seq',
      'Source Starting Position': 'epitope_start', 'Source Ending Position': 'epitope_end',
      'Source Accession': 'source_accession', 'Assigned Protein ID': 'parent_accession',
      'Assigned Protein Starting Position': 'parent_start', 
      'Assigned Protein Ending Position': 'parent_end'
    }).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })
    return exact_mappings

  def make_linear_mappings(self, non_exact_linear):
    self.blast_linear_peptides(non_exact_linear)
    blast_cols = [
      'Query', 'Subject', 'Query Sequence', 'Subject Sequence', '% Identity', 'Alignment Length',
      'Gaps', 'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-value'
    ]
    blast_results = pl.read_csv(
      self.blast_temp / 'blast_results.csv', separator=',', has_header=False, new_columns=blast_cols
    )
    top_alignments = non_exact_linear.with_columns(
      pl.col('Epitope ID').cast(pl.Int64).alias('Query')
    ).join(blast_results, how='left', on='Query', coalesce=True).filter(
      pl.col('Subject') == pl.col('Assigned Protein ID')
    ).group_by(['Source Accession', 'Epitope ID']).agg(pl.all().sort_by('% Identity').last())
    
    linear_mappings = top_alignments.with_columns(
      (pl.col('% Identity')).alias('identity_alignment'),
      (pl.col('% Identity')).alias('similarity_alignment'),
      (pl.col('Query Sequence').str.count_matches('-') / pl.col('Alignment Length')).alias('gaps_source_alignment'),
      (pl.col('Subject Sequence').str.count_matches('-') / pl.col('Alignment Length')).alias('gaps_parent_alignment'),
      (pl.col('Gaps') / pl.col('Alignment Length')).alias('all_gaps'),
      (pl.col('Query Sequence')).alias('source_alignment'),
      (pl.col('Subject Sequence')).alias('parent_alignment'),
      (pl.struct(['Query Sequence', 'Subject Sequence']).map_elements(
        lambda x: self.modify_parent_alignment(x['Query Sequence'], x['Subject Sequence']), 
        return_dtype=pl.String
      )).alias('parent_alignment_modified')
    ).select(
      'Epitope ID', 'Epitope Sequence', 'Source Starting Position', 'Source Ending Position',
      'Source Accession', 'Assigned Protein ID', 'Subject Sequence', 'Subject Start', 'Subject End',
      'identity_alignment', 'similarity_alignment', 'gaps_source_alignment', 'gaps_parent_alignment',
      'all_gaps', 'source_alignment', 'parent_alignment', 'parent_alignment_modified'
    ).rename({
      'Epitope ID': 'epitope_id', 'Epitope Sequence': 'epitope_seq', 
      'Source Starting Position': 'epitope_start', 'Source Ending Position': 'epitope_end',
      'Source Accession': 'source_accession', 'Assigned Protein ID': 'parent_accession',
      'Subject Sequence': 'parent_seq', 'Subject Start': 'parent_start', 'Subject End': 'parent_end',
     }).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64, 
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })
    return linear_mappings
  
  def blast_linear_peptides(self, non_exact_linear):
    self.blast_temp.mkdir(exist_ok=True)
    self.make_blast_db(non_exact_linear)
    self.make_epitope_file(non_exact_linear)
    cmd = [
      str(self.bin_path / 'blastp'),
      '-query', str(self.blast_temp / 'epitopes.fasta'),
      '-db', str(self.blast_temp / 'proteins.fasta'),
      '-outfmt', '10 qseqid sseqid qseq sseq pident length gaps qstart qend sstart send evalue',
      '-num_threads', str(self.num_threads),
      '-out', str(self.blast_temp / 'blast_results.csv')
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  def make_blast_db(self, non_exact_linear):
    proteins = non_exact_linear.unique(subset=['Assigned Protein ID'])
    with open(str(self.blast_temp / 'proteins.fasta'), 'w') as f:
      for row in proteins.iter_rows(named=True):
        f.write(f'>{row["Assigned Protein ID"]}\n{row["Assigned Protein Sequence"]}\n')
    cmd = [
      str(self.bin_path / 'makeblastdb'),
      '-in', str(self.blast_temp / 'proteins.fasta'),
      '-dbtype', 'prot',
      '-out', str(self.blast_temp / 'proteins.fasta')
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  def make_epitope_file(self, non_exact_linear):
    with open(str(self.blast_temp / 'epitopes.fasta'), 'w') as f:
      for row in non_exact_linear.iter_rows(named=True):
        f.write(f'>{int(row["Epitope ID"])}\n{row["Epitope Sequence"]}\n')

  def make_discontinuous_mappings(self, non_exact_discontinuous):
    parent_seqs, perc_identities, source_alignments, parent_alignments, parent_alignments_modified = [], [], [], [], []
    for row in non_exact_discontinuous.iter_rows(named=True):
      seq, parent_seq, perc_identity, parent_alignment = self.map_discontinuous_epitopes(
        row['Epitope Sequence'], row['Assigned Protein Sequence']
      )
      parent_seqs.append(parent_seq)
      perc_identities.append(perc_identity)
      parent_alignments.append(parent_alignment)

      collapsed_seq = ''.join([x[0] for x in seq.split(',')])
      source_alignments.append(collapsed_seq)
      parent_alignments_modified.append(self.modify_parent_alignment(collapsed_seq, parent_alignment))

    discontinous_mappings = non_exact_discontinuous.with_columns(
      pl.Series(name='parent_seq', values=parent_seqs),
      pl.Series(name='identity_alignment', values=perc_identities),
      pl.Series(name='similarity_alignment', values=perc_identities),
      (pl.lit(0.0)).alias('gaps_source_alignment'),
      (pl.lit(0.0)).alias('gaps_parent_alignment'),
      (pl.lit(0.0)).alias('all_gaps'),
      pl.Series(name='source_alignment', values=source_alignments),
      pl.Series(name='parent_alignment', values=parent_alignments),
      pl.Series(name='parent_alignment_modified', values=parent_alignments_modified)
    ).with_columns(
      pl.col('parent_seq').str.split(', ').list.get(0).str.slice(1, 5).alias('parent_start'),
      pl.col('parent_seq').str.split(', ').list.get(-1).str.slice(1, 5).alias('parent_end')
    ).select(
      'Epitope ID', 'Epitope Sequence', 'Source Starting Position', 'Source Ending Position',
      'Source Accession', 'Assigned Protein ID', 'parent_seq', 'parent_start', 'parent_end',
      'identity_alignment', 'similarity_alignment', 'gaps_source_alignment', 'gaps_parent_alignment',
      'all_gaps', 'source_alignment', 'parent_alignment', 'parent_alignment_modified'
    ).rename({
      'Epitope ID': 'epitope_id', 'Epitope Sequence': 'epitope_seq',
      'Source Starting Position': 'epitope_start', 'Source Ending Position': 'epitope_end',
      'Source Accession': 'source_accession', 'Assigned Protein ID': 'parent_accession'
    }).filter(
      pl.col('parent_alignment') != ''
    ).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })
    return discontinous_mappings
  
  def map_discontinuous_epitopes(self, seq, protein_seq):
    seq = seq.replace(',  ', ', ').replace(' ,', ', ').replace(', ', ',').replace(' ', ',').replace(',,', ',')
    seq = seq[0:-1] if seq[-1] == ',' else seq
    split_seq = seq.split(',')
    seq_len = len(split_seq)
    residue_matches = 0
    parent_seq = ''
    for residue in split_seq:
      try:
        pos = int(residue[1:])-1
        if residue[0] == protein_seq[pos]:
          parent_seq += f'{residue[0]}{pos+1}, '
          residue_matches += 1
        else:
          parent_seq += f'{protein_seq[pos]}{pos+1}, '
      except IndexError:
        parent_seq += f'{'X'}{pos+1}, '
      except ValueError:
        continue
    parent_seq = parent_seq[:-2]
    if parent_seq == '':
      return seq, parent_seq, 0.0, ''

    perc_identity = residue_matches / seq_len
    parent_alignment = ''.join([x[0] for x in parent_seq.split(', ')])
    return seq, parent_seq, perc_identity, parent_alignment

  def modify_parent_alignment(self, qseq, sseq):
    if len(qseq) != len(sseq):
      return ''
    seq = ''
    for i in range(len(qseq)):
      if qseq[i] != sseq[i]:
        seq += sseq[i].lower()
      else:
        seq += sseq[i]
    return seq

  def combine_and_write_mappings(self, exact_mappings, linear_mappings, discontinous_mappings):
    combined_mappings = pl.concat([exact_mappings, linear_mappings, discontinous_mappings])
    combined_mappings = combined_mappings.filter(pl.col('epitope_id').is_not_null())
    combined_mappings.write_csv(build_path / 'arborist' / 'epitope-mappings.tsv', separator='\t')

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-n', '--num_threads', type=int, default=1, help='Number of threads to use.')
  args = parser.parse_args()

  build_path = Path(__file__).parents[3] / 'build'
  assignments = pl.read_csv(build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t')
  source_data = pl.read_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')
  species_data = pl.read_csv(build_path / 'arborist' / 'all-species-data.tsv', separator='\t')

  all_parent_data = get_all_parent_data(assignments, source_data, species_data)
  make_source_parents(all_parent_data)
  make_parent_proteins(all_parent_data)

  epitope_mapping = EpitopeMapper(assignments, args.num_threads)
  epitope_mapping.make_epitope_mappings()