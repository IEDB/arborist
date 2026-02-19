import argparse
import polars as pl
from pathlib import Path
from Bio.Align import PairwiseAligner, substitution_matrices


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
    source_data, how='left', on='Source Accession'
  ).join(
    species_data, how='left', on='Species Taxon ID', coalesce=True
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
    'Organism ID_right', 'Organism Name_right', 'Species Taxon ID', 'Species Name', 'Proteome ID',
    'Proteome Label', 'Protein Strategy', 'Parent IRI', 'Parent Protein Database',
    'Parent Protein ID', 'Assigned Protein Length', 'Assigned Protein Sequence',
    'Source Assigned Gene'
  ).rename({
    'Source Accession': 'Accession', 'Assigned Protein Synonyms': 'Synonyms',
    'Organism ID_right': 'Taxon ID', 'Organism Name_right': 'Taxon Name', 'Species Taxon ID': 'Species ID',
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
  def __init__(self, assignments, source_data):
    self.assignments = assignments
    self.source_data = source_data
    self.source_seq_map = dict(
      source_data.filter(
        pl.col('Sequence').is_not_null() & (pl.col('Sequence') != '')
      ).select('Source Accession', 'Sequence').iter_rows()
    )
    self.aligner = self._init_aligner()
    self.alignment_cache = {}

  def _init_aligner(self):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = 'global'
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner

  def make_epitope_mappings(self):
    print(f'Total assignment rows: {self.assignments.height}')

    exact, non_exact = self._split_assignments()
    tier2_linear, tier2_discontinuous, tier3 = self._split_non_exact(non_exact)

    print(f'Tier 1 (exact match): {exact.height}')
    print(f'Tier 2 linear (source->parent): {tier2_linear.height}')
    print(f'Tier 2 discontinuous: {tier2_discontinuous.height}')
    print(f'Tier 3 (peptide->parent, no source seq): {tier3.height}')

    exact_mappings = self._make_exact_mappings(exact)
    print(f'Exact mappings produced: {exact_mappings.height}')

    tier2_linear_mappings = self._make_tier2_linear_mappings(tier2_linear)
    print(f'Tier 2 linear mappings produced: {tier2_linear_mappings.height}')

    tier2_disc_mappings = self._make_tier2_discontinuous_mappings(tier2_discontinuous)
    print(f'Tier 2 discontinuous mappings produced: {tier2_disc_mappings.height}')

    tier3_mappings = self._make_tier3_mappings(tier3)
    print(f'Tier 3 mappings produced: {tier3_mappings.height}')

    self._combine_and_write(exact_mappings, tier2_linear_mappings, tier2_disc_mappings, tier3_mappings)

  def _split_assignments(self):
    exact = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_not_null()) &
      (pl.col('Assigned Protein ID').is_not_null())
    )
    non_exact = self.assignments.filter(
      (pl.col('Assigned Protein Starting Position').is_null()) &
      (pl.col('Assigned Protein ID').is_not_null()) &
      (pl.col('Epitope ID').is_not_null())
    )
    return exact, non_exact

  def _split_non_exact(self, non_exact):
    has_source_seq = non_exact.with_columns(
      pl.col('Source Accession').is_in(list(self.source_seq_map.keys())).alias('_has_source')
    )

    with_source = has_source_seq.filter(pl.col('_has_source')).drop('_has_source')
    without_source = has_source_seq.filter(~pl.col('_has_source')).drop('_has_source')

    is_discontinuous = pl.col('Epitope Sequence').str.contains(r'[0-9]')

    tier2_linear = with_source.filter(~is_discontinuous)
    tier2_discontinuous = with_source.filter(is_discontinuous)

    tier3_linear = without_source.filter(~is_discontinuous)
    tier3_discontinuous = without_source.filter(is_discontinuous)

    if tier3_discontinuous.height > 0:
      print(f'  Tier 3 discontinuous (no source seq, unmappable): {tier3_discontinuous.height}')

    tier3 = tier3_linear
    return tier2_linear, tier2_discontinuous, tier3

  def _make_exact_mappings(self, exact):
    return exact.with_columns(
      pl.col('Epitope ID').cast(pl.String).str.replace(r"\.0$", "").alias('Epitope ID'),
    ).select(
      pl.col('Epitope ID').alias('epitope_id'),
      pl.col('Epitope Sequence').alias('epitope_seq'),
      pl.col('Source Starting Position').alias('epitope_start'),
      pl.col('Source Ending Position').alias('epitope_end'),
      pl.col('Source Accession').alias('source_accession'),
      pl.col('Assigned Protein ID').alias('parent_accession'),
      pl.col('Epitope Sequence').alias('parent_seq'),
      pl.col('Assigned Protein Starting Position').alias('parent_start'),
      pl.col('Assigned Protein Ending Position').alias('parent_end'),
      pl.lit(1.0).alias('identity_alignment'),
      pl.lit(1.0).alias('similarity_alignment'),
      pl.lit(0.0).alias('gaps_source_alignment'),
      pl.lit(0.0).alias('gaps_parent_alignment'),
      pl.lit(0.0).alias('all_gaps'),
      pl.col('Epitope Sequence').alias('source_alignment'),
      pl.col('Epitope Sequence').alias('parent_alignment'),
      pl.col('Epitope Sequence').alias('parent_alignment_modified'),
    ).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })

  def _get_source_parent_alignment(self, source_accession, parent_id, parent_seq):
    cache_key = (source_accession, parent_id)
    if cache_key in self.alignment_cache:
      return self.alignment_cache[cache_key]

    source_seq = self.source_seq_map.get(source_accession)
    if not source_seq or not parent_seq:
      self.alignment_cache[cache_key] = None
      return None

    try:
      alignments = self.aligner.align(source_seq, parent_seq)
      top = alignments[0]
      qseq = str(top[0])
      sseq = str(top[1])

      pos_map = {}
      qpos = 1
      spos = 1
      for qres, sres in zip(qseq, sseq):
        if qres != '-' and sres != '-':
          pos_map[qpos] = (spos, sres)
          qpos += 1
          spos += 1
        elif qres == '-':
          spos += 1
        elif sres == '-':
          pos_map[qpos] = None
          qpos += 1

      result = {
        'pos_map': pos_map,
        'qseq': qseq,
        'sseq': sseq,
      }
      self.alignment_cache[cache_key] = result
      return result
    except Exception as e:
      print(f'  WARNING: Alignment failed for {source_accession} -> {parent_id}: {e}')
      self.alignment_cache[cache_key] = None
      return None

  def _map_linear_via_source_alignment(self, row, alignment_data):
    pos_map = alignment_data['pos_map']
    epitope_seq = row['Epitope Sequence']
    src_start = int(row['Source Starting Position'])
    src_end = int(row['Source Ending Position'])

    parent_positions = []
    parent_residues = []
    source_residues = []

    for src_pos in range(src_start, src_end + 1):
      mapped = pos_map.get(src_pos)
      source_residues.append(epitope_seq[src_pos - src_start] if (src_pos - src_start) < len(epitope_seq) else '?')
      if mapped is not None:
        parent_positions.append(mapped[0])
        parent_residues.append(mapped[1])
      else:
        parent_positions.append(None)
        parent_residues.append('-')

    valid_positions = [p for p in parent_positions if p is not None]
    if not valid_positions:
      return None

    parent_start = min(valid_positions)
    parent_end = max(valid_positions)

    source_alignment = ''.join(source_residues)
    parent_alignment = ''.join(parent_residues)

    matches = sum(1 for s, p in zip(source_residues, parent_residues) if s == p and p != '-')
    total = len(source_residues)
    identity = matches / total if total > 0 else 0.0

    gaps_in_parent = sum(1 for p in parent_residues if p == '-')
    all_gaps = gaps_in_parent / total if total > 0 else 0.0

    parent_alignment_modified = self._modify_parent_alignment(source_alignment, parent_alignment)

    return {
      'parent_seq': parent_alignment,
      'parent_start': parent_start,
      'parent_end': parent_end,
      'identity_alignment': identity,
      'similarity_alignment': identity,
      'gaps_source_alignment': 0.0,
      'gaps_parent_alignment': all_gaps,
      'all_gaps': all_gaps,
      'source_alignment': source_alignment,
      'parent_alignment': parent_alignment,
      'parent_alignment_modified': parent_alignment_modified,
    }

  def _make_tier2_linear_mappings(self, tier2_linear):
    if tier2_linear.height == 0:
      return self._empty_mapping_df()

    pairs = tier2_linear.select('Source Accession', 'Assigned Protein ID', 'Assigned Protein Sequence').unique(
      subset=['Source Accession', 'Assigned Protein ID']
    )
    total_pairs = pairs.height
    print(f'  Aligning {total_pairs} unique source->parent pairs...')

    aligned_count = 0
    for i, pair_row in enumerate(pairs.iter_rows(named=True)):
      result = self._get_source_parent_alignment(
        pair_row['Source Accession'], pair_row['Assigned Protein ID'], pair_row['Assigned Protein Sequence']
      )
      if result is not None:
        aligned_count += 1
      if (i + 1) % 1000 == 0 or (i + 1) == total_pairs:
        print(f'    Aligned {i + 1}/{total_pairs} pairs ({aligned_count} successful)')

    results = []
    for row in tier2_linear.iter_rows(named=True):
      alignment_data = self.alignment_cache.get((row['Source Accession'], row['Assigned Protein ID']))
      if alignment_data is None:
        continue
      mapping = self._map_linear_via_source_alignment(row, alignment_data)
      if mapping is None:
        continue
      results.append({
        'epitope_id': int(str(row['Epitope ID']).replace('.0', '')),
        'epitope_seq': row['Epitope Sequence'],
        'epitope_start': int(row['Source Starting Position']),
        'epitope_end': int(row['Source Ending Position']),
        'source_accession': row['Source Accession'],
        'parent_accession': row['Assigned Protein ID'],
        **mapping,
      })

    if not results:
      return self._empty_mapping_df()

    return pl.DataFrame(results).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })

  def _make_tier2_discontinuous_mappings(self, tier2_disc):
    if tier2_disc.height == 0:
      return self._empty_mapping_df()

    pairs = tier2_disc.select('Source Accession', 'Assigned Protein ID', 'Assigned Protein Sequence').unique(
      subset=['Source Accession', 'Assigned Protein ID']
    )
    total_pairs = pairs.height
    print(f'  Aligning {total_pairs} unique source->parent pairs for discontinuous...')

    for pair_row in pairs.iter_rows(named=True):
      self._get_source_parent_alignment(
        pair_row['Source Accession'], pair_row['Assigned Protein ID'], pair_row['Assigned Protein Sequence']
      )

    results = []
    for row in tier2_disc.iter_rows(named=True):
      alignment_data = self.alignment_cache.get((row['Source Accession'], row['Assigned Protein ID']))
      mapping = self._map_discontinuous_via_alignment(row, alignment_data)
      if mapping is None:
        continue
      results.append({
        'epitope_id': int(str(row['Epitope ID']).replace('.0', '')),
        'epitope_seq': row['Epitope Sequence'],
        'epitope_start': int(row['Source Starting Position']) if row['Source Starting Position'] is not None else None,
        'epitope_end': int(row['Source Ending Position']) if row['Source Ending Position'] is not None else None,
        'source_accession': row['Source Accession'],
        'parent_accession': row['Assigned Protein ID'],
        **mapping,
      })

    if not results:
      return self._empty_mapping_df()

    return pl.DataFrame(results).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })

  def _map_discontinuous_via_alignment(self, row, alignment_data):
    seq = row['Epitope Sequence']
    parent_protein_seq = row['Assigned Protein Sequence']

    seq = seq.replace(',  ', ', ').replace(' ,', ', ').replace(', ', ',').replace(' ', ',').replace(',,', ',')
    seq = seq[:-1] if seq.endswith(',') else seq
    split_seq = seq.split(',')
    seq_len = len(split_seq)

    if alignment_data is not None:
      pos_map = alignment_data['pos_map']
      return self._map_discontinuous_with_pos_map(split_seq, seq_len, pos_map, parent_protein_seq)
    elif parent_protein_seq:
      return self._map_discontinuous_direct(split_seq, seq_len, parent_protein_seq)
    else:
      return None

  def _map_discontinuous_with_pos_map(self, split_seq, seq_len, pos_map, parent_protein_seq):
    residue_matches = 0
    parent_seq_parts = []
    source_residues = []
    parent_residues = []

    for residue in split_seq:
      try:
        aa = residue[0]
        src_pos = int(residue[1:])
      except (ValueError, IndexError):
        continue

      source_residues.append(aa)
      mapped = pos_map.get(src_pos)

      if mapped is not None:
        parent_pos, parent_aa = mapped
        parent_seq_parts.append(f'{parent_aa}{parent_pos}')
        parent_residues.append(parent_aa)
        if aa == parent_aa:
          residue_matches += 1
      else:
        try:
          parent_aa = parent_protein_seq[src_pos - 1] if parent_protein_seq else 'X'
        except IndexError:
          parent_aa = 'X'
        parent_seq_parts.append(f'{parent_aa}{src_pos}')
        parent_residues.append(parent_aa)

    if not parent_seq_parts:
      return None

    parent_seq_str = ', '.join(parent_seq_parts)
    perc_identity = residue_matches / seq_len if seq_len > 0 else 0.0
    source_alignment = ''.join(source_residues)
    parent_alignment = ''.join(parent_residues)
    parent_alignment_modified = self._modify_parent_alignment(source_alignment, parent_alignment)

    positions = []
    for part in parent_seq_parts:
      try:
        positions.append(int(part[1:]))
      except ValueError:
        continue

    if not positions:
      return None

    return {
      'parent_seq': parent_seq_str,
      'parent_start': min(positions),
      'parent_end': max(positions),
      'identity_alignment': perc_identity,
      'similarity_alignment': perc_identity,
      'gaps_source_alignment': 0.0,
      'gaps_parent_alignment': 0.0,
      'all_gaps': 0.0,
      'source_alignment': source_alignment,
      'parent_alignment': parent_alignment,
      'parent_alignment_modified': parent_alignment_modified,
    }

  def _map_discontinuous_direct(self, split_seq, seq_len, parent_protein_seq):
    residue_matches = 0
    parent_seq_parts = []
    source_residues = []
    parent_residues = []

    for residue in split_seq:
      try:
        aa = residue[0]
        pos = int(residue[1:]) - 1
      except (ValueError, IndexError):
        continue

      source_residues.append(aa)
      try:
        parent_aa = parent_protein_seq[pos]
      except IndexError:
        parent_aa = 'X'

      parent_seq_parts.append(f'{parent_aa}{pos + 1}')
      parent_residues.append(parent_aa)
      if aa == parent_aa:
        residue_matches += 1

    if not parent_seq_parts:
      return None

    parent_seq_str = ', '.join(parent_seq_parts)
    perc_identity = residue_matches / seq_len if seq_len > 0 else 0.0
    source_alignment = ''.join(source_residues)
    parent_alignment = ''.join(parent_residues)
    parent_alignment_modified = self._modify_parent_alignment(source_alignment, parent_alignment)

    positions = []
    for part in parent_seq_parts:
      try:
        positions.append(int(part[1:]))
      except ValueError:
        continue

    if not positions:
      return None

    return {
      'parent_seq': parent_seq_str,
      'parent_start': min(positions),
      'parent_end': max(positions),
      'identity_alignment': perc_identity,
      'similarity_alignment': perc_identity,
      'gaps_source_alignment': 0.0,
      'gaps_parent_alignment': 0.0,
      'all_gaps': 0.0,
      'source_alignment': source_alignment,
      'parent_alignment': parent_alignment,
      'parent_alignment_modified': parent_alignment_modified,
    }

  def _make_tier3_mappings(self, tier3):
    if tier3.height == 0:
      return self._empty_mapping_df()

    local_aligner = PairwiseAligner()
    local_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    local_aligner.mode = 'local'
    local_aligner.open_gap_score = -11
    local_aligner.extend_gap_score = -1

    total = tier3.height
    results = []
    failed = 0
    for i, row in enumerate(tier3.iter_rows(named=True)):
      epitope_seq = row['Epitope Sequence']
      parent_seq = row['Assigned Protein Sequence']

      if not epitope_seq or not parent_seq:
        failed += 1
        continue

      try:
        alignments = local_aligner.align(epitope_seq, parent_seq)
        if not alignments:
          failed += 1
          continue

        top = alignments[0]
        qseq = str(top[0])
        sseq = str(top[1])

        qseq_clean = ''
        sseq_clean = ''
        for qc, sc in zip(qseq, sseq):
          if qc != '-' or sc != '-':
            qseq_clean += qc
            sseq_clean += sc

        if not qseq_clean or qseq_clean == '-' * len(qseq_clean):
          failed += 1
          continue

        q_aligned = top.aligned
        if len(q_aligned[0]) == 0 or len(q_aligned[1]) == 0:
          failed += 1
          continue

        parent_start = int(q_aligned[1][0][0]) + 1
        parent_end = int(q_aligned[1][-1][-1])

        alignment_len = len(qseq_clean)
        matches = sum(1 for q, s in zip(qseq_clean, sseq_clean) if q == s and q != '-')
        identity = matches / alignment_len if alignment_len > 0 else 0.0
        gaps_source = sum(1 for c in qseq_clean if c == '-') / alignment_len if alignment_len > 0 else 0.0
        gaps_parent = sum(1 for c in sseq_clean if c == '-') / alignment_len if alignment_len > 0 else 0.0
        total_gaps = (sum(1 for c in qseq_clean if c == '-') + sum(1 for c in sseq_clean if c == '-')) / alignment_len if alignment_len > 0 else 0.0

        parent_alignment_modified = self._modify_parent_alignment(qseq_clean, sseq_clean)

        results.append({
          'epitope_id': int(str(row['Epitope ID']).replace('.0', '')),
          'epitope_seq': epitope_seq,
          'epitope_start': int(row['Source Starting Position']) if row['Source Starting Position'] is not None else None,
          'epitope_end': int(row['Source Ending Position']) if row['Source Ending Position'] is not None else None,
          'source_accession': row['Source Accession'],
          'parent_accession': row['Assigned Protein ID'],
          'parent_seq': sseq_clean,
          'parent_start': parent_start,
          'parent_end': parent_end,
          'identity_alignment': identity,
          'similarity_alignment': identity,
          'gaps_source_alignment': gaps_source,
          'gaps_parent_alignment': gaps_parent,
          'all_gaps': total_gaps,
          'source_alignment': qseq_clean,
          'parent_alignment': sseq_clean,
          'parent_alignment_modified': parent_alignment_modified,
        })
      except Exception as e:
        print(f'  WARNING: Tier 3 alignment failed for epitope {row["Epitope ID"]}: {e}')
        failed += 1
        continue

      if (i + 1) % 1000 == 0 or (i + 1) == total:
        print(f'    Tier 3 progress: {i + 1}/{total} ({failed} failed)')

    if not results:
      return self._empty_mapping_df()

    return pl.DataFrame(results).cast({
      'epitope_id': pl.Int64, 'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'parent_start': pl.Int64, 'parent_end': pl.Int64
    })

  def _modify_parent_alignment(self, qseq, sseq):
    if len(qseq) != len(sseq):
      return ''
    result = ''
    for i in range(len(qseq)):
      if qseq[i] != sseq[i]:
        result += sseq[i].lower()
      else:
        result += sseq[i]
    return result

  def _empty_mapping_df(self):
    return pl.DataFrame(schema={
      'epitope_id': pl.Int64, 'epitope_seq': pl.String,
      'epitope_start': pl.Int64, 'epitope_end': pl.Int64,
      'source_accession': pl.String, 'parent_accession': pl.String,
      'parent_seq': pl.String, 'parent_start': pl.Int64, 'parent_end': pl.Int64,
      'identity_alignment': pl.Float64, 'similarity_alignment': pl.Float64,
      'gaps_source_alignment': pl.Float64, 'gaps_parent_alignment': pl.Float64,
      'all_gaps': pl.Float64, 'source_alignment': pl.String,
      'parent_alignment': pl.String, 'parent_alignment_modified': pl.String,
    })

  def _combine_and_write(self, *mapping_dfs):
    combined = pl.concat([df for df in mapping_dfs if df.height > 0])
    combined = combined.filter(pl.col('epitope_id').is_not_null())

    col_order = [
      'epitope_id', 'epitope_seq', 'epitope_start', 'epitope_end', 'source_accession',
      'parent_accession', 'parent_seq', 'parent_start', 'parent_end', 'identity_alignment',
      'similarity_alignment', 'gaps_source_alignment', 'gaps_parent_alignment', 'all_gaps',
      'source_alignment', 'parent_alignment', 'parent_alignment_modified'
    ]
    combined = combined.select(col_order)

    print(f'Total combined mappings: {combined.height}')
    combined.write_csv(build_path / 'arborist' / 'epitope-mappings.tsv', separator='\t')
    print(f'Written to {build_path / "arborist" / "epitope-mappings.tsv"}')


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-b', '--build_path', type=str,
    default=Path(__file__).parents[3] / 'build'
  )
  args = parser.parse_args()

  build_path = Path(args.build_path)

  numeric_cols = [
    'Epitope ID', 'Species Taxon ID', 'Organism ID',
    'Source Starting Position', 'Source Ending Position',
    'Assigned Protein Starting Position', 'Assigned Protein Ending Position'
  ]

  print(f'Reading assignments from {build_path / "arborist" / "all-peptide-assignments.tsv"}')
  assignments = pl.read_csv(
    build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t',
    schema_overrides={col: pl.Float64 for col in numeric_cols}
  ).with_columns(pl.col(numeric_cols).cast(pl.Int64))

  print(f'Reading source data from {build_path / "arborist" / "all-source-data.tsv"}')
  source_data = pl.read_csv(build_path / 'arborist' / 'all-source-data.tsv', separator='\t')

  print(f'Reading species data from {build_path / "arborist" / "all-species-data.tsv"}')
  species_data = pl.read_csv(build_path / 'arborist' / 'all-species-data.tsv', separator='\t')

  print('Generating parent data...')
  all_parent_data = get_all_parent_data(assignments, source_data, species_data)
  make_source_parents(all_parent_data)
  print('Wrote source-parents.tsv')
  make_parent_proteins(all_parent_data)
  print('Wrote parent-proteins.tsv')

  print('Starting epitope mappings...')
  mapper = EpitopeMapper(assignments, source_data)
  mapper.make_epitope_mappings()
  print('Done.')
