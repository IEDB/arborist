import polars as pl
import pytest

from protein_tree import immunomebrowser

PARENT = 'MKTAYIAKQRQISFVKSHFSRQ'
SOURCE = 'MKTAYIAKQQISFVKSHFSRQ'  # PARENT without the R at parent position 10


def test_position_map_follows_both_gap_directions():
  # source has no residue opposite parent 'G' (gap in source), and its 'E'
  # sits opposite a gap in the parent
  position_map = immunomebrowser.PositionMap('AC-DE', 'ACGD-')
  assert position_map.get(1) == (1, 'A')
  assert position_map.get(2) == (2, 'C')
  assert position_map.get(3) == (4, 'D')   # parent position skips the source gap
  assert position_map.get(4) is None       # aligned to a gap in the parent
  assert position_map.get(0) is None
  assert position_map.get(5) is None       # beyond the source sequence


def _assignments():
  columns = {
    'Epitope ID': [1, 2, 3, 4, 5],
    'Epitope Sequence': ['KQRQ', 'QISF', 'QISF', 'F5, K8, H17', 'AYIA'],
    'Source Accession': ['S1', 'S1', 'S1', 'S1', 'NOSEQ'],
    'Source Starting Position': [8, 10, 10, None, 4],
    'Source Ending Position': [11, 13, 13, None, 7],
    'Assigned Protein ID': ['P1', 'P1', 'P1', 'P1', 'P1'],
    'Assigned Protein Starting Position': [8, None, None, None, None],
    'Assigned Protein Ending Position': [11, None, None, None, None],
  }
  return pl.DataFrame(columns, schema_overrides={
    'Source Starting Position': pl.Int64, 'Source Ending Position': pl.Int64,
    'Assigned Protein Starting Position': pl.Int64, 'Assigned Protein Ending Position': pl.Int64,
  })


def _run_mapper(tmp_path, monkeypatch, chunk_rows):
  build = tmp_path / f'build_{chunk_rows}'
  (build / 'arborist').mkdir(parents=True)
  monkeypatch.setattr(immunomebrowser, 'build_path', build, raising=False)
  monkeypatch.setattr(immunomebrowser, 'CHUNK_ROWS', chunk_rows)
  source_data = pl.DataFrame({'Source Accession': ['S1', 'NOSEQ'], 'Sequence': [SOURCE, None]},
                             schema_overrides={'Sequence': pl.String})
  protein_sequences = pl.DataFrame({'Assigned Protein ID': ['P1'], 'Assigned Protein Sequence': [PARENT]})
  mapper = immunomebrowser.EpitopeMapper(_assignments(), source_data, protein_sequences)
  mapper.make_epitope_mappings()
  return pl.read_csv(build / 'arborist' / 'epitope-mappings.tsv', separator='\t')


def test_tiers_map_through_the_source_alignment(tmp_path, monkeypatch):
  out = _run_mapper(tmp_path, monkeypatch, chunk_rows=50_000).sort('epitope_id')
  by_id = {row['epitope_id']: row for row in out.iter_rows(named=True)}

  assert by_id[1]['parent_start'] == 8 and by_id[1]['identity_alignment'] == 1.0  # exact
  # source positions 10-13 are parent positions 11-14 because of the deleted R
  assert (by_id[2]['parent_start'], by_id[2]['parent_end']) == (11, 14)
  assert by_id[2]['identity_alignment'] == 1.0
  assert by_id[3]['parent_start'] == 11            # same pair, served from the cache
  assert by_id[4]['parent_accession'] == 'P1'      # discontinuous, via alignment
  assert by_id[4]['parent_start'] < by_id[4]['parent_end']
  assert (by_id[5]['parent_start'], by_id[5]['parent_end']) == (4, 7)  # tier 3: no source sequence


def test_chunking_does_not_change_rows_or_order(tmp_path, monkeypatch):
  one_row_chunks = _run_mapper(tmp_path, monkeypatch, chunk_rows=1)
  one_chunk = _run_mapper(tmp_path, monkeypatch, chunk_rows=50_000)
  assert one_row_chunks.equals(one_chunk)


def test_protein_sequences_collapse_rows_and_refuse_conflicts(tmp_path):
  header = 'Assigned Protein ID\tAssigned Protein Sequence\n'
  good = tmp_path / 'good.tsv'
  good.write_text(header + 'P1\tAAAA\nP1\tAAAA\nP2\tCCCC\n\t\n')
  collapsed = immunomebrowser.read_protein_sequences(good)
  assert sorted(collapsed['Assigned Protein ID'].to_list()) == ['P1', 'P2']

  conflicting = tmp_path / 'conflicting.tsv'
  conflicting.write_text(header + 'P1\tAAAA\nP1\tCCCC\n')
  with pytest.raises(ValueError, match='more than one'):
    immunomebrowser.read_protein_sequences(conflicting)
