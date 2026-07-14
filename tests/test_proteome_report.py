import polars as pl
import pytest

from protein_tree.report_proteomes import (
  build_report,
  classify_species,
  summarize,
  STATUS_SELECTED,
  STATUS_ORPHANS,
  STATUS_EMPTY,
  STATUS_MISSING,
)

ACTIVE_COLUMNS = [
  'Species Key', 'Species ID', 'Species Label', 'Active Taxa', 'Group', 'Epitope Count'
]
SPECIES_DATA_COLUMNS = [
  'Species Taxon ID', 'Species Name', 'Proteome ID', 'Proteome Taxon ID',
  'Proteome Type', 'Proteome Label'
]


@pytest.fixture
def build_path(tmp_path):
  return tmp_path / 'build'


def _write_active_species(build_path, rows):
  path = build_path / 'arborist' / 'active-species.tsv'
  path.parent.mkdir(parents=True, exist_ok=True)
  pl.DataFrame(rows, schema=ACTIVE_COLUMNS, orient='row').write_csv(path, separator='\t')


def _species_dir(build_path, taxon_id):
  d = build_path / 'species' / str(taxon_id)
  d.mkdir(parents=True, exist_ok=True)
  return d


def _write_species_data(species_dir, taxon_id, name, proteome_id, proteome_type, label=''):
  row = [[str(taxon_id), name, proteome_id, '', proteome_type, label]]
  pl.DataFrame(row, schema=SPECIES_DATA_COLUMNS, orient='row').write_csv(
    species_dir / 'species-data.tsv', separator='\t'
  )


def _write_proteome_list(species_dir, n):
  schema = {'Proteome ID': pl.String, 'Proteome Type': pl.String}
  rows = [{'Proteome ID': f'UP{i:09d}', 'Proteome Type': 'Reference'} for i in range(n)]
  pl.DataFrame(rows, schema=schema).write_csv(species_dir / 'proteome-list.tsv', separator='\t')


def _row(report, species_id):
  return report.filter(pl.col('Species ID') == species_id).row(0, named=True)


def test_selected(build_path):
  _write_active_species(build_path, [['9606-h', '9606', 'Homo sapiens', '9606', 'vertebrate', '100']])
  d = _species_dir(build_path, 9606)
  (d / 'proteome.fasta').write_text('>sp|P12345|X_HUMAN\nMKVLA\n')
  _write_species_data(d, 9606, 'Homo sapiens', 'UP000005640',
                      'Reference and representative proteome', 'Homo sapiens')
  _write_proteome_list(d, 3)

  row = _row(build_report(build_path), '9606')
  assert row['Status'] == STATUS_SELECTED
  assert row['Proteome ID'] == 'UP000005640'
  assert row['Group'] == 'vertebrate'
  assert row['Species Name'] == 'Homo sapiens'
  assert row['FASTA Bytes'] > 0
  assert row['Candidate Count'] == 3
  assert row['Reason'] == 'ok'


def test_orphans(build_path):
  _write_active_species(build_path, [['x', '12345', 'Gaunavirus GA1', '12345', 'virus', '5']])
  d = _species_dir(build_path, 12345)
  (d / 'proteome.fasta').write_text('>orphan1\nMKV\n')
  _write_species_data(d, 12345, 'Gaunavirus GA1', '', 'Orphans', 'Gaunavirus GA1')

  row = _row(build_report(build_path), '12345')
  assert row['Status'] == STATUS_ORPHANS
  assert row['FASTA Bytes'] > 0


def test_empty_no_candidates(build_path):
  _write_active_species(build_path, [['x', '111', 'Empty species', '111', 'bacterium', '2']])
  d = _species_dir(build_path, 111)
  (d / 'proteome.fasta').write_text('')  # 0 bytes
  _write_proteome_list(d, 0)             # header only -> 0 candidates

  row = _row(build_report(build_path), '111')
  assert row['Status'] == STATUS_EMPTY
  assert row['FASTA Bytes'] == 0
  assert row['Candidate Count'] == 0
  assert row['Reason'] == 'no-candidates'


def test_empty_fetch_empty(build_path):
  _write_active_species(build_path, [['x', '222', 'Fetch failed', '222', 'bacterium', '2']])
  d = _species_dir(build_path, 222)
  (d / 'proteome.fasta').write_text('')  # 0 bytes despite candidates existing
  _write_proteome_list(d, 5)

  row = _row(build_report(build_path), '222')
  assert row['Status'] == STATUS_EMPTY
  assert row['Candidate Count'] == 5
  assert row['Reason'] == 'fetch-empty'


def test_missing(build_path):
  _write_active_species(build_path, [['x', '333', 'No dir species', '333', 'bacterium', '2']])
  # no species dir created

  row = _row(build_report(build_path), '333')
  assert row['Status'] == STATUS_MISSING
  assert row['FASTA Bytes'] == 0
  assert row['Reason'] == 'no species dir'


def test_selected_without_metadata(build_path):
  _write_active_species(build_path, [['x', '444', 'No metadata', '444', 'bacterium', '2']])
  d = _species_dir(build_path, 444)
  (d / 'proteome.fasta').write_text('>a\nMK\n')  # content but no species-data.tsv

  row = _row(build_report(build_path), '444')
  assert row['Status'] == STATUS_SELECTED
  assert row['Reason'] == 'no-metadata'


def test_row_count_and_empties_isolable(build_path):
  _write_active_species(build_path, [
    ['a', '9606', 'Homo sapiens', '9606', 'vertebrate', '100'],
    ['b', '111', 'Empty one', '111', 'bacterium', '2'],
    ['c', '333', 'Missing one', '333', 'bacterium', '2'],
  ])
  sel = _species_dir(build_path, 9606)
  (sel / 'proteome.fasta').write_text('>a\nMK\n')
  _write_species_data(sel, 9606, 'Homo sapiens', 'UP000005640', 'Reference proteome')
  empty = _species_dir(build_path, 111)
  (empty / 'proteome.fasta').write_text('')

  report = build_report(build_path)
  assert report.height == 3  # one row per active species

  flagged = set(report.filter(pl.col('Status').is_in([STATUS_EMPTY, STATUS_MISSING]))['Species ID'])
  assert flagged == {'111', '333'}


def test_summarize_counts(build_path):
  _write_active_species(build_path, [
    ['a', '9606', 'Homo sapiens', '9606', 'vertebrate', '100'],
    ['b', '111', 'Empty one', '111', 'bacterium', '2'],
  ])
  sel = _species_dir(build_path, 9606)
  (sel / 'proteome.fasta').write_text('>a\nMK\n')
  _write_species_data(sel, 9606, 'Homo sapiens', 'UP000005640', 'Reference proteome')
  _species_dir(build_path, 111)  # dir but no fasta -> EMPTY

  summary = summarize(build_report(build_path))
  assert '1 selected' in summary
  assert '1 EMPTY' in summary


def test_classify_species_missing_dir(tmp_path):
  result = classify_species(tmp_path / 'does-not-exist')
  assert result['Status'] == STATUS_MISSING
