"""The promotion gate: exit codes must reflect what actually changed."""

import polars as pl
import pytest

from protein_tree.alert import Thresholds
from protein_tree.diff_metrics import (
  EXIT_CHANGED,
  EXIT_ESCALATED,
  EXIT_IDENTICAL,
  compare_metrics,
  compare_selection,
  delta_frame,
  latest_archive,
  main,
  verdict,
)

ID_COLUMNS = ['Species ID', 'Species Name', 'Metric']
REPORT_COLUMNS = ['Species ID', 'Species Name', 'Status', 'Proteome ID']


def sheet(rows, weeks=('2026-W29', '2026-W30')):
  """rows: (species_id, name, metric, baseline, current)."""
  return pl.DataFrame(
    [{'Species ID': r[0], 'Species Name': r[1], 'Metric': r[2],
      weeks[0]: str(r[3]), weeks[1]: str(r[4])} for r in rows],
    schema={c: pl.String for c in list(ID_COLUMNS) + list(weeks)},
  )


def write_report(path, rows):
  path.parent.mkdir(parents=True, exist_ok=True)
  pl.DataFrame(
    [dict(zip(REPORT_COLUMNS, row)) for row in rows],
    schema={c: pl.String for c in REPORT_COLUMNS},
  ).write_csv(path, separator='\t')
  return path


@pytest.fixture
def thresholds():
  return Thresholds()


def test_identical_numbers_are_promotable(thresholds):
  deltas, weeks = compare_metrics(sheet([
    ('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 5000),
    ('TOTAL', '', 'epitopes_assigned', 5000, 5000),
  ]), thresholds)

  assert deltas == []
  assert weeks == ['2026-W29', '2026-W30']
  assert verdict(deltas, []) == EXIT_IDENTICAL


def test_small_drop_holds_but_does_not_escalate(thresholds):
  deltas, _ = compare_metrics(sheet([
    ('1234', 'Small species', 'epitopes_assigned', 50, 48),
  ]), thresholds)

  assert len(deltas) == 1
  assert deltas[0]['Delta'] == -2
  assert deltas[0]['Severity'] == 'INFO'
  assert verdict(deltas, []) == EXIT_CHANGED


def test_the_2026_02_collapse_escalates(thresholds):
  """T. cruzi lost ~9% of assigned rows while the rate held; that must escalate."""
  deltas, _ = compare_metrics(sheet([
    ('5693', 'Trypanosoma cruzi', 'epitopes_assigned', 22000, 20000),
  ]), thresholds)

  assert deltas[0]['Severity'] == 'CRITICAL'
  assert verdict(deltas, []) == EXIT_ESCALATED


def test_total_drop_uses_the_aggregate_threshold(thresholds):
  """5% on TOTAL is CRITICAL even though the same ratio per-species would not be."""
  deltas, _ = compare_metrics(sheet([
    ('TOTAL', '', 'epitopes_assigned', 100000, 94000),
  ]), thresholds)

  assert deltas[0]['Severity'] == 'CRITICAL'
  assert verdict(deltas, []) == EXIT_ESCALATED


def test_a_rise_is_reported_but_never_escalates(thresholds):
  deltas, _ = compare_metrics(sheet([
    ('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 9000),
  ]), thresholds)

  assert deltas[0]['Delta'] == 4000
  assert deltas[0]['Severity'] == 'INFO'
  assert verdict(deltas, []) == EXIT_CHANGED


def test_single_week_sheet_has_no_baseline(thresholds):
  one_week = pl.DataFrame(
    [{'Species ID': '9606', 'Species Name': 'Homo sapiens',
      'Metric': 'epitopes_assigned', '2026-W30': '5000'}],
    schema={c: pl.String for c in ID_COLUMNS + ['2026-W30']},
  )
  deltas, weeks = compare_metrics(one_week, thresholds)

  assert deltas == []
  assert len(weeks) == 1


def test_selected_to_empty_escalates(tmp_path):
  baseline = write_report(tmp_path / 'reports' / 'proteome-selection-report-2026-06-14.tsv',
                          [('9606', 'Homo sapiens', 'SELECTED', 'UP000005640')])
  current = write_report(tmp_path / 'arborist' / 'proteome-selection-report.tsv',
                         [('9606', 'Homo sapiens', 'EMPTY', '')])

  changes = compare_selection(current, baseline)

  status = next(c for c in changes if c['Change'] == 'status')
  assert status['Severity'] == 'CRITICAL'
  assert verdict([], changes) == EXIT_ESCALATED


def test_recovery_to_selected_does_not_escalate(tmp_path):
  """EMPTY -> SELECTED is good news; it should not block a promotion."""
  baseline = write_report(tmp_path / 'reports' / 'proteome-selection-report-2026-06-14.tsv',
                          [('9606', 'Homo sapiens', 'EMPTY', '')])
  current = write_report(tmp_path / 'arborist' / 'proteome-selection-report.tsv',
                         [('9606', 'Homo sapiens', 'SELECTED', 'UP000005640')])

  changes = compare_selection(current, baseline)

  assert [c['Change'] for c in changes] == ['status', 'proteome_id']
  assert next(c for c in changes if c['Change'] == 'status')['Severity'] == 'INFO'
  assert verdict([], changes) == EXIT_CHANGED


def test_proteome_id_flip_is_noted_not_escalated(tmp_path):
  """New UPIDs are normal churn. If one cost us anything, the counts escalate."""
  baseline = write_report(tmp_path / 'reports' / 'proteome-selection-report-2026-06-14.tsv',
                          [('5693', 'T. cruzi', 'SELECTED', 'UP000002296')])
  current = write_report(tmp_path / 'arborist' / 'proteome-selection-report.tsv',
                         [('5693', 'T. cruzi', 'SELECTED', 'UP000694941')])

  changes = compare_selection(current, baseline)

  assert len(changes) == 1
  assert changes[0]['Change'] == 'proteome_id'
  assert changes[0]['Severity'] == 'INFO'
  assert verdict([], changes) == EXIT_CHANGED


def test_new_species_is_not_a_regression(tmp_path):
  """A species that did not exist last build has nothing to regress from."""
  baseline = write_report(tmp_path / 'reports' / 'proteome-selection-report-2026-06-14.tsv',
                          [('9606', 'Homo sapiens', 'SELECTED', 'UP000005640')])
  current = write_report(tmp_path / 'arborist' / 'proteome-selection-report.tsv',
                         [('9606', 'Homo sapiens', 'SELECTED', 'UP000005640'),
                          ('1234', 'New species', 'EMPTY', '')])

  assert compare_selection(current, baseline) == []


def test_missing_baseline_report_is_not_a_change(tmp_path):
  current = write_report(tmp_path / 'arborist' / 'proteome-selection-report.tsv',
                         [('9606', 'Homo sapiens', 'SELECTED', 'UP000005640')])
  assert compare_selection(current, None) == []


def test_latest_archive_picks_the_newest(tmp_path):
  reports = tmp_path / 'reports'
  for day in ('2026-05-01', '2026-07-26', '2026-06-14'):
    write_report(reports / f'proteome-selection-report-{day}.tsv', [])

  assert latest_archive(reports).name == 'proteome-selection-report-2026-07-26.tsv'
  assert latest_archive(tmp_path / 'nope') is None


def test_delta_frame_is_typed_when_empty():
  frame = delta_frame([])
  assert frame.height == 0
  assert frame.columns == ['Species ID', 'Species Name', 'Metric',
                           'Baseline', 'Current', 'Delta', 'Severity']


def test_main_exits_identical_and_writes_an_artifact(tmp_path):
  stats = tmp_path / 'stats.tsv'
  sheet([
    ('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 5000),
  ]).write_csv(stats, separator='\t')
  build = tmp_path / 'build'
  output = build / 'reports' / 'diff-metrics.tsv'

  code = main(['-b', str(build), '-s', str(stats)])

  assert code == EXIT_IDENTICAL
  assert output.exists()
  assert pl.read_csv(output, separator='\t').height == 0


def test_main_exits_escalated_on_a_collapse(tmp_path):
  stats = tmp_path / 'stats.tsv'
  sheet([
    ('5693', 'Trypanosoma cruzi', 'epitopes_assigned', 22000, 20000),
  ]).write_csv(stats, separator='\t')
  build = tmp_path / 'build'

  code = main(['-b', str(build), '-s', str(stats)])

  assert code == EXIT_ESCALATED
  written = pl.read_csv(build / 'reports' / 'diff-metrics.tsv', separator='\t')
  assert written['Severity'].to_list() == ['CRITICAL']
  assert written['Delta'].to_list() == [-2000]


def test_main_without_a_stats_sheet_fails_closed(tmp_path):
  """No numbers means no evidence; that must never read as promotable."""
  assert main(['-b', str(tmp_path / 'build'), '-s', str(tmp_path / 'missing.tsv')]) == EXIT_ESCALATED


def test_upid_flip_that_actually_cost_us_still_escalates(thresholds):
  """The UPID change is INFO, but the assignment collapse it caused is not."""
  deltas, _ = compare_metrics(sheet([
    ('5693', 'Trypanosoma cruzi', 'epitopes_assigned', 22000, 20000),
  ]), thresholds)
  changes = [{'Species ID': '5693', 'Change': 'proteome_id', 'Severity': 'INFO'}]

  assert verdict(deltas, changes) == EXIT_ESCALATED
