import datetime

import polars as pl
import pytest

from protein_tree.stats import (
  compute_week_stats,
  upsert_stats,
  iso_week,
  METRICS,
  ID_COLUMNS,
)

ASSIGNMENT_COLUMNS = {
  'Species Taxon ID': pl.String,
  'Species Name': pl.String,
  'Source Accession': pl.String,
  'Source Assigned Protein ID': pl.String,
  'ARC Assignment': pl.String,
  'Epitope ID': pl.String,
  'Assigned Protein ID': pl.String,
  'Assigned Protein Starting Position': pl.String,
}


def _assignments(rows):
  return pl.DataFrame(rows, schema=ASSIGNMENT_COLUMNS, orient='row')


def _val(long, species_id, metric):
  return long.filter(
    (pl.col('Species ID') == species_id) & (pl.col('Metric') == metric)
  ).row(0, named=True)['value']


def test_metrics_two_levels_and_arc():
  # rows: taxon, name, source, source_pid, arc, epitope, assigned_pid, position
  rows = [
    # source S1 assigned; epitope E1 located exactly
    ['9606', 'Homo sapiens', 'S1', 'P1', None, 'E1', 'P1', '5'],
    # same source; epitope E2 has parent (fallback) but NO position -> assigned, not located
    ['9606', 'Homo sapiens', 'S1', 'P1', None, 'E2', 'P1', None],
    # source S2 NOT assigned, no ARC; epitope E3 unassigned -> lost
    ['9606', 'Homo sapiens', 'S2', None, None, 'E3', None, None],
    # source S3 assigned only via ARC; epitope E4 has no protein id -> assigned via ARC
    ['9606', 'Homo sapiens', 'S3', None, 'TCR_alpha', 'E4', None, None],
  ]
  long = compute_week_stats(_assignments(rows))

  assert _val(long, '9606', 'source_antigens_total') == 3        # S1,S2,S3
  assert _val(long, '9606', 'source_antigens_assigned') == 2     # S1 (protein), S3 (ARC)
  assert _val(long, '9606', 'epitopes_total') == 4               # E1..E4
  assert _val(long, '9606', 'epitopes_assigned') == 3            # E1,E2 (pid), E4 (ARC); E3 lost
  assert _val(long, '9606', 'epitopes_located_exact') == 1       # only E1


def test_empty_string_counts_as_unassigned():
  rows = [
    ['9606', 'Homo sapiens', 'S1', '', '', 'E1', '', ''],   # empty strings, not null
  ]
  long = compute_week_stats(_assignments(rows))
  assert _val(long, '9606', 'source_antigens_assigned') == 0
  assert _val(long, '9606', 'epitopes_assigned') == 0
  assert _val(long, '9606', 'epitopes_located_exact') == 0


def test_total_is_column_sum():
  rows = [
    ['9606', 'Homo sapiens', 'S1', 'P1', None, 'E1', 'P1', '5'],
    ['111', 'Spirosoma', 'S9', 'P9', None, 'E9', 'P9', '2'],
  ]
  long = compute_week_stats(_assignments(rows))
  for metric in METRICS:
    total = _val(long, 'TOTAL', metric)
    per_species = sum(
      long.filter((pl.col('Metric') == metric) & (pl.col('Species ID') != 'TOTAL'))['value']
    )
    assert total == per_species


def test_epitopes_total_is_unique_epitope_id():
  # same Epitope ID under two sources -> counted once
  rows = [
    ['9606', 'Homo sapiens', 'S1', 'P1', None, 'E1', 'P1', '5'],
    ['9606', 'Homo sapiens', 'S2', 'P2', None, 'E1', 'P2', '7'],
  ]
  long = compute_week_stats(_assignments(rows))
  assert _val(long, '9606', 'epitopes_total') == 1
  assert _val(long, '9606', 'source_antigens_total') == 2


# ---- upsert / master sheet ----

def _week_long(species_id, values):
  """Build a minimal week-long frame for one species across all metrics."""
  return pl.DataFrame(
    {
      'Species ID': [species_id] * len(METRICS),
      'Species Name': ['X'] * len(METRICS),
      'Metric': METRICS,
      'value': [values[m] for m in METRICS],
    },
    schema={'Species ID': pl.String, 'Species Name': pl.String, 'Metric': pl.String, 'value': pl.Int64},
  )


def _cell(sheet, species_id, metric, week):
  return sheet.filter(
    (pl.col('Species ID') == species_id) & (pl.col('Metric') == metric)
  ).row(0, named=True)[week]


def test_upsert_fresh_creates_week_column():
  wk = _week_long('9606', {m: i for i, m in enumerate(METRICS)})
  sheet = upsert_stats(None, '2026-W01', wk)
  assert '2026-W01' in sheet.columns
  assert _cell(sheet, '9606', 'epitopes_total', '2026-W01') == METRICS.index('epitopes_total')


def test_upsert_new_week_appends_and_keeps_old():
  s1 = upsert_stats(None, '2026-W01', _week_long('9606', {m: 10 for m in METRICS}))
  s2 = upsert_stats(s1, '2026-W02', _week_long('9606', {m: 20 for m in METRICS}))
  assert '2026-W01' in s2.columns and '2026-W02' in s2.columns
  assert _cell(s2, '9606', 'epitopes_total', '2026-W01') == 10
  assert _cell(s2, '9606', 'epitopes_total', '2026-W02') == 20


def test_upsert_same_week_overwrites_no_duplicate():
  s1 = upsert_stats(None, '2026-W02', _week_long('9606', {m: 10 for m in METRICS}))
  s2 = upsert_stats(s1, '2026-W02', _week_long('9606', {m: 99 for m in METRICS}))
  assert s2.columns.count('2026-W02') == 1
  assert _cell(s2, '9606', 'epitopes_total', '2026-W02') == 99


def test_upsert_species_churn_outer_pivot():
  s1 = upsert_stats(None, '2026-W01', _week_long('108', {m: 1 for m in METRICS}))
  # week 2 introduces a different species; 108 absent this week
  s2 = upsert_stats(s1, '2026-W02', _week_long('222', {m: 2 for m in METRICS}))
  species = set(s2['Species ID'])
  assert {'108', '222'}.issubset(species)
  # 108 has no W02 value -> null; 222 has no W01 value -> null
  assert _cell(s2, '108', 'epitopes_total', '2026-W02') is None
  assert _cell(s2, '222', 'epitopes_total', '2026-W01') is None


def test_iso_week_format():
  assert iso_week(datetime.date(2026, 1, 5)) == '2026-W02'
  assert iso_week(datetime.date(2026, 7, 13)) == '2026-W29'
