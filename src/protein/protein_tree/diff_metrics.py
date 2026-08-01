#!/usr/bin/env python3
"""Compare a fresh build's numbers against the previous one, and exit nonzero.

The gate `alert.py` cannot be. alert.py grades drops but its Makefile leaf always
exits 0 by design, and it keys on assignment *rate*, which held at 99.4-99.6%
straight through the 2026_02 collapse that lost ~9% of assigned rows. A promotion
gate has to key on volume and on selection-status transitions, and it has to fail.

Two comparisons:
  1. protein-tree-weekly-stats.tsv -- the new ISO-week column vs the one before it,
     five metrics per species plus TOTAL, graded with alert.py's own thresholds.
  2. proteome-selection-report.tsv vs its most recent archive in build/reports/ --
     Status regressions (SELECTED -> ORPHANS/EMPTY/MISSING) and Proteome ID flips.

Exit codes:
   0  identical. Every metric equal, no status or UPID change.
  10  changed, within thresholds. Needs a human.
  20  escalated. A WARNING/CRITICAL drop, a status regression, or a UPID flip.

Writes a per-species delta TSV so the verdict is auditable rather than asserted.
"""

import argparse
import sys
from pathlib import Path

import polars as pl

from protein_tree.alert import (
  OVERVIEW_METRICS,
  SEVERITY_CRITICAL,
  SEVERITY_INFO,
  SEVERITY_WARNING,
  Thresholds,
  classify_severity,
  week_columns,
)
from protein_tree.stats import REPO_ROOT

EXIT_IDENTICAL = 0
EXIT_CHANGED = 10
EXIT_ESCALATED = 20

STATUS_RANK = {'SELECTED': 0, 'ORPHANS': 1, 'EMPTY': 2, 'MISSING': 3}
REPORT_NAME = 'proteome-selection-report.tsv'

DELTA_COLUMNS = [
  'Species ID', 'Species Name', 'Metric', 'Baseline', 'Current', 'Delta', 'Severity',
]


def _to_int(value) -> int:
  if value is None or value == '':
    return 0
  try:
    return int(value)
  except (ValueError, TypeError):
    return 0


def compare_metrics(sheet: pl.DataFrame, thresholds: Thresholds):
  """Per-species metric deltas between the last two ISO-week columns."""
  weeks = week_columns(sheet)
  if len(weeks) < 2:
    return [], weeks
  baseline_week, current_week = weeks[-2], weeks[-1]

  deltas = []
  for row in sheet.iter_rows(named=True):
    if row['Metric'] not in OVERVIEW_METRICS:
      continue
    baseline = _to_int(row[baseline_week])
    current = _to_int(row[current_week])
    if baseline == current:
      continue
    severity = classify_severity(
      baseline, current, thresholds, is_total=row['Species ID'] == 'TOTAL'
    )
    deltas.append({
      'Species ID': row['Species ID'],
      'Species Name': row['Species Name'],
      'Metric': row['Metric'],
      'Baseline': baseline,
      'Current': current,
      'Delta': current - baseline,
      # a rise is a change worth reporting, but never an escalation
      'Severity': severity or SEVERITY_INFO,
    })
  return deltas, [baseline_week, current_week]


def latest_archive(reports_dir: Path):
  """The most recent archived selection report, by filename date."""
  if not reports_dir.exists():
    return None
  archives = sorted(reports_dir.glob(f'{Path(REPORT_NAME).stem}-*.tsv'))
  return archives[-1] if archives else None


def compare_selection(current_report: Path, baseline_report: Path):
  """Status regressions and Proteome ID flips, keyed by species."""
  if not current_report.exists() or baseline_report is None or not baseline_report.exists():
    return []
  current = pl.read_csv(current_report, separator='\t', infer_schema_length=0)
  baseline = pl.read_csv(baseline_report, separator='\t', infer_schema_length=0)
  before = {row['Species ID']: row for row in baseline.iter_rows(named=True)}

  changes = []
  for row in current.iter_rows(named=True):
    was = before.get(row['Species ID'])
    if was is None:
      continue
    # an empty TSV cell parses as null; treat it as the empty string it was
    old_status, new_status = was.get('Status') or '', row.get('Status') or ''
    old_upid, new_upid = was.get('Proteome ID') or '', row.get('Proteome ID') or ''
    if old_status != new_status:
      # Only a move toward less proteome is gate-worthy; recovery is good news.
      regressed = STATUS_RANK.get(new_status, 9) > STATUS_RANK.get(old_status, 9)
      changes.append({
        'Species ID': row['Species ID'],
        'Species Name': row.get('Species Name', ''),
        'Change': 'status',
        'Baseline': old_status,
        'Current': new_status,
        'Severity': SEVERITY_CRITICAL if regressed else SEVERITY_INFO,
      })
    if old_upid != new_upid:
      changes.append({
        'Species ID': row['Species ID'],
        'Species Name': row.get('Species Name', ''),
        'Change': 'proteome_id',
        'Baseline': old_upid,
        'Current': new_upid,
        'Severity': SEVERITY_WARNING,
      })
  return changes


def delta_frame(deltas) -> pl.DataFrame:
  """Typed even when empty, so an identical build still writes a real artifact."""
  schema = {
    'Species ID': pl.String, 'Species Name': pl.String, 'Metric': pl.String,
    'Baseline': pl.Int64, 'Current': pl.Int64, 'Delta': pl.Int64, 'Severity': pl.String,
  }
  return pl.DataFrame(deltas, schema=schema).select(DELTA_COLUMNS)


def verdict(deltas, selection_changes) -> int:
  escalated = {SEVERITY_WARNING, SEVERITY_CRITICAL}
  if any(d['Severity'] in escalated for d in deltas):
    return EXIT_ESCALATED
  if any(c['Severity'] in escalated for c in selection_changes):
    return EXIT_ESCALATED
  if deltas or selection_changes:
    return EXIT_CHANGED
  return EXIT_IDENTICAL


def summarize(deltas, selection_changes, weeks) -> str:
  if not weeks or len(weeks) < 2:
    return 'no baseline week to compare against'
  critical = sum(1 for d in deltas if d['Severity'] == SEVERITY_CRITICAL)
  warning = sum(1 for d in deltas if d['Severity'] == SEVERITY_WARNING)
  status = sum(1 for c in selection_changes if c['Change'] == 'status')
  upid = sum(1 for c in selection_changes if c['Change'] == 'proteome_id')
  return (
    f'{weeks[0]} -> {weeks[1]}: {len(deltas)} metric deltas '
    f'({critical} CRITICAL, {warning} WARNING), '
    f'{status} status changes, {upid} proteome ID flips'
  )


def main(argv=None) -> int:
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-b', '--build_path', type=Path, default=REPO_ROOT / 'build')
  parser.add_argument('-s', '--stats', type=Path,
                      default=REPO_ROOT / 'protein-tree-weekly-stats.tsv')
  parser.add_argument('-o', '--output', type=Path,
                      help='per-species delta TSV (default build/reports/diff-metrics.tsv)')
  args = parser.parse_args(argv)

  if not args.stats.exists():
    print(f'diff-metrics: no stats sheet at {args.stats}', file=sys.stderr)
    return EXIT_ESCALATED

  sheet = pl.read_csv(args.stats, separator='\t', infer_schema_length=0)
  thresholds = Thresholds.from_env()
  deltas, weeks = compare_metrics(sheet, thresholds)

  reports_dir = args.build_path / 'reports'
  selection_changes = compare_selection(
    args.build_path / 'arborist' / REPORT_NAME, latest_archive(reports_dir)
  )

  output = args.output or reports_dir / 'diff-metrics.tsv'
  output.parent.mkdir(parents=True, exist_ok=True)
  delta_frame(deltas).write_csv(output, separator='\t')

  code = verdict(deltas, selection_changes)
  print(f'diff-metrics: {summarize(deltas, selection_changes, weeks)}', file=sys.stderr)
  for change in selection_changes:
    if change['Severity'] in (SEVERITY_WARNING, SEVERITY_CRITICAL):
      print(f"  {change['Severity']} {change['Species ID']} {change['Change']}: "
            f"{change['Baseline']} -> {change['Current']}", file=sys.stderr)
  print(f'diff-metrics: wrote {output}, exit {code}', file=sys.stderr)
  return code


if __name__ == '__main__':
  sys.exit(main())
