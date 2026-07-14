#!/usr/bin/env python3
"""Weekly parent-assignment alert, derived from the master stats sheet.

Runs at the END of `make protein`, after stats.py. Compares the two most recent
ISO-week columns of protein-tree-weekly-stats.tsv and grades week-to-week DROPS in
assigned counts by a severity that scales with the number:

  - a 1-epitope species dropping to 0 is INFO (logged, not escalated),
  - a big species losing ~10% is CRITICAL ("human -10%").

It always sends a weekly digest email (traffic-light subject), prints an stderr
banner, and writes an HTML report. Advisory only: ALWAYS exits 0, so it can never
break or cancel the run.
"""

import argparse
import html as html_lib
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from protein_tree import notify_email
from protein_tree.notify_email import STATUS_GREEN, STATUS_WARNING, STATUS_CRITICAL

REPO_ROOT = Path(__file__).parents[3]
ID_COLUMNS = ['Species ID', 'Species Name', 'Metric']
ALERT_METRICS = ['epitopes_assigned', 'source_antigens_assigned']
OVERVIEW_METRICS = [
  'source_antigens_total', 'source_antigens_assigned',
  'epitopes_total', 'epitopes_assigned', 'epitopes_located_exact',
]
EMPTY_STATUSES = {'EMPTY', 'MISSING'}

SEVERITY_INFO = 'INFO'
SEVERITY_WARNING = 'WARNING'
SEVERITY_CRITICAL = 'CRITICAL'
ESCALATED = (SEVERITY_WARNING, SEVERITY_CRITICAL)


@dataclass
class Thresholds:
  hard_abs: int = 1000     # per-species absolute drop -> CRITICAL
  big: int = 1000          # "large" baseline
  crit_pct: float = 10.0   # % drop on a large baseline -> CRITICAL
  warn_abs: int = 200      # per-species absolute drop -> WARNING
  mid: int = 100           # "moderate" baseline
  warn_pct: float = 25.0   # % drop on a moderate baseline -> WARNING
  agg_abs: int = 1000      # TOTAL absolute drop -> CRITICAL
  agg_pct: float = 5.0     # TOTAL % drop -> CRITICAL

  @classmethod
  def from_env(cls, env=None):
    env = env if env is not None else os.environ
    def num(name, default, cast):
      raw = env.get(name)
      if raw is None or raw == '':
        return default
      try:
        return cast(raw)
      except ValueError:
        return default
    return cls(
      hard_abs=num('ALERT_HARD_ABS', cls.hard_abs, int),
      big=num('ALERT_BIG', cls.big, int),
      crit_pct=num('ALERT_CRIT_PCT', cls.crit_pct, float),
      warn_abs=num('ALERT_WARN_ABS', cls.warn_abs, int),
      mid=num('ALERT_MID', cls.mid, int),
      warn_pct=num('ALERT_WARN_PCT', cls.warn_pct, float),
      agg_abs=num('ALERT_AGG_ABS', cls.agg_abs, int),
      agg_pct=num('ALERT_AGG_PCT', cls.agg_pct, float),
    )


def classify_severity(prev: int, curr: int, thresholds: Thresholds, is_total: bool = False):
  """Grade a drop from prev->curr. Returns a severity or None (no drop)."""
  drop = prev - curr
  if drop <= 0:
    return None
  pct = (drop / prev * 100.0) if prev > 0 else 100.0
  if is_total:
    if drop >= thresholds.agg_abs or pct >= thresholds.agg_pct:
      return SEVERITY_CRITICAL
    return SEVERITY_INFO
  if drop >= thresholds.hard_abs or (prev >= thresholds.big and pct >= thresholds.crit_pct):
    return SEVERITY_CRITICAL
  if drop >= thresholds.warn_abs or (prev >= thresholds.mid and pct >= thresholds.warn_pct):
    return SEVERITY_WARNING
  return SEVERITY_INFO


def _to_int(value) -> int:
  if value is None or value == '':
    return 0
  try:
    return int(value)
  except (ValueError, TypeError):
    return 0


def week_columns(sheet: pl.DataFrame) -> list:
  return sorted(c for c in sheet.columns if c not in ID_COLUMNS)


def load_empty_species(build_path: Path) -> set:
  """Species flagged EMPTY/MISSING in the proteome-selection report, if present."""
  report = build_path / 'arborist' / 'proteome-selection-report.tsv'
  if not report.exists():
    return set()
  df = pl.read_csv(report, separator='\t', infer_schema_length=0)
  return set(df.filter(pl.col('Status').is_in(list(EMPTY_STATUSES)))['Species ID'])


def compute_findings(sheet: pl.DataFrame, thresholds: Thresholds, empty_species: set = None):
  """Return (findings, [prev_week, curr_week]). Empty findings when <2 weeks."""
  weeks = week_columns(sheet)
  if len(weeks) < 2:
    return [], weeks
  prev_week, curr_week = weeks[-2], weeks[-1]
  empty_species = empty_species or set()

  findings = []
  for row in sheet.iter_rows(named=True):
    if row['Metric'] not in ALERT_METRICS:
      continue
    is_total = row['Species ID'] == 'TOTAL'
    prev = _to_int(row[prev_week])
    curr = _to_int(row[curr_week])
    severity = classify_severity(prev, curr, thresholds, is_total=is_total)
    if severity is None:
      continue
    drop = prev - curr
    findings.append({
      'species_id': row['Species ID'],
      'species_name': row['Species Name'],
      'metric': row['Metric'],
      'prev': prev,
      'curr': curr,
      'delta': -drop,
      'pct': (drop / prev * 100.0) if prev > 0 else 100.0,
      'severity': severity,
      'empty_proteome': row['Species ID'] in empty_species,
      'is_total': is_total,
    })
  return findings, [prev_week, curr_week]


def overall_status(findings: list) -> str:
  severities = {f['severity'] for f in findings}
  if SEVERITY_CRITICAL in severities:
    return STATUS_CRITICAL
  if SEVERITY_WARNING in severities:
    return STATUS_WARNING
  return STATUS_GREEN


def _scope(finding: dict) -> str:
  if finding['is_total']:
    return 'TOTAL'
  return f"{finding['species_name']} ({finding['species_id']})"


def format_finding(finding: dict) -> str:
  """One-line 'species (id) - metric -N% (prev -> curr) [SEVERITY]'."""
  line = (
    f"{_scope(finding)} - {finding['metric']} -{finding['pct']:.0f}% "
    f"({finding['prev']} -> {finding['curr']}) [{finding['severity']}]"
  )
  if finding['empty_proteome']:
    line += ' [empty proteome]'
  return line


def render_banner(status: str, findings: list, week: str) -> str:
  escalated = [f for f in findings if f['severity'] in ESCALATED]
  marker = {
    STATUS_CRITICAL: 'ARBORIST-ALERT[CRITICAL]',
    STATUS_WARNING: 'ARBORIST-ALERT[WARNING]',
    STATUS_GREEN: 'ARBORIST-OK',
  }[status]
  lines = [f'{marker}: Arborist {week} - {status.upper()}']
  for finding in escalated:
    lines.append('  ' + format_finding(finding))
  return '\n'.join(lines)


def _overview_rows(sheet: pl.DataFrame, weeks: list) -> list:
  """TOTAL this-week vs last-week for the overview metrics."""
  if not weeks:
    return []
  curr_week = weeks[-1]
  prev_week = weeks[-2] if len(weeks) > 1 else None
  totals = sheet.filter(pl.col('Species ID') == 'TOTAL')
  by_metric = {r['Metric']: r for r in totals.iter_rows(named=True)}
  rows = []
  for metric in OVERVIEW_METRICS:
    row = by_metric.get(metric)
    if row is None:
      continue
    curr = _to_int(row[curr_week])
    prev = _to_int(row[prev_week]) if prev_week else None
    rows.append((metric, prev, curr))
  return rows


def render_html(sheet: pl.DataFrame, findings: list, weeks: list, status: str, week: str) -> str:
  esc = html_lib.escape
  escalated = [f for f in findings if f['severity'] in ESCALATED]
  parts = [
    '<html><body style="font-family:sans-serif">',
    f'<h2>Arborist prod run - {esc(week)} - {esc(status.upper())}</h2>',
  ]

  parts.append('<h3>Overall stats</h3><table border="1" cellpadding="4" cellspacing="0">')
  parts.append('<tr><th>Metric</th><th>Previous</th><th>This week</th><th>Delta</th></tr>')
  for metric, prev, curr in _overview_rows(sheet, weeks):
    delta = '' if prev is None else f'{curr - prev:+d}'
    prev_txt = '' if prev is None else str(prev)
    parts.append(
      f'<tr><td>{esc(metric)}</td><td>{prev_txt}</td><td>{curr}</td><td>{delta}</td></tr>'
    )
  parts.append('</table>')

  if escalated:
    parts.append(f'<h3>Findings ({len(escalated)})</h3><ul>')
    for finding in escalated:
      parts.append(f'<li>{esc(format_finding(finding))}</li>')
    parts.append('</ul>')
  else:
    parts.append('<p>No week-to-week warnings or criticals.</p>')

  parts.append('</body></html>')
  return '\n'.join(parts)


def _send_digest(status, findings, weeks, week, sheet, env):
  """Build + send the digest email. Raises on send failure (caller swallows)."""
  html_body = render_html(sheet, findings, weeks, status, week)
  subject = notify_email.subject_for_status(status, week)
  message = notify_email.build_message(
    subject, html_body,
    to=env.get('ALERT_EMAIL_TO', notify_email.DEFAULT_EMAIL_TO),
    sender=env.get('ALERT_EMAIL_FROM', notify_email.DEFAULT_EMAIL_FROM),
  )
  host = env.get('ALERT_SMTP_HOST', notify_email.DEFAULT_SMTP_HOST)
  port = int(env.get('ALERT_SMTP_PORT', notify_email.DEFAULT_SMTP_PORT) or notify_email.DEFAULT_SMTP_PORT)
  notify_email.send_email(message, host=host, port=port)


def run(args, env=None) -> None:
  env = env if env is not None else os.environ
  stats_file = Path(args.stats_file)
  build_path = Path(args.build_path)
  html_out = Path(args.html_out)

  if not stats_file.exists():
    print(f'alert: no stats sheet at {stats_file}, nothing to compare', file=sys.stderr)
    return

  sheet = pl.read_csv(stats_file, separator='\t', infer_schema_length=0)
  thresholds = Thresholds.from_env(env)
  empty_species = load_empty_species(build_path)
  findings, weeks = compute_findings(sheet, thresholds, empty_species)

  if len(weeks) < 2:
    print('alert: baseline established (only one week of stats), no comparison', file=sys.stderr)

  week = weeks[-1] if weeks else '?'
  status = overall_status(findings)

  print(render_banner(status, findings, week), file=sys.stderr)

  html_out.parent.mkdir(parents=True, exist_ok=True)
  html_out.write_text(render_html(sheet, findings, weeks, status, week))

  if env.get('ALERT_EMAIL', '1') == '1':
    _send_digest(status, findings, weeks, week, sheet, env)
    print(f'alert: emailed weekly digest ({status})', file=sys.stderr)
  else:
    print('alert: ALERT_EMAIL disabled, skipping email', file=sys.stderr)


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('-b', '--build_path', type=str, default=str(REPO_ROOT / 'build'),
                      help='Path for all Arborist build files (for the proteome report).')
  parser.add_argument('-o', '--stats_file', type=str,
                      default=str(REPO_ROOT / 'protein-tree-weekly-stats.tsv'),
                      help='Master weekly stats sheet to read.')
  parser.add_argument('--html_out', type=str,
                      default=str(REPO_ROOT / 'build' / 'arborist' / 'protein-tree-alert.html'),
                      help='Where to write the HTML digest.')
  args = parser.parse_args()

  # Advisory only: log and exit 0 on ANY error so the weekly run never breaks.
  try:
    run(args)
  except Exception as error:  # noqa: BLE001 - deliberately non-breaking
    print(f'alert: non-fatal error, alert skipped: {error}', file=sys.stderr)
  return 0


if __name__ == '__main__':
  sys.exit(main())
