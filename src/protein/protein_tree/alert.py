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
from protein_tree.notify_email import STATUS_GREEN, STATUS_WARNING, STATUS_ORANGE, STATUS_CRITICAL
from protein_tree import skipped_proteomes

REPO_ROOT = Path(__file__).parents[3]
ID_COLUMNS = ['Species ID', 'Species Name', 'Metric']
# Per-species drop grading. epitopes_located_exact is graded too: an index/matcher failure empties
# it while epitopes_assigned stays flat (2026-11-01: 1,728,041 -> 329,248 went out as "green").
ALERT_METRICS = ['epitopes_assigned', 'source_antigens_assigned', 'epitopes_located_exact']
OVERVIEW_METRICS = [
  'source_antigens_total', 'source_antigens_assigned',
  'epitopes_total', 'epitopes_assigned', 'epitopes_located_exact',
]
EMPTY_STATUSES = {'EMPTY', 'MISSING'}

SEVERITY_INFO = 'INFO'
SEVERITY_WARNING = 'WARNING'    # shown as YELLOW
SEVERITY_ORANGE = 'ORANGE'
SEVERITY_CRITICAL = 'CRITICAL'  # shown as RED
ESCALATED = (SEVERITY_WARNING, SEVERITY_ORANGE, SEVERITY_CRITICAL)
LEVEL_LABEL = {
  SEVERITY_INFO: 'INFO', SEVERITY_WARNING: 'YELLOW', SEVERITY_ORANGE: 'ORANGE', SEVERITY_CRITICAL: 'RED',
}
RATIO_METRIC = ('epitopes_located_exact', 'epitopes_assigned')  # numerator, denominator


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
  # TOTAL metrics: week-over-week relative change, rises and drops alike
  yellow_pct: float = 2.0    # >= -> YELLOW (WARNING)
  orange_pct: float = 10.0   # >= -> ORANGE
  red_pct: float = 25.0      # >= -> RED (CRITICAL)
  ratio_red_pct: float = 10.0  # located_exact / assigned ratio falling by more than this -> RED

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
      yellow_pct=num('ALERT_YELLOW_PCT', cls.yellow_pct, float),
      orange_pct=num('ALERT_ORANGE_PCT', cls.orange_pct, float),
      red_pct=num('ALERT_RED_PCT', cls.red_pct, float),
      ratio_red_pct=num('ALERT_RATIO_RED_PCT', cls.ratio_red_pct, float),
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


def classify_change(prev: int, curr: int, thresholds: Thresholds):
  """Grade a week-over-week relative change prev->curr, rises and drops alike (TOTAL metrics).
  Returns YELLOW(WARNING) / ORANGE / RED(CRITICAL), or None below the yellow threshold."""
  if prev <= 0:
    return None
  pct = abs(curr - prev) / prev * 100.0
  if pct >= thresholds.red_pct:
    return SEVERITY_CRITICAL
  if pct >= thresholds.orange_pct:
    return SEVERITY_ORANGE
  if pct >= thresholds.yellow_pct:
    return SEVERITY_WARNING
  return None


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


def compute_findings(sheet: pl.DataFrame, thresholds: Thresholds, empty_species: set = None,
                     include_total: bool = True):
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
    if is_total and not include_total:
      continue
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


def _opt_int(value):
  """int(value), or None when the cell is absent/empty/non-numeric (unlike _to_int, which maps those to 0)."""
  if value is None or value == '':
    return None
  try:
    return int(value)
  except (ValueError, TypeError):
    return None


def _total_finding(metric, prev, curr, severity, kind='change', pct=None):
  if pct is None:
    pct = (abs(curr - prev) / prev * 100.0) if prev else 100.0
  return {
    'species_id': 'TOTAL', 'species_name': '', 'metric': metric,
    'prev': prev, 'curr': curr, 'delta': curr - prev if kind != 'ratio' else 0,
    'pct': pct, 'severity': severity, 'empty_proteome': False, 'is_total': True, 'kind': kind,
  }


def compute_overview_findings(sheet: pl.DataFrame, weeks: list, thresholds: Thresholds) -> list:
  """Grade every OVERVIEW_METRICS TOTAL on week-over-week relative change (rises and drops),
  flag a metric missing from the sheet, and flag a collapse of located_exact / assigned."""
  if len(weeks) < 2:
    return []
  prev_week, curr_week = weeks[-2], weeks[-1]
  totals = {
    r['Metric']: r for r in sheet.filter(pl.col('Species ID') == 'TOTAL').iter_rows(named=True)
  }
  findings, values = [], {}
  for metric in OVERVIEW_METRICS:
    row = totals.get(metric)
    prev = _opt_int(row[prev_week]) if row else None
    curr = _opt_int(row[curr_week]) if row else None
    values[metric] = (prev, curr)
    if curr is None:
      findings.append(_total_finding(metric, prev or 0, 0, SEVERITY_CRITICAL, kind='missing'))
      continue
    if prev is None:
      continue  # first appearance of the metric: nothing to compare against
    severity = classify_change(prev, curr, thresholds)
    if severity:
      findings.append(_total_finding(metric, prev, curr, severity))

  num, den = RATIO_METRIC
  (pn, cn), (pd, cd) = values[num], values[den]
  if None not in (pn, cn, pd, cd) and pd > 0 and cd > 0:
    prev_ratio, curr_ratio = pn / pd, cn / cd
    if prev_ratio > 0:
      drop = (prev_ratio - curr_ratio) / prev_ratio * 100.0
      if drop > thresholds.ratio_red_pct:
        findings.append(_total_finding(
          f'{num}/{den}', prev_ratio, curr_ratio, SEVERITY_CRITICAL, kind='ratio', pct=drop))
  return findings


def overall_status(findings: list) -> str:
  severities = {f['severity'] for f in findings}
  if SEVERITY_CRITICAL in severities:
    return STATUS_CRITICAL
  if SEVERITY_ORANGE in severities:
    return STATUS_ORANGE
  if SEVERITY_WARNING in severities:
    return STATUS_WARNING
  return STATUS_GREEN


def _scope(finding: dict) -> str:
  if finding['is_total']:
    return 'TOTAL'
  return f"{finding['species_name']} ({finding['species_id']})"


def format_finding(finding: dict) -> str:
  """One line naming the scope, the metric that tripped, by how much, and the level."""
  level = LEVEL_LABEL[finding['severity']]
  kind = finding.get('kind', 'change')
  if kind == 'missing':
    line = f"{_scope(finding)} - {finding['metric']} MISSING from the stats sheet this week (previous {finding['prev']}) [{level}]"
  elif kind == 'ratio':
    line = (
      f"{_scope(finding)} - {finding['metric']} ratio {finding['prev']:.4f} -> {finding['curr']:.4f} "
      f"(-{finding['pct']:.1f}%) [{level}]"
    )
  else:
    sign = '-' if finding['delta'] < 0 else '+'
    line = (
      f"{_scope(finding)} - {finding['metric']} {sign}{finding['pct']:.1f}% "
      f"({finding['prev']} -> {finding['curr']}, {finding['delta']:+d}) [{level}]"
    )
  if finding['empty_proteome']:
    line += ' [empty proteome]'
  return line


def render_banner(status: str, findings: list, week: str) -> str:
  escalated = [f for f in findings if f['severity'] in ESCALATED]
  marker = {
    STATUS_CRITICAL: 'ARBORIST-ALERT[CRITICAL]',
    STATUS_ORANGE: 'ARBORIST-ALERT[ORANGE]',
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


def render_html(sheet: pl.DataFrame, findings: list, weeks: list, status: str, week: str, skips: list = None) -> str:
  esc = html_lib.escape
  escalated = [f for f in findings if f['severity'] in ESCALATED]
  parts = [
    '<html><body style="font-family:sans-serif">',
    f'<h2>Arborist prod run - {esc(week)} - {esc(status.upper())}</h2>',
  ]

  parts.append('<h3>Overall stats</h3><table border="1" cellpadding="4" cellspacing="0">')
  parts.append('<tr><th>Metric</th><th>Previous</th><th>This week</th><th>Delta</th><th>Change</th><th>Level</th></tr>')
  levels = {f['metric']: f['severity'] for f in findings if f['is_total']}
  overview = _overview_rows(sheet, weeks)
  for metric, prev, curr in overview:
    delta = '' if prev is None else f'{curr - prev:+d}'
    prev_txt = '' if prev is None else str(prev)
    pct_txt = '' if not prev else f'{(curr - prev) / prev * 100:+.2f}%'
    level = LEVEL_LABEL.get(levels.get(metric), 'ok')
    parts.append(
      f'<tr><td>{esc(metric)}</td><td>{prev_txt}</td><td>{curr}</td><td>{delta}</td>'
      f'<td>{pct_txt}</td><td>{esc(level)}</td></tr>'
    )
  by_metric = {m: (p, c) for m, p, c in overview}
  num, den = RATIO_METRIC
  if num in by_metric and den in by_metric:
    (pn, cn), (pd, cd) = by_metric[num], by_metric[den]
    prev_ratio = f'{pn / pd:.4f}' if pn is not None and pd else ''
    curr_ratio = f'{cn / cd:.4f}' if cd else ''
    level = LEVEL_LABEL.get(levels.get(f'{num}/{den}'), 'ok')
    parts.append(
      f'<tr><td>{esc(num)} / {esc(den)}</td><td>{prev_ratio}</td><td>{curr_ratio}</td><td></td><td></td>'
      f'<td>{esc(level)}</td></tr>'
    )
  parts.append('</table>')

  if escalated:
    parts.append(f'<h3>Findings ({len(escalated)})</h3><ul>')
    for finding in escalated:
      parts.append(f'<li>{esc(format_finding(finding))}</li>')
    parts.append('</ul>')
  elif skips:
    parts.append('<p>No metric findings; the WARNING comes from the skipped UniProt proteomes below.</p>')
  else:
    parts.append('<p>No week-to-week warnings or criticals.</p>')

  if skips:
    parts.append(f'<h3>Skipped UniProt proteomes ({len(skips)})</h3>')
    parts.append('<p>These proteome fetches failed after all retries and were skipped instead of failing the build. '
                 'They are retried on the next build.</p>')
    parts.append('<table border="1" cellpadding="4" cellspacing="0">')
    parts.append('<tr><th>Species</th><th>Proteome</th><th>Outcome</th><th>Error</th></tr>')
    for row in skips:
      parts.append(
        f"<tr><td>{esc(row['Species Name'])} ({esc(row['Species ID'])})</td><td>{esc(row['Proteome ID'] or '-')}</td>"
        f"<td>{esc(row['Outcome'])}</td><td>{esc(row['Error'])}</td></tr>"
      )
    parts.append('</table>')

  parts.append('</body></html>')
  return '\n'.join(parts)


def _send_digest(status, findings, weeks, week, sheet, env, skips=None):
  """Build + send the digest email. Raises on send failure (caller swallows)."""
  html_body = render_html(sheet, findings, weeks, status, week, skips)
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
  findings, weeks = compute_findings(sheet, thresholds, empty_species, include_total=False)
  # TOTAL rows are graded by week-over-week relative change (all overview metrics), not by
  # the per-species absolute-drop rules, so +-100 epitopes on 2.4M stays green.
  findings = compute_overview_findings(sheet, weeks, thresholds) + findings

  if len(weeks) < 2:
    print('alert: baseline established (only one week of stats), no comparison', file=sys.stderr)

  week = weeks[-1] if weeks else '?'
  status = overall_status(findings)
  # A skipped proteome silently loses that taxon's assignments for the week, so
  # a run with skips is never plain green.
  skips = skipped_proteomes.load_skips(build_path)
  if skips and status == STATUS_GREEN:
    status = STATUS_WARNING

  print(render_banner(status, findings, week), file=sys.stderr)
  if skips:
    print(skipped_proteomes.format_skips(skips), file=sys.stderr)

  html_out.parent.mkdir(parents=True, exist_ok=True)
  html_out.write_text(render_html(sheet, findings, weeks, status, week, skips))

  if env.get('ALERT_EMAIL', '1') == '1':
    _send_digest(status, findings, weeks, week, sheet, env, skips)
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
