import sys
import types

import polars as pl
import pytest

from protein_tree import alert, notify_email
from protein_tree.alert import (
  Thresholds,
  classify_severity,
  compute_findings,
  overall_status,
  SEVERITY_INFO,
  SEVERITY_WARNING,
  SEVERITY_CRITICAL,
)
from protein_tree.notify_email import STATUS_GREEN, STATUS_WARNING, STATUS_CRITICAL

TH = Thresholds()
WEEKS = ('2026-W01', '2026-W02')


def _sheet(rows, weeks=WEEKS):
  """rows: (species_id, species_name, metric, prev, curr) with str/None cells."""
  data = {'Species ID': [], 'Species Name': [], 'Metric': [], weeks[0]: [], weeks[1]: []}
  for sid, name, metric, prev, curr in rows:
    data['Species ID'].append(sid)
    data['Species Name'].append(name)
    data['Metric'].append(metric)
    data[weeks[0]].append(None if prev is None else str(prev))
    data[weeks[1]].append(None if curr is None else str(curr))
  schema = {'Species ID': pl.String, 'Species Name': pl.String, 'Metric': pl.String,
            weeks[0]: pl.String, weeks[1]: pl.String}
  return pl.DataFrame(data, schema=schema)


# ---- severity tiers (scale with the number) ----

def test_tiny_species_losing_only_epitope_is_info():
  assert classify_severity(1, 0, TH) == SEVERITY_INFO


def test_human_scale_ten_percent_is_critical():
  # 48210 -> 42400 : ~12% and abs > 1000
  assert classify_severity(48210, 42400, TH) == SEVERITY_CRITICAL


def test_ten_percent_on_large_baseline_is_critical_even_under_abs():
  # drop 500 (< hard_abs) but prev >= big and pct >= 10
  assert classify_severity(5000, 4500, TH) == SEVERITY_CRITICAL


def test_moderate_percent_drop_is_warning():
  # prev 300 -> 200 : 33% on a moderate baseline
  assert classify_severity(300, 200, TH) == SEVERITY_WARNING


def test_absolute_drop_is_warning_even_at_low_percent():
  # drop 250 (>= warn_abs) though only 0.25%
  assert classify_severity(100000, 99750, TH) == SEVERITY_WARNING


def test_increase_and_no_change_are_not_findings():
  assert classify_severity(100, 120, TH) is None
  assert classify_severity(100, 100, TH) is None


def test_total_uses_aggregate_thresholds():
  assert classify_severity(100000, 94000, TH, is_total=True) == SEVERITY_CRITICAL  # 6% and >1000
  assert classify_severity(1000, 940, TH, is_total=True) == SEVERITY_CRITICAL       # 6% >= agg_pct
  assert classify_severity(100000, 99900, TH, is_total=True) == SEVERITY_INFO       # 0.1%, <1000


# ---- findings + overall status ----

def test_compute_findings_and_status():
  sheet = _sheet([
    ('TOTAL', '', 'epitopes_assigned', 68000, 60000),        # -8000 total -> CRITICAL
    ('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 4500),  # -10% large -> CRITICAL
    ('500', 'Some sp', 'epitopes_assigned', 300, 200),       # -33% moderate -> WARNING
    ('7', 'Tiny sp', 'epitopes_assigned', 1, 0),             # -> INFO
    ('9', 'Growing sp', 'epitopes_assigned', 100, 150),      # increase -> no finding
  ])
  findings, weeks = compute_findings(sheet, TH)
  assert weeks == ['2026-W01', '2026-W02']
  by_species = {f['species_id']: f['severity'] for f in findings}
  assert by_species['TOTAL'] == SEVERITY_CRITICAL
  assert by_species['9606'] == SEVERITY_CRITICAL
  assert by_species['500'] == SEVERITY_WARNING
  assert by_species['7'] == SEVERITY_INFO
  assert '9' not in by_species  # increase produced no finding
  assert overall_status(findings) == STATUS_CRITICAL


def test_status_warning_when_no_critical():
  sheet = _sheet([('500', 'Some sp', 'epitopes_assigned', 300, 200)])
  findings, _ = compute_findings(sheet, TH)
  assert overall_status(findings) == STATUS_WARNING


def test_status_green_when_only_info():
  sheet = _sheet([('7', 'Tiny sp', 'epitopes_assigned', 1, 0)])
  findings, _ = compute_findings(sheet, TH)
  assert overall_status(findings) == STATUS_GREEN


def test_single_week_is_baseline_no_findings():
  sheet = pl.DataFrame(
    {'Species ID': ['TOTAL'], 'Species Name': [''], 'Metric': ['epitopes_assigned'],
     '2026-W01': ['100']},
    schema={'Species ID': pl.String, 'Species Name': pl.String, 'Metric': pl.String, '2026-W01': pl.String},
  )
  findings, weeks = compute_findings(sheet, TH)
  assert findings == []
  assert len(weeks) < 2


def test_empty_proteome_annotation():
  sheet = _sheet([('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 0)])
  findings, _ = compute_findings(sheet, TH, empty_species={'9606'})
  assert findings[0]['empty_proteome'] is True
  assert '[empty proteome]' in alert.format_finding(findings[0])


# ---- non-breaking guarantees ----

def _args(tmp_path, sheet):
  stats_file = tmp_path / 'stats.tsv'
  sheet.write_csv(stats_file, separator='\t')
  return types.SimpleNamespace(
    stats_file=str(stats_file),
    build_path=str(tmp_path / 'build'),
    html_out=str(tmp_path / 'alert.html'),
  )


def test_run_writes_html_without_email(tmp_path):
  sheet = _sheet([
    ('TOTAL', '', 'epitopes_assigned', 68000, 60000),
    ('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 4500),
  ])
  alert.run(_args(tmp_path, sheet), env={'ALERT_EMAIL': '0'})
  html = (tmp_path / 'alert.html').read_text()
  assert 'CRITICAL' in html
  assert 'Homo sapiens (9606)' in html


def test_main_exit_zero_on_generic_error(monkeypatch):
  monkeypatch.setattr(sys, 'argv', ['alert.py'])
  monkeypatch.setattr(alert, 'run', lambda *a, **k: (_ for _ in ()).throw(RuntimeError('boom')))
  assert alert.main() == 0


def test_main_exit_zero_when_email_send_fails(tmp_path, monkeypatch):
  sheet = _sheet([
    ('TOTAL', '', 'epitopes_assigned', 68000, 60000),
    ('9606', 'Homo sapiens', 'epitopes_assigned', 5000, 4500),
  ])
  args = _args(tmp_path, sheet)
  monkeypatch.setattr(sys, 'argv',
                      ['alert.py', '-o', args.stats_file, '-b', args.build_path,
                       '--html_out', args.html_out])
  monkeypatch.setenv('ALERT_EMAIL', '1')

  def boom(*a, **k):
    raise OSError('relay down')

  monkeypatch.setattr(notify_email, 'send_email', boom)
  assert alert.main() == 0
  # HTML digest still written before the (failed) send
  assert (tmp_path / 'alert.html').exists()
