import html as html_lib
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


# ---- every overview metric is graded on week-over-week relative change ----

def _totals(weeks, **series):
  """TOTAL rows for the overview metrics; series maps metric -> (prev, curr)."""
  return _sheet([('TOTAL', '', m, p, c) for m, (p, c) in series.items()], weeks=weeks)


# The real TOTAL rows of the weekly stats sheet (prod, 2026-11-01).
W40_W41 = dict(
  source_antigens_total=(241502, 241503), source_antigens_assigned=(238476, 238477),
  epitopes_total=(2393207, 2393285), epitopes_assigned=(2381226, 2381304),
  epitopes_located_exact=(1727972, 1728041),
)
W41_W44 = dict(
  source_antigens_total=(241503, 241534), source_antigens_assigned=(238477, 238505),
  epitopes_total=(2393285, 2393386), epitopes_assigned=(2381304, 2381431),
  epitopes_located_exact=(1728041, 329248),   # -81%: shipped as "green" on 2026-11-01
)


def _overview(series, **overrides):
  merged = {**series, **overrides}
  sheet = _totals(WEEKS, **merged)
  weeks = alert.week_columns(sheet)
  return alert.compute_overview_findings(sheet, weeks, TH)


def test_replayed_w44_located_exact_collapse_is_red():
  findings = _overview(W41_W44)
  by_metric = {f['metric']: f for f in findings}
  assert by_metric['epitopes_located_exact']['severity'] == SEVERITY_CRITICAL
  # only the collapse and the ratio trip; the +31/+28/+101/+127 moves stay green
  assert set(by_metric) == {'epitopes_located_exact', 'epitopes_located_exact/epitopes_assigned'}
  assert by_metric['epitopes_located_exact/epitopes_assigned']['severity'] == SEVERITY_CRITICAL
  assert overall_status(findings) == STATUS_CRITICAL
  line = alert.format_finding(by_metric['epitopes_located_exact'])
  assert 'epitopes_located_exact' in line and '-80.9%' in line and '1728041 -> 329248' in line and '[RED]' in line
  assert notify_email.subject_for_status(overall_status(findings), '2026-W44') == \
    'Arborist prod run is RED with a CRITICAL - 2026-W44'


def test_normal_week_is_green():
  assert _overview(W40_W41) == []


def test_normal_week_with_a_hundred_epitopes_either_way_is_green():
  assert _overview(W41_W44, epitopes_located_exact=(1728041, 1728041 + 100)) == []
  assert _overview(W41_W44, epitopes_located_exact=(1728041, 1728041 - 100),
                   epitopes_assigned=(2381304, 2381304 - 100)) == []


@pytest.mark.parametrize('prev,curr,expected', [
  (1000, 1019, None),                  # +1.9%
  (1000, 1020, SEVERITY_WARNING),      # +2%  yellow
  (1000, 1099, SEVERITY_WARNING),
  (1000, 1100, 'ORANGE'),              # +10% orange
  (1000, 1249, 'ORANGE'),
  (1000, 1250, SEVERITY_CRITICAL),     # +25% red
  (1000, 980, SEVERITY_WARNING),       # drops grade the same as rises
  (1000, 900, 'ORANGE'),
  (1000, 750, SEVERITY_CRITICAL),
  (1000, 0, SEVERITY_CRITICAL),
  (0, 50, None),                       # no baseline to compare against
])
def test_relative_change_levels(prev, curr, expected):
  assert alert.classify_change(prev, curr, TH) == expected


def test_every_overview_metric_is_graded_not_just_assigned():
  for metric in alert.OVERVIEW_METRICS:
    findings = _overview(W40_W41, **{metric: (1000, 700)})
    assert any(f['metric'] == metric and f['severity'] == SEVERITY_CRITICAL for f in findings), metric


def test_orange_status_and_banner_and_subject():
  findings = _overview(W40_W41, epitopes_total=(2393207, 2393207 + 250000))   # +10.4%
  assert overall_status(findings) == 'orange'
  assert 'ORANGE' in alert.render_banner('orange', findings, '2026-W02')
  assert notify_email.subject_for_status('orange', '2026-W02') == \
    'Arborist prod run is ORANGE with a WARNING - 2026-W02'


def test_ratio_collapse_is_red_even_when_each_total_moves_less_than_red():
  # located -20% (orange alone) while assigned is flat: the ratio falls 20% -> RED
  findings = _overview(W40_W41, epitopes_located_exact=(1000000, 800000), epitopes_assigned=(1400000, 1400000))
  by_metric = {f['metric']: f['severity'] for f in findings}
  assert by_metric['epitopes_located_exact'] == 'ORANGE'
  assert by_metric['epitopes_located_exact/epitopes_assigned'] == SEVERITY_CRITICAL
  assert overall_status(findings) == STATUS_CRITICAL


def test_ratio_threshold_is_configurable():
  th = Thresholds.from_env({'ALERT_RATIO_RED_PCT': '30'})
  sheet = _totals(WEEKS, **{**W40_W41, 'epitopes_located_exact': (1000000, 800000), 'epitopes_assigned': (1400000, 1400000)})
  findings = alert.compute_overview_findings(sheet, alert.week_columns(sheet), th)
  assert 'epitopes_located_exact/epitopes_assigned' not in {f['metric'] for f in findings}


def test_thresholds_come_from_env():
  th = Thresholds.from_env({'ALERT_YELLOW_PCT': '5', 'ALERT_ORANGE_PCT': '20', 'ALERT_RED_PCT': '50'})
  assert alert.classify_change(1000, 1030, th) is None
  assert alert.classify_change(1000, 1200, th) == 'ORANGE'
  assert alert.classify_change(1000, 1500, th) == SEVERITY_CRITICAL


def test_missing_metric_row_is_red():
  series = {k: v for k, v in W40_W41.items() if k != 'epitopes_located_exact'}
  findings = _overview(series)
  missing = [f for f in findings if f['kind'] == 'missing']
  assert [f['metric'] for f in missing] == ['epitopes_located_exact']
  assert missing[0]['severity'] == SEVERITY_CRITICAL
  assert 'MISSING' in alert.format_finding(missing[0])


def test_missing_current_value_is_red():
  sheet = _totals(WEEKS, **W40_W41)
  sheet = sheet.with_columns(
    pl.when(pl.col('Metric') == 'epitopes_total').then(None).otherwise(pl.col(WEEKS[1])).alias(WEEKS[1]))
  findings = alert.compute_overview_findings(sheet, list(WEEKS), TH)
  assert [(f['metric'], f['kind']) for f in findings] == [('epitopes_total', 'missing')]


def test_run_replay_writes_red_html_naming_the_metric(tmp_path):
  sheet = _totals(WEEKS, **W41_W44)
  alert.run(_args(tmp_path, sheet), env={'ALERT_EMAIL': '0'})
  html = html_lib.unescape((tmp_path / 'alert.html').read_text())
  assert 'CRITICAL' in html
  assert 'epitopes_located_exact -80.9% (1728041 -> 329248, -1398793) [RED]' in html
  assert 'ratio' in html


def test_run_replay_sends_red_subject(tmp_path, monkeypatch):
  sent = []
  monkeypatch.setattr(notify_email, 'send_email', lambda message, **k: sent.append(message))
  alert.run(_args(tmp_path, _totals(WEEKS, **W41_W44)), env={'ALERT_EMAIL': '1'})
  assert sent[0]['Subject'] == 'Arborist prod run is RED with a CRITICAL - 2026-W02'


def test_run_normal_week_is_green_subject(tmp_path, monkeypatch):
  sent = []
  monkeypatch.setattr(notify_email, 'send_email', lambda message, **k: sent.append(message))
  alert.run(_args(tmp_path, _totals(WEEKS, **W40_W41)), env={'ALERT_EMAIL': '1'})
  assert sent[0]['Subject'] == 'Arborist prod run is green - 2026-W02'


def test_total_is_not_graded_by_the_absolute_drop_rule_in_run(tmp_path, monkeypatch):
  # a 1,500-epitope drop on 2.4M (0.06%) used to be CRITICAL via agg_abs; it is green by relative change
  sent = []
  monkeypatch.setattr(notify_email, 'send_email', lambda message, **k: sent.append(message))
  series = {**W40_W41, 'epitopes_assigned': (2381226, 2381226 - 1500)}
  alert.run(_args(tmp_path, _totals(WEEKS, **series)), env={'ALERT_EMAIL': '1'})
  assert sent[0]['Subject'] == 'Arborist prod run is green - 2026-W02'


def test_per_species_located_exact_collapse_is_critical():
  sheet = _sheet([('9606', 'Homo sapiens', 'epitopes_located_exact', 1193335, 6476)])
  findings, _ = compute_findings(sheet, TH)
  assert findings[0]['severity'] == SEVERITY_CRITICAL
