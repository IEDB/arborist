"""Guards on ops/uniprot-cron/refresh-run.sh that do not need arborist-dev."""

import re
import subprocess
from datetime import date, timedelta
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent.parent / 'ops' / 'uniprot-cron' / 'refresh-run.sh'
SOURCE = SCRIPT.read_text()


def test_lookback_uses_the_expression_these_tests_model():
  """`date -d 'last sunday'` is ambiguous on a Sunday; pin the %w arithmetic."""
  assert 'date +%w' in SOURCE
  assert re.search(r'week \* 7\s*\)\)\s*days', SOURCE)
  assert 'last sunday' not in SOURCE


def sunday_from_shell(reference: date, week: int) -> str:
  """Run the script's own date expression with `date` pinned to a reference day."""
  weekday = int(subprocess.run(['date', '-d', reference.isoformat(), '+%w'],
                               capture_output=True, text=True, check=True).stdout)
  offset = weekday + week * 7
  out = subprocess.run(['date', '-d', f'{reference.isoformat()} -{offset} days', '+%Y%m%d'],
                       capture_output=True, text=True, check=True)
  return out.stdout.strip()


def most_recent_sunday(reference: date, week: int) -> date:
  return reference - timedelta(days=reference.isoweekday() % 7 + week * 7)


@pytest.mark.parametrize('reference', [
  date(2026, 8, 2),   # a Sunday: must resolve to itself, not the week before
  date(2026, 8, 3),   # Monday
  date(2026, 7, 31),  # Friday -- Dan's real invocation used 20260726
  date(2027, 1, 1),   # year boundary
])
@pytest.mark.parametrize('week', [0, 1, 3])
def test_datestamp_resolves_to_a_sunday(reference, week):
  resolved = sunday_from_shell(reference, week)
  expected = most_recent_sunday(reference, week)
  assert resolved == expected.strftime('%Y%m%d')
  assert date.fromisoformat(
    f'{resolved[:4]}-{resolved[4:6]}-{resolved[6:]}'
  ).isoweekday() == 7


def test_lookback_walks_backwards_not_forwards():
  """A missing snapshot must degrade to older data, never to a future Sunday."""
  reference = date(2026, 7, 31)
  stamps = [sunday_from_shell(reference, week) for week in range(4)]
  assert stamps == sorted(stamps, reverse=True)
  assert stamps[0] == '20260726'


def test_species_wipe_is_a_literal_path():
  """An unexpanded variable in this rm is how prod dies. It must be a literal."""
  wipes = re.findall(r'^\s*rm -rf .*$', SOURCE, re.MULTILINE)
  assert wipes, 'expected the species wipe to be present'
  for line in wipes:
    assert '$' not in line, f'variable in destructive rm: {line.strip()}'
    assert '/mnt/data/arborist/build/species/' in line


def test_fails_closed_off_the_build_host():
  """Exact full hostname before the delete -- no glob, no short-name match."""
  assert 'BUILD_HOST_FQDN=arborist-dev.lji.org' in SOURCE
  assert 'hostname -s' not in SOURCE
  assert 'arborist-dev*' not in SOURCE
  assert SOURCE.index('BUILD_HOST_FQDN') < SOURCE.index('rm -rf')


def test_build_includes_proteome_before_protein():
  """`make weekly` omits proteome; using it would reassign against stale FASTA."""
  target_line = next(line for line in SOURCE.splitlines()
                     if line.startswith('make ') and 'protein' in line)
  targets = target_line.split()
  assert 'proteome' in targets
  assert targets.index('proteome') < targets.index('protein')
  assert not re.search(r'^make weekly$', SOURCE, re.MULTILINE)
