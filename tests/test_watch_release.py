"""Offline coverage for the monthly UniProt release watch (ops/uniprot-cron)."""

import importlib.util
import json
import socket
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent.parent / 'ops' / 'uniprot-cron' / 'watch-release.py'


@pytest.fixture(scope='module')
def watch():
  spec = importlib.util.spec_from_file_location('watch_release', SCRIPT)
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


RELDATE_BODY = (
  'UniProt Knowledgebase Release 2026_03 consists of:\n'
  'UniProtKB/Swiss-Prot Release 2026_03 of 05-Aug-2026\n'
  'UniProtKB/TrEMBL Release 2026_03 of 05-Aug-2026\n'
)

METALINK_BODY = (
  '<metalink xmlns="http://www.metalinker.org/" version="3.0">\n'
  '  <version>2026_03</version>\n</metalink>\n'
)


class FakeResponse:
  def __init__(self, body='', headers=None, code=200):
    self._body = body.encode()
    self.headers = headers or {}
    self.code = code

  def read(self):
    return self._body


def wire(watch, monkeypatch, *, reldate, rest='2026_03', metalink=METALINK_BODY,
         last_modified='Wed, 05 Aug 2026 20:00:00 GMT', etag='"new-etag"'):
  """Route the three signal fetches to canned responses; nothing hits the network."""
  def fake_open(request):
    url = request.full_url
    if 'reldate.txt' in url:
      if reldate is None:  # server says 304
        return FakeResponse(code=304)
      return FakeResponse(reldate, {'ETag': etag, 'Last-Modified': last_modified})
    if 'rest.uniprot.org' in url:
      return FakeResponse(headers={'x-uniprot-release': rest})
    if 'RELEASE.metalink' in url:
      return FakeResponse(metalink)
    raise AssertionError(f'unexpected fetch: {url}')

  monkeypatch.setattr(watch, '_open', fake_open)
  monkeypatch.setattr(watch, 'notify', lambda *a, **k: None)


def run(watch, state_dir, *extra):
  return watch.main(['--state-dir', str(state_dir), '--allow-any-host', '--no-notify', *extra])


def read_state(state_dir):
  return json.loads((state_dir / 'state.json').read_text())


def seed(state_dir, release='2026_02'):
  """A watch that has already adopted a baseline release."""
  state_dir.mkdir(parents=True, exist_ok=True)
  (state_dir / 'state.json').write_text(json.dumps({
    'last_seen_release': release, 'confirmed_at': '2026-06-13T00:00:00Z',
  }))


def test_first_run_adopts_current_release_without_alarm(watch, tmp_path, monkeypatch):
  """Fresh install must not report today's release as new, nor arm a refresh."""
  wire(watch, monkeypatch, reldate=RELDATE_BODY)
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 10, tzinfo=timezone.utc))

  assert run(watch, tmp_path) == 0

  state = read_state(tmp_path)
  assert state['last_seen_release'] == '2026_03'
  assert not list(tmp_path.glob('RELEASE-CONFIRMED-*'))


def test_wrong_host_fails_closed(watch, tmp_path, monkeypatch, capsys):
  monkeypatch.setattr(socket, 'gethostname', lambda: 'asahidake')
  assert watch.main(['--state-dir', str(tmp_path)]) == watch.EXIT_WRONG_HOST
  assert 'refusing' in capsys.readouterr().err
  assert not (tmp_path / 'state.json').exists()


def test_settled_release_confirms_and_writes_marker(watch, tmp_path, monkeypatch):
  seed(tmp_path)
  wire(watch, monkeypatch, reldate=RELDATE_BODY)
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 10, tzinfo=timezone.utc))

  assert run(watch, tmp_path) == 0

  state = read_state(tmp_path)
  assert state['last_seen_release'] == '2026_03'
  assert state['confirmed_at'] == '2026-08-10T00:00:00Z'
  assert state['pending_release'] is None
  assert (tmp_path / 'RELEASE-CONFIRMED-2026_03').exists()


def test_fresh_release_is_sighted_not_confirmed(watch, tmp_path, monkeypatch):
  """Inside the 48h window a release is announced, never acted on."""
  seed(tmp_path)
  wire(watch, monkeypatch, reldate=RELDATE_BODY)
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 5, 23, 0, tzinfo=timezone.utc))

  assert run(watch, tmp_path) == 0

  state = read_state(tmp_path)
  assert state['pending_release'] == '2026_03'
  assert state['last_seen_release'] == '2026_02'
  assert not (tmp_path / 'RELEASE-CONFIRMED-2026_03').exists()


def test_two_poll_fallback_confirms_without_last_modified(watch, tmp_path, monkeypatch):
  """No usable Last-Modified: the stored first-sighting timestamp settles it."""
  seed(tmp_path)
  wire(watch, monkeypatch, reldate=RELDATE_BODY, last_modified='')
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 5, 23, 0, tzinfo=timezone.utc))
  assert run(watch, tmp_path) == 0
  assert read_state(tmp_path)['pending_release'] == '2026_03'

  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 8, 23, 0, tzinfo=timezone.utc))
  assert run(watch, tmp_path) == 0
  state = read_state(tmp_path)
  assert state['last_seen_release'] == '2026_03'
  assert state['first_seen_at'] == '2026-08-05T23:00:00Z'


def test_source_disagreement_never_confirms(watch, tmp_path, monkeypatch):
  """A mirror mid-sync (REST still on the old release) must not trigger anything."""
  seed(tmp_path)
  wire(watch, monkeypatch, reldate=RELDATE_BODY, rest='2026_02')
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 10, tzinfo=timezone.utc))

  assert run(watch, tmp_path) == 0

  state = read_state(tmp_path)
  assert state['last_seen_release'] == '2026_02'
  assert state['pending_release'] is None
  assert not list(tmp_path.glob('RELEASE-CONFIRMED-*'))


def test_known_release_is_a_no_op(watch, tmp_path, monkeypatch):
  (tmp_path / 'state.json').write_text(json.dumps({
    'last_seen_release': '2026_03', 'reldate_etag': '"old"',
    'confirmed_at': '2026-08-07T00:00:00Z',
  }))
  wire(watch, monkeypatch, reldate=None)  # 304 on the stored ETag
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 20, tzinfo=timezone.utc))

  assert run(watch, tmp_path) == 0

  state = read_state(tmp_path)
  assert state['last_seen_release'] == '2026_03'
  assert state['last_poll_at'] == '2026-08-20T00:00:00Z'
  assert not list(tmp_path.glob('RELEASE-CONFIRMED-*'))


def test_fetch_failure_reports_and_exits_nonzero(watch, tmp_path, monkeypatch):
  """Silence must mean broken: a failed poll records the error and exits 70."""
  def boom(request):
    raise watch.FetchError('ftp.uniprot.org: timed out')

  monkeypatch.setattr(watch, '_open', boom)
  monkeypatch.setattr(watch, 'notify', lambda *a, **k: None)

  assert run(watch, tmp_path) == watch.EXIT_FETCH_FAILED
  assert 'timed out' in read_state(tmp_path)['last_poll_error']


def test_dry_run_writes_nothing(watch, tmp_path, monkeypatch):
  state_dir = tmp_path / 'state'
  seed(state_dir)
  before = (state_dir / 'state.json').read_text()
  wire(watch, monkeypatch, reldate=RELDATE_BODY)
  monkeypatch.setattr(watch, 'utcnow',
                      lambda: datetime(2026, 8, 10, tzinfo=timezone.utc))

  assert run(watch, state_dir, '--dry-run') == 0
  assert (state_dir / 'state.json').read_text() == before
  assert not list(state_dir.glob('RELEASE-CONFIRMED-*'))


def test_deployment_date_is_not_a_signal(watch):
  """x-api-deployment-date redeploys independently; it must not appear anywhere."""
  source = SCRIPT.read_text()
  assert 'x-uniprot-release' in source
  assert 'x-api-deployment-date' not in source.split('"""', 2)[2]


@pytest.mark.parametrize('stored,released,expected', [
  ({'last_seen_release': '2026_03'}, None, 'known'),
  ({}, datetime(2026, 8, 5, tzinfo=timezone.utc), 'confirmed'),
  ({}, datetime(2026, 8, 9, 12, tzinfo=timezone.utc), 'sighted'),
  ({'pending_release': '2026_03', 'pending_first_seen_at': '2026-08-09T00:00:00Z'},
   None, 'settling'),
  ({'pending_release': '2026_03', 'pending_first_seen_at': '2026-08-07T00:00:00Z'},
   None, 'confirmed'),
])
def test_decide_matrix(watch, stored, released, expected):
  now = datetime(2026, 8, 10, tzinfo=timezone.utc)
  assert watch.decide(stored, '2026_03', now, 48, released) == expected


def test_release_mail_does_not_go_to_the_weekly_digest_list(watch):
  """Randi is on the prod weekly stats, not dev UniProt plumbing."""
  assert watch.WATCH_EMAIL_TO == 'dmarrama@lji.org'
  assert 'rvita' not in watch.WATCH_EMAIL_TO
  code = [line for line in SCRIPT.read_text().splitlines()
          if not line.lstrip().startswith('#')]
  assert not any('DEFAULT_EMAIL_TO' in line for line in code)
