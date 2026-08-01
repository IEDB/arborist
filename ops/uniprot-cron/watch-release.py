#!/usr/bin/env python3
"""Detect a settled UniProt release and say so. Builds nothing, promotes nothing.

Three independent signals must agree before a release counts:
  1. FTP reldate.txt          (conditional GET on the stored ETag)
  2. rest.uniprot.org header  x-uniprot-release  (what select_proteome.py actually reads)
  3. reference_proteomes/RELEASE.metalink <version>
x-api-deployment-date is deliberately ignored: it moves without a data release.

A new string is only *confirmed* after it has held for --stability-hours (default 48),
which keeps a half-mirrored release from triggering anything. On confirmation the
script writes RELEASE-CONFIRMED-<release> in the state dir; that marker is the
handoff to the (not yet built) refresh runner. Every terminal outcome notifies --
silence must mean the timer is broken, not that nothing happened.
"""

import argparse
import json
import os
import re
import socket
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

BUILD_HOST = 'arborist-dev'
STATE_DIR = Path('/var/lib/arborist-uniprot-watch')
STABILITY_HOURS = 48
HTTP_TIMEOUT = 30
USER_AGENT = 'arborist-uniprot-watch/1 (+dmarrama@lji.org)'

# Dev release plumbing, not the weekly digest -- deliberately narrower than
# notify_email.DEFAULT_EMAIL_TO. Randi gets the prod weekly stats, not this.
WATCH_EMAIL_TO = 'dmarrama@lji.org'

RELDATE_URL = ('https://ftp.uniprot.org/pub/databases/uniprot/current_release/'
               'knowledgebase/complete/reldate.txt')
REST_URL = 'https://rest.uniprot.org/uniprotkb/P05067.json'
METALINK_URL = ('https://ftp.uniprot.org/pub/databases/uniprot/current_release/'
                'knowledgebase/reference_proteomes/RELEASE.metalink')

RELEASE_RE = re.compile(r'Release (\d{4}_\d{2})')
SWISSPROT_DATE_RE = re.compile(r'Swiss-Prot Release \d{4}_\d{2} of (\S+)')
METALINK_VERSION_RE = re.compile(r'<version>\s*(\d{4}_\d{2})\s*</version>')

EXIT_WRONG_HOST = 78
EXIT_FETCH_FAILED = 70


def utcnow() -> datetime:
  return datetime.now(timezone.utc)


def iso(moment: datetime) -> str:
  return moment.strftime('%Y-%m-%dT%H:%M:%SZ')


def parse_iso(value):
  if not value:
    return None
  return datetime.strptime(value, '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)


def parse_http_date(value):
  """RFC 2822 Last-Modified -> aware datetime; None when absent or malformed."""
  if not value:
    return None
  try:
    moment = parsedate_to_datetime(value)
  except (TypeError, ValueError):
    return None
  return moment if moment.tzinfo else moment.replace(tzinfo=timezone.utc)


class FetchError(RuntimeError):
  pass


def _open(request):
  try:
    return urllib.request.urlopen(request, timeout=HTTP_TIMEOUT)
  except urllib.error.HTTPError as error:
    if error.code == 304:
      return error
    raise FetchError(f'{request.full_url}: HTTP {error.code}') from error
  except (urllib.error.URLError, socket.timeout, OSError) as error:
    raise FetchError(f'{request.full_url}: {error}') from error


def fetch_reldate(etag: str) -> dict:
  """Conditional GET; 304 means the release string is unchanged since last poll."""
  headers = {'User-Agent': USER_AGENT}
  if etag:
    headers['If-None-Match'] = etag
  response = _open(urllib.request.Request(RELDATE_URL, headers=headers))
  if getattr(response, 'code', None) == 304:
    return {'not_modified': True}
  body = response.read().decode('utf-8', 'replace')
  release = RELEASE_RE.search(body)
  if not release:
    raise FetchError(f'reldate.txt has no release string: {body!r}')
  date = SWISSPROT_DATE_RE.search(body)
  return {
    'not_modified': False,
    'release': release.group(1),
    'release_date': date.group(1) if date else None,
    'etag': response.headers.get('ETag', ''),
    'last_modified': response.headers.get('Last-Modified', ''),
  }


def fetch_rest_release() -> str:
  request = urllib.request.Request(REST_URL, headers={'User-Agent': USER_AGENT}, method='HEAD')
  release = _open(request).headers.get('x-uniprot-release')
  if not release:
    raise FetchError('rest.uniprot.org sent no x-uniprot-release header')
  return release.strip()


def fetch_metalink_version() -> str:
  request = urllib.request.Request(METALINK_URL, headers={'User-Agent': USER_AGENT})
  body = _open(request).read().decode('utf-8', 'replace')
  version = METALINK_VERSION_RE.search(body)
  if not version:
    raise FetchError('RELEASE.metalink has no <version>')
  return version.group(1)


def load_state(path: Path) -> dict:
  if not path.exists():
    return {}
  return json.loads(path.read_text())


def save_state(path: Path, state: dict) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  temp = path.with_suffix('.json.tmp')
  temp.write_text(json.dumps(state, indent=2, sort_keys=True) + '\n')
  temp.replace(path)


def send_email(subject: str, body: str, repo: Path) -> str:
  """Reuse the weekly digest relay. Never fatal -- returns a status string."""
  sys.path.insert(0, str(repo / 'src' / 'protein'))
  try:
    from protein_tree import notify_email
  except ImportError as error:
    return f'email skipped (no protein_tree at {repo}: {error})'
  to = os.environ.get('WATCH_EMAIL_TO', WATCH_EMAIL_TO)
  html = '<pre>' + body.replace('&', '&amp;').replace('<', '&lt;') + '</pre>'
  try:
    message = notify_email.build_message(subject, html, to=to, text_body=body)
    notify_email.send_email(message)
  except Exception as error:  # observability must not break the watcher
    return f'email FAILED: {error}'
  return f'email sent to {to}'


def send_ntfy(subject: str, body: str, priority: str) -> str:
  """Optional push. Configured by NTFY_URL (+NTFY_TOPIC); silently off when unset."""
  url = os.environ.get('NTFY_URL')
  topic = os.environ.get('NTFY_TOPIC', 'jobs')
  if not url:
    return 'ntfy skipped (NTFY_URL unset)'
  command = ['curl', '-fsS', '-m', '15',
             '-H', f'Title: {subject}', '-H', f'Priority: {priority}', '-H', 'Tags: dna',
             '--data-binary', body, f'{url.rstrip("/")}/{topic}']
  user = os.environ.get('NTFY_USER')
  password = os.environ.get('NTFY_PASSWORD')
  if user and password:
    command[1:1] = ['-u', f'{user}:{password}']
  result = subprocess.run(command, capture_output=True, text=True)
  return 'ntfy sent' if result.returncode == 0 else f'ntfy FAILED: {result.stderr.strip()}'


def notify(subject: str, body: str, args, priority: str = 'default') -> None:
  log(subject)
  log(body)
  if args.no_notify:
    log('notifications suppressed (--no-notify)')
    return
  log(send_email(subject, body, args.repo))
  log(send_ntfy(subject, body, priority))


def log(message: str) -> None:
  print(f'[{iso(utcnow())}] {message}', flush=True)


def decide(state: dict, release: str, now: datetime, stability_hours: int,
           released_at: datetime = None) -> str:
  """known | sighted | settling | confirmed -- the outcomes for an agreed release.

  A weekly timer cannot rely on two polls 48h apart, so the release's own
  Last-Modified settles it when available; the two-poll path is the fallback.
  """
  if release == state.get('last_seen_release'):
    return 'known'
  window = timedelta(hours=stability_hours)
  if released_at and now - released_at >= window:
    return 'confirmed'
  if state.get('pending_release') != release:
    return 'sighted'
  first_seen = parse_iso(state.get('pending_first_seen_at'))
  if first_seen and now - first_seen >= window:
    return 'confirmed'
  return 'settling'


def parse_args(argv):
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-s', '--state-dir', type=Path, default=STATE_DIR)
  parser.add_argument('-r', '--repo', type=Path,
                      default=Path(os.environ.get('ARBORIST_REPO', '/home/dmarrama/lji/arborist')),
                      help='arborist checkout, for the notify_email relay')
  parser.add_argument('--stability-hours', type=int, default=STABILITY_HOURS,
                      help='how long a new release string must hold before it counts')
  parser.add_argument('--no-notify', action='store_true', help='log only; no email, no ntfy')
  parser.add_argument('--dry-run', action='store_true', help='do not write state or markers')
  parser.add_argument('--allow-any-host', action='store_true',
                      help='bypass the arborist-dev gate (testing only)')
  return parser.parse_args(argv)


def main(argv=None) -> int:
  args = parse_args(argv)

  host = socket.gethostname().split('.')[0]
  if host != BUILD_HOST and not args.allow_any_host:
    print(f'refusing: not {BUILD_HOST} (this is {host})', file=sys.stderr)
    return EXIT_WRONG_HOST

  state_path = args.state_dir / 'state.json'
  state = load_state(state_path)
  fresh_state = not state
  now = utcnow()

  released_at = None
  try:
    reldate = fetch_reldate(state.get('reldate_etag', ''))
    if reldate['not_modified']:
      reldate_release = state.get('last_seen_release')
      log(f'reldate.txt 304 not-modified (release {reldate_release})')
      if not reldate_release:
        raise FetchError('reldate.txt returned 304 but no release is stored yet')
    else:
      reldate_release = reldate['release']
      state['reldate_etag'] = reldate['etag']
      state['reldate_last_modified'] = reldate['last_modified']
      state['last_seen_release_date'] = reldate['release_date']
    released_at = parse_http_date(state.get('reldate_last_modified', ''))
    rest_release = fetch_rest_release()
    metalink_release = fetch_metalink_version()
  except FetchError as error:
    state['last_poll_at'] = iso(now)
    state['last_poll_error'] = str(error)
    if not args.dry_run:
      save_state(state_path, state)
    notify('Arborist UniProt watch FAILED',
           f'Release poll failed, no release decision made.\n\n{error}\n\n'
           f'Host: {host}\nState: {state_path}\n'
           'The scheduled check did not complete -- fix the fetch or the box.',
           args, priority='high')
    return EXIT_FETCH_FAILED

  state['last_poll_at'] = iso(now)
  state.pop('last_poll_error', None)
  state['rest_release_header'] = rest_release
  state['metalink_version'] = metalink_release
  log(f'reldate={reldate_release} rest={rest_release} metalink={metalink_release}')

  sources = {reldate_release, rest_release, metalink_release}
  if len(sources) != 1:
    state['pending_release'] = None
    state['pending_first_seen_at'] = None
    if not args.dry_run:
      save_state(state_path, state)
    notify(f'Arborist UniProt watch: sources disagree ({", ".join(sorted(sources))})',
           f'reldate.txt={reldate_release}\nrest x-uniprot-release={rest_release}\n'
           f'RELEASE.metalink={metalink_release}\n\n'
           'Likely a mirror mid-sync. No refresh triggered; the next poll re-checks.',
           args)
    return 0

  release = sources.pop()
  verdict = decide(state, release, now, args.stability_hours, released_at)

  if fresh_state:
    # First run on a new box: adopt whatever UniProt is on now. Announcing the
    # current release as "new" would be a lie and would arm the refresh path.
    state['last_seen_release'] = release
    state['first_seen_at'] = iso(released_at or now)
    state['confirmed_at'] = iso(now)
    state['last_notified_at'] = iso(now)
    state['pending_release'] = None
    state['pending_first_seen_at'] = None
    if not args.dry_run:
      save_state(state_path, state)
    notify(f'Arborist UniProt watch: baseline seeded at {release}',
           f'Watch installed on {host}; adopting current release {release} as the baseline.\n'
           'No refresh triggered. The next distinct, settled release will be reported.',
           args, priority='low')
    return 0

  if verdict == 'known':
    log(f'no new release (still {release})')
    if not args.dry_run:
      save_state(state_path, state)
    notify(f'Arborist UniProt watch: no new release ({release})',
           f'UniProt is still on {release}; last confirmed {state.get("confirmed_at")}.\n'
           'Nothing to refresh. This mail proves the watch ran.',
           args, priority='low')
    return 0

  if verdict == 'sighted':
    state['pending_release'] = release
    state['pending_first_seen_at'] = iso(now)
    state['last_notified_at'] = iso(now)
    if not args.dry_run:
      save_state(state_path, state)
    notify(f'Arborist UniProt watch: {release} sighted (unconfirmed)',
           f'All three sources report {release} (was {state.get("last_seen_release")}).\n'
           f'Holding for {args.stability_hours}h of stability before confirming.\n'
           'Nothing was built. Re-run the watch after the window to confirm.',
           args)
    return 0

  if verdict == 'settling':
    first_seen = state.get('pending_first_seen_at')
    log(f'{release} still settling since {first_seen}; no notification')
    if not args.dry_run:
      save_state(state_path, state)
    return 0

  previous = state.get('last_seen_release')
  state['last_seen_release'] = release
  state['first_seen_at'] = state.get('pending_first_seen_at') or iso(released_at or now)
  state['confirmed_at'] = iso(now)
  state['last_notified_at'] = iso(now)
  state['pending_release'] = None
  state['pending_first_seen_at'] = None
  marker = args.state_dir / f'RELEASE-CONFIRMED-{release}'
  if not args.dry_run:
    save_state(state_path, state)
    marker.write_text(f'confirmed_at={iso(now)}\nprevious={previous}\n')
  age = f'released {iso(released_at)}' if released_at else 'release date unknown'
  notify(f'Arborist UniProt watch: {release} CONFIRMED',
         f'UniProt {release} is settled ({age}; previous {previous}).\n'
         f'Marker: {marker}\n\n'
         'This watcher builds nothing. Run the refresh on arborist-dev:\n\n'
         '  sudo nohup /opt/arborist-uniprot/bin/refresh-run.sh \\\n'
         '    > /var/log/arborist-uniprot/refresh.log 2>&1 &\n\n'
         'It takes hours and promotes nothing. Watch it with\n'
         '  tail -f /var/log/arborist-uniprot/refresh.log',
         args, priority='high')
  return 0


if __name__ == '__main__':
  sys.exit(main())
