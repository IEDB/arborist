#!/usr/bin/env python3
"""Urgent failure email for the weekly Arborist run (called by arborist.sh).

Unlike alert.py (an advisory digest that always exits 0), this reports that the
run itself broke: a non-zero make exit (a make child killed by the kernel OOM killer
shows up as make exit 2 plus a "Killed" line in the log), a failed pre-check, or
missing outputs. It goes to Dan only by
default. If the relay refuses the mail this exits 2 and says so on stderr, so the
failure is at least in the build log; arborist.sh still exits with the build's
own status.
"""

import argparse
import os
import platform
import re
import signal
import subprocess
import sys
from datetime import datetime
from email.message import EmailMessage
from pathlib import Path

from protein_tree import notify_email

DEFAULT_FAILURE_TO = 'dmarrama@lji.org'
LOG_TAIL_LINES = 60


def describe_exit(code: int) -> str:
  if code > 128:
    try:
      name = signal.Signals(code - 128).name
    except ValueError:
      name = f'signal {code - 128}'
    hint = ' (usually the kernel OOM killer)' if code - 128 == signal.SIGKILL else ''
    return f'exit {code} = killed by {name}{hint}'
  return f'exit {code}'


def log_tail(path: Path, lines: int = LOG_TAIL_LINES) -> str:
  try:
    text = path.read_text(errors='replace')
  except OSError as error:
    return f'(could not read {path}: {error})'
  return '\n'.join(text.splitlines()[-lines:]) or '(log is empty)'


def recent_oom_kills() -> str:
  """Best effort: kernel OOM lines, if dmesg is readable by this user."""
  try:
    output = subprocess.run(
      ['dmesg', '-T'], capture_output=True, text=True, timeout=10, check=True
    ).stdout
  except (OSError, subprocess.SubprocessError):
    return '(dmesg not readable here)'
  kills = [line for line in output.splitlines() if 'Out of memory: Killed process' in line]
  return '\n'.join(kills[-3:]) or '(no OOM kills in dmesg)'


def build_failure_message(*, stage: str, exit_code: int, log_path: Path, to: str,
                          sender: str, detail: str = '') -> EmailMessage:
  host = platform.node()
  message = EmailMessage()
  message['Subject'] = f'URGENT: Arborist weekly build FAILED on {host}: {stage} ({describe_exit(exit_code)})'
  message['From'] = sender
  message['To'] = to
  message['X-Priority'] = '1'
  message['Importance'] = 'high'
  body = [
    f'The Arborist weekly build failed at {datetime.now().astimezone():%Y-%m-%d %H:%M:%S %Z}.',
    '',
    f'Host:   {host}',
    f'Dir:    {os.getcwd()}',
    f'Stage:  {stage}',
    f'Result: {describe_exit(exit_code)}',
  ]
  tail = log_tail(log_path)
  if re.search(r'\*\*\* \[.*\] Killed', tail):
    detail = (detail + '\n' if detail else '') + \
      'make reports a child process was KILLED (SIGKILL): check the OOM lines below.'
  if detail:
    body += ['', detail]
  body += [
    '', 'Kernel OOM kills (most recent):', recent_oom_kills(),
    '', f'Last {LOG_TAIL_LINES} lines of {log_path}:', '', tail,
  ]
  message.set_content('\n'.join(body))
  return message


def main(argv=None) -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--stage', required=True, help='What was running, e.g. "make weekly".')
  parser.add_argument('--exit-code', type=int, required=True)
  parser.add_argument('--log', type=Path, required=True)
  parser.add_argument('--detail', default='')
  args = parser.parse_args(argv)

  env = os.environ
  message = build_failure_message(
    stage=args.stage, exit_code=args.exit_code, log_path=args.log,
    to=env.get('ARBORIST_FAILURE_EMAIL', DEFAULT_FAILURE_TO),
    sender=env.get('ALERT_EMAIL_FROM', notify_email.DEFAULT_EMAIL_FROM),
    detail=args.detail,
  )
  try:
    notify_email.send_email(
      message,
      host=env.get('ALERT_SMTP_HOST', notify_email.DEFAULT_SMTP_HOST),
      port=int(env.get('ALERT_SMTP_PORT', notify_email.DEFAULT_SMTP_PORT)),
    )
  except Exception as error:
    print(f'notify_failure: COULD NOT SEND the failure email ({message["Subject"]}): {error}',
          file=sys.stderr)
    return 2
  print(f'notify_failure: sent "{message["Subject"]}" to {message["To"]}', file=sys.stderr)
  return 0


if __name__ == '__main__':
  sys.exit(main())
