import smtplib

from protein_tree import notify_email
from protein_tree.notify_email import STATUS_GREEN, STATUS_WARNING, STATUS_CRITICAL


def test_subject_wording():
  assert notify_email.subject_for_status(STATUS_GREEN, '2026-W29') == \
    'Arborist prod run is green - 2026-W29'
  assert notify_email.subject_for_status(STATUS_WARNING, '2026-W29') == \
    'Arborist prod run is green with WARNING(s) - 2026-W29'
  assert notify_email.subject_for_status(STATUS_CRITICAL, '2026-W29') == \
    'Arborist prod run is RED with a CRITICAL - 2026-W29'


def test_build_message_headers_and_html():
  msg = notify_email.build_message('Subj', '<p>hi</p>', to='a@b.org', sender='c@d.org')
  assert msg['To'] == 'a@b.org'
  assert msg['From'] == 'c@d.org'
  assert msg['Subject'] == 'Subj'
  html_parts = [p.get_content() for p in msg.walk() if p.get_content_type() == 'text/html']
  assert html_parts and '<p>hi</p>' in html_parts[0]


class _FakeSMTP:
  instances = []
  starttls_available = False

  def __init__(self, host, port, timeout=None):
    self.host = host
    self.port = port
    self.events = []
    self.sent = None
    _FakeSMTP.instances.append(self)

  def __enter__(self):
    return self

  def __exit__(self, *exc):
    return False

  def ehlo(self):
    self.events.append('ehlo')

  def has_extn(self, name):
    return name == 'starttls' and _FakeSMTP.starttls_available

  def starttls(self):
    self.events.append('starttls')

  def send_message(self, message):
    self.events.append('send')
    self.sent = message


def test_send_email_plain_when_no_starttls(monkeypatch):
  _FakeSMTP.instances = []
  _FakeSMTP.starttls_available = False
  monkeypatch.setattr(smtplib, 'SMTP', _FakeSMTP)

  message = notify_email.build_message('s', '<p>x</p>')
  notify_email.send_email(message, host='smtp.lji.org', port=25)

  inst = _FakeSMTP.instances[-1]
  assert (inst.host, inst.port) == ('smtp.lji.org', 25)
  assert 'send' in inst.events
  assert 'starttls' not in inst.events  # relay offers none -> plain (matches prod)


def test_send_email_uses_starttls_when_offered(monkeypatch):
  _FakeSMTP.instances = []
  _FakeSMTP.starttls_available = True
  monkeypatch.setattr(smtplib, 'SMTP', _FakeSMTP)

  notify_email.send_email(notify_email.build_message('s', '<p>x</p>'))
  assert 'starttls' in _FakeSMTP.instances[-1].events
