#!/usr/bin/env python3
"""Send the weekly Arborist digest email via the internal SMTP relay.

Verified path (arborist prod box): smtp.lji.org:25, an unauthenticated internal
relay that does NOT advertise STARTTLS. STARTTLS is therefore used only if offered,
and no credentials are sent. Callers must swallow failures -- observability email
must never break the weekly run.
"""

import smtplib
from email.message import EmailMessage

DEFAULT_SMTP_HOST = 'smtp.lji.org'
DEFAULT_SMTP_PORT = 25
DEFAULT_EMAIL_TO = 'dmarrama@lji.org'
DEFAULT_EMAIL_FROM = 'arborist@lji.org'

STATUS_GREEN = 'green'
STATUS_WARNING = 'warning'
STATUS_CRITICAL = 'critical'


def subject_for_status(status: str, week: str) -> str:
  """Traffic-light subject line; keyword surfaced when not green."""
  if status == STATUS_CRITICAL:
    return f'Arborist prod run is RED with a CRITICAL - {week}'
  if status == STATUS_WARNING:
    return f'Arborist prod run is green with WARNING(s) - {week}'
  return f'Arborist prod run is green - {week}'


def build_message(subject: str, html_body: str, *, to: str = DEFAULT_EMAIL_TO,
                  sender: str = DEFAULT_EMAIL_FROM, text_body: str = None) -> EmailMessage:
  message = EmailMessage()
  message['Subject'] = subject
  message['From'] = sender
  message['To'] = to
  message.set_content(text_body or 'This is the Arborist weekly digest; view in an HTML client.')
  message.add_alternative(html_body, subtype='html')
  return message


def send_email(message: EmailMessage, *, host: str = DEFAULT_SMTP_HOST,
               port: int = DEFAULT_SMTP_PORT, timeout: int = 30) -> None:
  """Send via the relay. Raises on failure; the caller decides how to handle it."""
  with smtplib.SMTP(host, port, timeout=timeout) as smtp:
    smtp.ehlo()
    if smtp.has_extn('starttls'):
      smtp.starttls()
      smtp.ehlo()
    smtp.send_message(message)
