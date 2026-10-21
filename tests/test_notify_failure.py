from protein_tree import notify_failure


def test_describe_exit_names_the_signal():
  assert notify_failure.describe_exit(2) == 'exit 2'
  assert 'SIGKILL' in notify_failure.describe_exit(137)
  assert 'OOM' in notify_failure.describe_exit(137)


def test_make_killed_child_is_called_out_with_the_log_tail(tmp_path):
  log = tmp_path / 'run.log'
  log.write_text('building\nmake: *** [Makefile:523: build/arborist/epitope-mappings.tsv] Killed\n')
  message = notify_failure.build_failure_message(
    stage='make weekly', exit_code=2, log_path=log, to='dan@example.org', sender='arborist@example.org',
  )
  body = message.get_content()
  assert message['To'] == 'dan@example.org'
  assert message['X-Priority'] == '1'
  assert message['Subject'].startswith('URGENT: Arborist weekly build FAILED')
  assert 'make weekly' in message['Subject']
  assert 'KILLED' in body
  assert 'epitope-mappings.tsv] Killed' in body


def test_unsendable_alert_exits_nonzero_and_says_so(tmp_path, monkeypatch, capsys):
  log = tmp_path / 'run.log'
  log.write_text('boom\n')

  def refuse(*args, **kwargs):
    raise ConnectionRefusedError('relay down')

  monkeypatch.setattr(notify_failure.notify_email, 'send_email', refuse)
  code = notify_failure.main(['--stage', 'make weekly', '--exit-code', '2', '--log', str(log)])
  assert code == 2
  assert 'COULD NOT SEND' in capsys.readouterr().err
