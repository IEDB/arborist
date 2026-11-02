"""Offline tests: a UniProt proteome fetch that fails for good must be SKIPPED,
never crash the build (2026-11-01: 8x HTTP 500 on one candidate proteome left no
FASTA, and PEPMatch's Rust preprocess panicked on the missing file).

No network and no pepmatch: the HTTP session is scripted, time.sleep is a no-op
and the match count is stubbed where the real Matcher would be needed.
"""

import types

import polars as pl
import pytest
from requests.exceptions import ChunkedEncodingError, HTTPError

from protein_tree import alert, skipped_proteomes
from protein_tree import select_proteome as sp_mod
from protein_tree.select_proteome import (
  MAX_HTTP_ATTEMPTS, ProteomeFetchError, ProteomeSelector, proteome_complete,
)

TAXON = 2485925
FASTA = '>P1\nMPEPTIDEAAAA\n'


class _Resp:
  """requests.Response stand-in. `chunks` may end with an Exception to raise mid-stream."""

  def __init__(self, status=200, chunks=(FASTA.encode(),)):
    self.status_code = status
    self.chunks = list(chunks)

  def __enter__(self):
    return self

  def __exit__(self, *exc):
    return False

  def raise_for_status(self):
    if self.status_code >= 400:
      raise HTTPError(f'{self.status_code} Server Error', response=self)

  def iter_content(self, chunk_size=None):
    for chunk in self.chunks:
      if isinstance(chunk, Exception):
        raise chunk
      yield chunk


class _Session:
  """Serves scripted responses in order; the last one repeats once the script runs out."""

  def __init__(self, *script):
    self.script = list(script)
    self.calls = []

  def get(self, url, **kwargs):
    self.calls.append(url)
    index = min(len(self.calls) - 1, len(self.script) - 1)
    return self.script[index]


@pytest.fixture(autouse=True)
def _no_sleep(monkeypatch):
  monkeypatch.setattr(sp_mod.time, 'sleep', lambda s: None)


def _selector(species_path, proteome_list, session=None):
  sel = ProteomeSelector.__new__(ProteomeSelector)
  sel.taxon_id = TAXON
  sel.species_name = 'Oscillospiraceae bacterium'
  sel.group = 'bacterium'
  sel.peptides = pl.DataFrame({'Sequence': []})
  sel.species_path = species_path
  sel.proteome_list = proteome_list
  sel.session = session
  sel._preprocess_proteome = lambda: (species_path / 'proteome_5mers.pepidx').write_text('FRESH')
  sel._fetch_fragment_data = lambda pid: None
  sel._fetch_synonym_data = lambda pid: (species_path / 'synonym-data.json').write_text('{}')
  return sel


def _candidates(*ids):
  return pl.DataFrame([{
    'Proteome Type': 'Non Reference proteome', 'Proteome ID': pid,
    'Proteome Taxon ID': TAXON, 'Proteome Label': 'Lbl',
  } for pid in ids])


def _fetch_stub(sel, failing=(), counts=None):
  """Replace the HTTP fetch: `failing` ids raise like exhausted retries; others write a FASTA."""
  def fetch(pid):
    if pid in failing:
      raise ProteomeFetchError(pid, MAX_HTTP_ATTEMPTS, HTTPError('500 Server Error'))
    (sel.species_path / f'{pid}.fasta').write_text(FASTA)
  sel._fetch_proteome_file = fetch
  sel._get_match_count = lambda linear, disc, pid: (counts or {}).get(pid, 0)


def _dir_names(path):
  return sorted(p.name for p in path.iterdir())


# ---- _fetch_proteome_file: exhausted retries, partial streams, recovery ----

def test_all_retries_fail_raises_and_leaves_no_file(silo):
  path = silo.species_path(TAXON)
  session = _Session(_Resp(status=500))
  sel = _selector(path, None, session)

  with pytest.raises(ProteomeFetchError) as raised:
    sel._fetch_proteome_file('UP001055208')

  assert raised.value.proteome_id == 'UP001055208'
  assert raised.value.attempts == MAX_HTTP_ATTEMPTS
  assert len(session.calls) == MAX_HTTP_ATTEMPTS
  assert _dir_names(path) == []  # neither UP001055208.fasta nor a .part leftover


def test_partial_stream_never_leaves_a_fasta(silo):
  path = silo.species_path(TAXON)
  # Every attempt delivers a first chunk, then the connection drops mid-stream.
  session = _Session(_Resp(chunks=[b'>P1\nMPEP', ChunkedEncodingError('connection broken')]))
  sel = _selector(path, None, session)

  with pytest.raises(ProteomeFetchError):
    sel._fetch_proteome_file('UP1')

  assert _dir_names(path) == []  # a truncated FASTA must not look like a fetched one


def test_500_then_success_writes_the_fasta(silo):
  path = silo.species_path(TAXON)
  session = _Session(_Resp(status=500), _Resp(status=500), _Resp())
  sel = _selector(path, None, session)

  sel._fetch_proteome_file('UP1')

  assert len(session.calls) == 3
  assert (path / 'UP1.fasta').read_text() == FASTA
  assert _dir_names(path) == ['UP1.fasta']  # the .part file was renamed away


def test_non_retryable_http_error_fails_without_retrying(silo):
  path = silo.species_path(TAXON)
  session = _Session(_Resp(status=404))
  sel = _selector(path, None, session)

  with pytest.raises(ProteomeFetchError) as raised:
    sel._fetch_proteome_file('UP1')

  assert raised.value.attempts == 1
  assert len(session.calls) == 1


# ---- _proteome_tiebreak / select(): drop the candidate, keep the build going ----

def test_only_one_candidate_fails_best_remaining_wins_and_is_recorded(silo):
  path = silo.species_path(TAXON)
  sel = _selector(path, _candidates('UPA', 'UPB', 'UPC'))
  _fetch_stub(sel, failing={'UPC'}, counts={'UPA': 1, 'UPB': 5, 'UPC': 99})

  sel.select()

  assert (path / 'proteome.fasta').read_text() == FASTA
  assert (path / 'species-data.tsv').exists()
  meta = pl.read_csv(path / 'species-data.tsv', separator='\t')
  assert meta['Proteome ID'][0] == 'UPB'  # the failed UPC would have won; it is dropped, not scored
  assert not (path / 'UPC.fasta').exists()
  rows = skipped_proteomes.load_skips(silo.build)
  assert [(r['Proteome ID'], r['Outcome']) for r in rows] == [('UPC', skipped_proteomes.OUTCOME_CANDIDATE_DROPPED)]
  assert 'after 8 attempts' in rows[0]['Error']


def test_skip_is_logged_in_the_documented_format(silo, capsys):
  path = silo.species_path(TAXON)
  sel = _selector(path, _candidates('UPA', 'UPB'))
  _fetch_stub(sel, failing={'UPB'}, counts={'UPA': 1})

  sel.select()

  err = capsys.readouterr().err
  assert (f'[skip] proteome UPB for taxon {TAXON} (Oscillospiraceae bacterium): '
          'fetch failed after 8 attempts: HTTPError: 500 Server Error') in err


def test_all_candidates_fail_skips_the_taxon_without_raising(silo):
  path = silo.species_path(TAXON)
  (path / 'proteome_5mers.pepidx').write_text('STALE')
  sel = _selector(path, _candidates('UPA', 'UPB'))
  _fetch_stub(sel, failing={'UPA', 'UPB'})

  sel.select()  # must not raise

  # Same state the existing 'Proteome is empty, skipping' branch expects.
  assert (path / 'proteome.fasta').read_text() == ''
  assert not (path / 'species-data.tsv').exists()
  assert not list(path.glob('*.pepidx'))
  assert not list(path.glob('UP*.fasta*'))
  assert silo.assign.check_for_skips(TAXON) is True
  outcomes = [(r['Proteome ID'], r['Outcome']) for r in skipped_proteomes.load_skips(silo.build)]
  assert outcomes == [
    ('UPA', skipped_proteomes.OUTCOME_CANDIDATE_DROPPED),
    ('UPB', skipped_proteomes.OUTCOME_CANDIDATE_DROPPED),
    ('', skipped_proteomes.OUTCOME_TAXON_SKIPPED),
  ]


def test_single_candidate_fetch_failure_skips_the_taxon(silo):
  path = silo.species_path(TAXON)
  sel = _selector(path, _candidates('UPONLY'))
  _fetch_stub(sel, failing={'UPONLY'})

  sel.select()

  assert (path / 'proteome.fasta').read_text() == ''
  assert not (path / 'species-data.tsv').exists()
  rows = skipped_proteomes.load_skips(silo.build)
  assert [(r['Proteome ID'], r['Outcome']) for r in rows] == [('UPONLY', skipped_proteomes.OUTCOME_TAXON_SKIPPED)]


def test_candidate_list_fetch_failure_skips_the_taxon(silo, monkeypatch):
  path = silo.species_path(TAXON)
  monkeypatch.setattr(ProteomeSelector, '_get_candidate_proteomes',
                      lambda self: (_ for _ in ()).throw(HTTPError('503 Server Error')))

  sel = ProteomeSelector(TAXON, 'Oscillospiraceae bacterium', 'bacterium', pl.DataFrame({'Sequence': []}), silo.build)
  sel.select()  # constructor and select() must not raise

  assert (path / 'proteome.fasta').read_text() == ''
  assert [r['Outcome'] for r in skipped_proteomes.load_skips(silo.build)] == [skipped_proteomes.OUTCOME_TAXON_SKIPPED]


def test_missing_fasta_never_reaches_pepmatch(silo, monkeypatch):
  path = silo.species_path(TAXON)
  sel = _selector(path, None)

  def matcher(**kwargs):
    raise AssertionError('PEPMatch must not be called for a FASTA that is missing or empty')

  monkeypatch.setattr(sp_mod, 'Matcher', matcher)
  assert sel._get_match_count(['MPEPTIDE'], [], 'UPNOFILE') == 0
  (path / 'UPEMPTY.fasta').write_text('')
  assert sel._get_match_count(['MPEPTIDE'], [], 'UPEMPTY') == 0


# ---- a failure is never cached ----

def test_skipped_taxon_is_retried_and_recovers(silo):
  path = silo.species_path(TAXON)
  sel = _selector(path, _candidates('UPA'))
  _fetch_stub(sel, failing={'UPA'})
  sel.select()
  assert skipped_proteomes.has_skips(path)
  assert proteome_complete(TAXON, silo.build) is False

  _fetch_stub(sel)  # UniProt is back
  sel.select()

  assert (path / 'proteome.fasta').read_text() == FASTA
  assert (path / 'species-data.tsv').exists()
  assert not skipped_proteomes.has_skips(path)  # record cleared by the clean run


def test_dropped_candidate_keeps_the_species_unfinished_for_resume(silo):
  path = silo.species_path(TAXON)
  sel = _selector(path, _candidates('UPA', 'UPB'))
  _fetch_stub(sel, failing={'UPB'}, counts={'UPA': 1})
  sel.select()
  (path / 'fragment-data.json').write_text('{}')
  # Every output a finished species has is on disk...
  assert (path / 'species-data.tsv').exists() and (path / 'synonym-data.json').exists()
  assert (path / 'proteome_5mers.pepidx').exists()
  # ...yet a species that dropped a candidate must be redone, not resumed as complete.
  assert proteome_complete(TAXON, silo.build) is False

  skipped_proteomes.clear_skips(path)
  assert proteome_complete(TAXON, silo.build) is True


def test_assign_reselects_a_species_with_a_skip_record_even_if_its_fasta_exists(silo, monkeypatch):
  path = silo.species_path(TAXON)
  (path / 'proteome.fasta').write_text(FASTA)
  skipped_proteomes.record_skip(path, TAXON, 'X', 'UPC', skipped_proteomes.OUTCOME_CANDIDATE_DROPPED, 'boom')
  selected = []

  class FakeSelector:
    def __init__(self, *args):
      pass

    def select(self):
      selected.append(True)

  class FakeFetcher:
    def __init__(self, build_path):
      pass

    def get_peptides_for_species(self, all_peptides, taxa):
      return None

  monkeypatch.setattr(silo.assign, 'ProteomeSelector', FakeSelector)
  monkeypatch.setattr(silo.assign, 'DataFetcher', FakeFetcher)
  monkeypatch.setattr(silo.assign, 'all_peptides', None, raising=False)

  silo.assign.check_for_proteome(TAXON, [TAXON], 'X', 'bacterium')
  assert selected == [True]

  skipped_proteomes.clear_skips(path)  # no record and a real FASTA: nothing to redo
  silo.assign.check_for_proteome(TAXON, [TAXON], 'X', 'bacterium')
  assert selected == [True]


# ---- visibility: end-of-run summary and the weekly digest ----

def test_summary_and_digest_list_every_skip(silo, tmp_path):
  path = silo.species_path(TAXON)
  skipped_proteomes.record_skip(path, TAXON, 'Oscillospiraceae bacterium', 'UP001055208',
                                skipped_proteomes.OUTCOME_CANDIDATE_DROPPED,
                                'fetch failed after 8 attempts: HTTPError: 500 <b>boom</b>')
  skips = skipped_proteomes.load_skips(silo.build)

  summary = skipped_proteomes.format_skips(skips)
  assert 'UP001055208' in summary and str(TAXON) in summary

  sheet = pl.DataFrame(
    {'Species ID': ['TOTAL'], 'Species Name': [''], 'Metric': ['epitopes_assigned'], '2026-W01': ['10'], '2026-W02': ['10']})
  html = alert.render_html(sheet, [], ['2026-W01', '2026-W02'], 'warning', '2026-W02', skips)
  assert 'UP001055208' in html
  assert 'Oscillospiraceae bacterium' in html
  assert '&lt;b&gt;boom' in html  # error text is escaped, not injected into the email


def test_run_with_skips_is_a_warning_not_green(silo, tmp_path):
  skipped_proteomes.record_skip(silo.species_path(TAXON), TAXON, 'X', 'UPC',
                                skipped_proteomes.OUTCOME_TAXON_SKIPPED, 'boom')
  stats = tmp_path / 'stats.tsv'
  metrics = alert.OVERVIEW_METRICS
  pl.DataFrame({'Species ID': ['TOTAL'] * len(metrics), 'Species Name': [''] * len(metrics), 'Metric': metrics,
                '2026-W01': ['10'] * len(metrics), '2026-W02': ['10'] * len(metrics)}).write_csv(stats, separator='\t')
  args = types.SimpleNamespace(stats_file=str(stats), build_path=str(silo.build), html_out=str(tmp_path / 'alert.html'))

  alert.run(args, env={'ALERT_EMAIL': '0'})

  html = (tmp_path / 'alert.html').read_text()
  assert 'WARNING' in html and 'UPC' in html
