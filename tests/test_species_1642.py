"""Offline regression for the Listeria innocua (taxon 1642) crash.

Bus brief #124: `make protein` died with

  FileNotFoundError: build/species/1642/proteome.fasta
  at assign.py check_for_skips -> proteome.fasta.stat().st_size

Root cause: upstream Makefile stages created build/species/1642/ for a
newly-active species before proteome selection produced proteome.fasta.
`check_for_proteome` only ran on-demand selection when the species DIRECTORY
was absent, so it skipped 1642; `check_for_skips` then `.stat()`ed the missing
file and killed the whole run.

These tests freeze that exact state (a species dir with no proteome.fasta) in
tests/fixtures/species_1642/ and drive assign.py against it offline: no MySQL,
no organism rebuild, no network (ProteomeSelector is spied so on-demand
selection never touches UniProt).
"""

from types import SimpleNamespace
from pathlib import Path

import polars as pl
import pytest

import protein_tree.assign as assign

FIXTURE_BUILD = Path(__file__).parent / 'fixtures' / 'species_1642' / 'build'
TAXON = 1642
SPECIES_DIR = FIXTURE_BUILD / 'species' / str(TAXON)


class _SpySelector:
  """Stand-in for ProteomeSelector that records calls and never hits the network.

  select() deliberately produces NOTHING -> the worst-case dir-but-no-fasta
  hole. The fix must survive that (skip cleanly), so the spy models it.
  """

  instances = []

  def __init__(self, taxon_id, species_name, group, peptides, build_path):
    self.taxon_id = taxon_id
    self.selected = False
    _SpySelector.instances.append(self)

  def select(self):
    self.selected = True


@pytest.fixture
def assign_env(monkeypatch):
  data_fetcher = assign.DataFetcher(FIXTURE_BUILD)
  monkeypatch.setattr(assign, 'build_path', FIXTURE_BUILD, raising=False)
  monkeypatch.setattr(
    assign, 'active_species',
    # Active Taxa is comma-delimited text in production (e.g. "9593, 9595"); a
    # single-row fixture would otherwise infer Int64 and break the .split().
    pl.read_csv(
      FIXTURE_BUILD / 'arborist' / 'active-species.tsv', separator='\t',
      schema_overrides={'Active Taxa': pl.String},
    ),
    raising=False,
  )
  monkeypatch.setattr(assign, 'data_fetcher', data_fetcher, raising=False)
  monkeypatch.setattr(assign, 'all_peptides', data_fetcher.get_all_peptides(), raising=False)
  monkeypatch.setattr(assign, 'all_sources', data_fetcher.get_all_sources(), raising=False)
  monkeypatch.setattr(assign, 'args', SimpleNamespace(num_threads=1), raising=False)
  monkeypatch.setattr(assign, 'ProteomeSelector', _SpySelector)
  _SpySelector.instances = []
  yield
  # Keep the committed fixture pristine: drop anything a run may have created.
  for junk in SPECIES_DIR.iterdir():
    if junk.name != '.gitkeep':
      junk.unlink()


def test_fixture_reproduces_failing_state():
  assert SPECIES_DIR.is_dir()
  assert not (SPECIES_DIR / 'proteome.fasta').exists()


def test_check_for_skips_missing_proteome_returns_true(assign_env):
  # Previously raised FileNotFoundError on .stat() and killed the whole run.
  assert not (SPECIES_DIR / 'proteome.fasta').exists()
  assert assign.check_for_skips(TAXON) is True


def test_check_for_skips_empty_proteome_returns_true(assign_env):
  (SPECIES_DIR / 'proteome.fasta').touch()
  assert assign.check_for_skips(TAXON) is True


def test_check_for_proteome_triggers_selection_when_missing(assign_env):
  assign.check_for_proteome(TAXON, [TAXON], 'Listeria innocua', 'bacterium')
  assert len(_SpySelector.instances) == 1
  assert _SpySelector.instances[0].selected is True


def test_check_for_proteome_skips_selection_when_present(assign_env):
  (SPECIES_DIR / 'proteome.fasta').write_text('>fake\nMSIINFEKL\n')
  assign.check_for_proteome(TAXON, [TAXON], 'Listeria innocua', 'bacterium')
  assert _SpySelector.instances == []


def test_do_assignments_1642_no_crash_and_skips_cleanly(assign_env):
  # End-to-end: selection produces nothing (offline), so the run must skip 1642
  # cleanly and continue -- never raise.
  assign.do_assignments(TAXON)
  assert len(_SpySelector.instances) == 1, 'on-demand selection should have been attempted'
  assert not (SPECIES_DIR / 'peptide-assignments.tsv').exists(), 'should skip, not assign'
