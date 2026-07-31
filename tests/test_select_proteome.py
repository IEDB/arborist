"""Siloed, offline tests for ProteomeSelector.select() index hygiene.

Regression guard for the stale `.pepidx` bug: assign.py reuses
`proteome_5mers.pepidx` whenever it exists, so every select() path must leave
the k-mer index consistent with the current proteome.fasta -- rebuilt when the
fasta has records, absent when it is empty, and never a survivor from a prior
selection. No network, no pepmatch, no bins: fetches and preprocessing are
stubbed and the fasta is written directly.
"""

import polars as pl

from protein_tree.select_proteome import ProteomeSelector

TAXON = 424242


class _FakeBatch:
  def __init__(self, text):
    self.text = text


def _selector(species_path, proteome_list):
  sel = ProteomeSelector.__new__(ProteomeSelector)
  sel.taxon_id = TAXON
  sel.species_name = 'Testus specius'
  sel.group = 'bacterium'
  sel.peptides = pl.DataFrame({'Sequence': []})
  sel.species_path = species_path
  sel.proteome_list = proteome_list
  # Stand in for pepmatch: writing a FRESH index lets a test tell a genuine
  # rebuild apart from a surviving STALE one.
  sel._preprocess_proteome = lambda: (species_path / 'proteome_5mers.pepidx').write_text('FRESH')
  return sel


def _seed_stale(species_path):
  (species_path / 'proteome.fasta').write_text('>OLD\nMOLD\n')
  (species_path / 'proteome_5mers.pepidx').write_text('STALE')


def _index(species_path):
  p = species_path / 'proteome_5mers.pepidx'
  return p.read_text() if p.exists() else None


def _fasta(species_path):
  return (species_path / 'proteome.fasta').read_text()


def _single_candidate_list():
  return pl.DataFrame([{
    'Proteome Type': 'Reference proteome', 'Proteome ID': 'UP1',
    'Proteome Taxon ID': 123, 'Proteome Label': 'Lbl',
  }])


def test_orphan_path_rebuilds_index_and_rewrites_fasta(silo):
  sp = silo.species_path(TAXON)
  _seed_stale(sp)
  sel = _selector(sp, pl.DataFrame({'Proteome Type': []}))  # empty -> orphan path
  sel._get_batches = lambda url: [_FakeBatch('>NEW\nMNEW\n')]

  sel.select()

  assert _fasta(sp) == '>NEW\nMNEW\n'  # not appended onto the stale content
  assert _index(sp) == 'FRESH'         # index rebuilt to match new fasta


def test_orphan_path_empty_fasta_leaves_no_index(silo):
  sp = silo.species_path(TAXON)
  _seed_stale(sp)
  sel = _selector(sp, pl.DataFrame({'Proteome Type': []}))
  sel._get_batches = lambda url: []  # taxon has zero orphans

  sel.select()

  assert _fasta(sp) == ''
  assert _index(sp) is None  # stale index must not survive an empty fasta


def test_single_candidate_purges_stale_and_leftover_indexes(silo):
  sp = silo.species_path(TAXON)
  _seed_stale(sp)
  (sp / 'UPOLD_5mers.pepidx').write_text('STALE2')  # leftover candidate index
  sel = _selector(sp, _single_candidate_list())
  sel._fetch_proteome_file = lambda pid: (sp / f'{pid}.fasta').write_text('>P\nMPEP\n')
  sel._fetch_fragment_data = lambda pid: None
  sel._fetch_synonym_data = lambda pid: None

  sel.select()

  assert _fasta(sp) == '>P\nMPEP\n'
  assert _index(sp) == 'FRESH'
  assert not (sp / 'UPOLD_5mers.pepidx').exists()


def test_empty_fasta_fallback_rebuilds_from_orphans(silo):
  sp = silo.species_path(TAXON)
  _seed_stale(sp)
  sel = _selector(sp, _single_candidate_list())
  # Selected proteome fetch returns nothing -> empty fasta -> orphan fallback.
  sel._fetch_proteome_file = lambda pid: (sp / f'{pid}.fasta').write_text('')
  sel._fetch_fragment_data = lambda pid: None
  sel._fetch_synonym_data = lambda pid: None
  sel._get_batches = lambda url: [_FakeBatch('>ORPH\nMORPH\n')]

  sel.select()

  assert _fasta(sp) == '>ORPH\nMORPH\n'
  assert _index(sp) == 'FRESH'
