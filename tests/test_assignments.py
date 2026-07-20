import os
import sys
import shutil
import subprocess
from pathlib import Path

import pytest
import polars as pl

TESTS_DIR = Path(__file__).parent
SRC_PROTEIN = TESTS_DIR.parent / 'src' / 'protein'
ASSIGN_PY = SRC_PROTEIN / 'protein_tree' / 'assign.py'
BIN = TESTS_DIR.parent / 'bin'

species = [108, 9593, 9606, 3052236]


def _aligners_available():
  tools = ('makeblastdb', 'blastp', 'mmseqs')
  return all((BIN / t).exists() or shutil.which(t) for t in tools)


# End-to-end: runs a real assign.py against the committed golden fixtures in
# tests/build and diffs the output against each species' expected-assignments.tsv.
# Needs the external aligners (blastp/mmseqs2) + ARC, so it is marked `e2e`:
# bare `pytest` deselects it; run `pytest -m e2e` on arborist-dev. The skipif is
# belt-and-suspenders so a stray e2e selection without bins skips, not fails.
pytestmark = [
  pytest.mark.e2e,
  pytest.mark.skipif(
    not _aligners_available(),
    reason='e2e: needs blastp/mmseqs2/ARC binaries (run `pytest -m e2e` on arborist-dev)',
  ),
]


@pytest.fixture(scope='session')
def species_path() -> Path:
  return TESTS_DIR / 'build' / 'species'


@pytest.fixture(scope='session', autouse=True)
def generated_run(species_path):
  # Generate assignments for the golden species in-place (native pytest -- no
  # Makefile wrapper), then let the tests diff against the committed goldens.
  env = {**os.environ}
  env['PYTHONPATH'] = os.pathsep.join(filter(None, [str(SRC_PROTEIN), env.get('PYTHONPATH', '')]))
  subprocess.run(
    [sys.executable, str(ASSIGN_PY), '-b', 'build', '-n', '4'],
    cwd=TESTS_DIR, env=env, check=True,
  )
  yield
  files_to_remove = [
    'peptide-assignments.tsv',
    'proteome.db',
    'source-data.tsv',
    'peptide-matches.tsv',
    'arc-results.tsv',
  ]
  for taxon in species:
    species_dir = species_path / str(taxon)
    for file in files_to_remove:
      file_path = species_dir / file
      if file_path.exists():
        file_path.unlink()


def _cols_to_dict(df, key_col, value_col):
  return dict(df.select(key_col, value_col).iter_rows())


def _check_column_matches(species_path, column_name):
  for taxon in species:
    generated_df = pl.read_csv(species_path / str(taxon) / 'peptide-assignments.tsv', separator='\t')
    expected_df = pl.read_csv(species_path / str(taxon) / 'expected-assignments.tsv', separator='\t')

    generated_dict = _cols_to_dict(generated_df, 'Epitope Sequence', column_name)
    expected_dict = _cols_to_dict(expected_df, 'Epitope Sequence', column_name)
    assert generated_dict == expected_dict, f"Mismatch in column '{column_name}' for taxon {taxon}"


def test_generated_files_exist(species_path):
  for taxon in species:
    assert (species_path / str(taxon) / 'peptide-assignments.tsv').exists(), 'Assignment file not generated'


def test_assignment_output(species_path):
  _check_column_matches(species_path, 'Assigned Protein ID')


def test_fragment_output(species_path):
  _check_column_matches(species_path, 'Assigned Protein Fragments')


def test_synonym_output(species_path):
  _check_column_matches(species_path, 'Assigned Protein Synonyms')
