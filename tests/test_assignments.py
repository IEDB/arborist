import pytest
import polars as pl
from pathlib import Path

species = [108, 9593, 9606, 3052236]

@pytest.fixture(scope='session')
def species_path() -> Path:
  return Path(__file__).parent / 'build' / 'species'

@pytest.fixture(scope='session', autouse=True)
def cleanup_files(species_path):
  yield
  files_to_remove = [
    'peptide-assignments.tsv',
    'proteome.db',
    'source-data.tsv',
    'peptide-matches.tsv',
    'arc-results.tsv'
  ]
  
  for taxon in species:
    species_dir = species_path / str(taxon)
    for file in files_to_remove:
      file_path = species_dir / file
      if file_path.exists():
        file_path.unlink()

def _cols_to_dict(df, key_col, value_col):
  return dict(df.select(key_col,value_col).iter_rows())

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
