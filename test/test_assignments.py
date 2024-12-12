import pytest
import polars as pl
from pathlib import Path

species = [108, 3052236]

@pytest.fixture
def species_path() -> Path:
  return Path(__file__).parent / 'build' / 'species'

@pytest.fixture(autouse=True)
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

def test_assignment_output(species_path):
  for taxon in species:
    species_dir = species_path / str(taxon)
    
    generated_path = species_dir / 'peptide-assignments.tsv'
    assert generated_path.exists(), f"Generated assignments file not found for species {taxon}"
    generated_df = pl.read_csv(generated_path, separator='\t')
    
    expected_path = species_dir / 'expected-assignments.tsv'
    assert expected_path.exists(), f"Expected assignments file not found for species {taxon}"
    expected_df = pl.read_csv(expected_path, separator='\t')
    
    generated_df = generated_df.sort(generated_df.columns)
    expected_df = expected_df.sort(expected_df.columns)
    
    assert generated_df.columns == expected_df.columns, \
        f"Column mismatch for species {taxon}"
    
    assert generated_df.equals(expected_df)
