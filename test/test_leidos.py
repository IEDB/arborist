import pytest
import polars as pl
from pathlib import Path


@pytest.fixture
def files_path() -> Path:
  return Path(__file__).parents[1] / 'build' / 'proteins' / 'latest'

def test_source_parent_cols(files_path):
  cols = [
    'Source ID', 'Accession', 'Database', 'Name', 'Aliases', 'Synonyms', 'Taxon ID',
    'Taxon Name', 'Species ID', 'Species Label', 'Proteome ID', 'Proteome Label', 
    'Protein Strategy', 'Parent IRI', 'Parent Protein Database', 
    'Parent Protein Accession', 'Parent Sequence Length', 'Sequence', 'Parent Protein Gene'
  ]
  df = pl.read_csv(files_path / 'source-parents.tsv', separator='\t')
  assert df.columns == cols

def test_parent_protein_cols(files_path):
  cols = [
    'Accession', 'Database', 'Name', 'Title', 'Proteome ID', 'Proteome Label', 'Sequence'
  ]
  df = pl.read_csv(files_path / 'parent-proteins.tsv', separator='\t')
  assert df.columns == cols