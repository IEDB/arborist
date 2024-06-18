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

def test_unique_source_parent_ids(files_path):
  df = pl.read_csv(files_path / 'source-parents.tsv', separator='\t')
  assert df['Source ID'].n_unique() == df.shape[0]

def test_unique_parent_protein_ids(files_path):
  df = pl.read_csv(files_path / 'parent-proteins.tsv', separator='\t')
  assert df['Accession'].n_unique() == df.shape[0]

def test_parent_protein_cols(files_path):
  cols = [
    'Accession', 'Database', 'Name', 'Title', 'Proteome ID', 'Proteome Label', 'Sequence'
  ]
  df = pl.read_csv(files_path / 'parent-proteins.tsv', separator='\t')
  assert df.columns == cols

def test_no_empty_epitope_mapping_ids(files_path):
  df = pl.read_csv(files_path / 'epitope-mappings.tsv', separator='\t')
  assert df['epitope_id'].null_count() == 0