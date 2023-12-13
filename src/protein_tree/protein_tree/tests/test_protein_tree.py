import pytest
import pandas as pd
from pathlib import Path

import os
import sys
import glob

# add path to parent directory to sys.path so that we can import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from select_proteome import ProteomeSelector
from assign_gene_protein import GeneAndProteinAssigner


data_path = Path(__file__).parent / "data"


@pytest.fixture(
  params=[
    [106, 'Runella slithyformis', False, 'UP000000493'],
    [10042, 'Peromyscus maniculatus', True, 'UP000504601'],
    [334205, 'Nupapillomavirus 1', False, 'UP000006367']
  ]
)
def organism(request):
  return request.param


@pytest.fixture
def species_path(organism) -> Path:
  taxon_id, species_name, _, _ = organism
  return data_path / "species" / f"{taxon_id}-{species_name.replace(' ', '-')}"


@pytest.fixture
def epitopes(species_path) -> Path:
  return species_path / "epitopes.tsv"


@pytest.fixture
def sources(species_path) -> Path:
  return species_path / "sources.tsv"


@pytest.fixture(scope='session', autouse=True)
def cleanup_after_all_tests():
  yield

  file_names = ['*proteome*', '*ARC_results*']
  files_to_remove = []
  for name in file_names:
    path = data_path / 'species' / '**' / name
    files_to_remove.extend(glob.glob(str(path), recursive=True))

  for file in files_to_remove:
    os.remove(file)


def test_select_proteome(species_path, epitopes, organism):
  taxon_id, _, _, proteome_id = organism 

  epitopes_df = pd.read_csv(epitopes, sep='\t')
  Selector = ProteomeSelector(taxon_id, species_path)
  proteome_data = Selector.select_best_proteome(epitopes_df)
  Selector.proteome_to_tsv()

  assert proteome_data[0] == proteome_id


def test_assignments(species_path, epitopes, sources, organism):
  taxon_id, _, is_vertebrate, _ = organism

  epitopes_df = pd.read_csv(epitopes, sep='\t')
  sources_df = pd.read_csv(sources, sep='\t')

  Assigner = GeneAndProteinAssigner(
    taxon_id, species_path, is_vertebrate, num_threads=1, data_path=data_path, bin_path='/usr/bin'
  )
  _, epitope_assignments, source_assignments = Assigner.assign(sources_df, epitopes_df)
  
  epitopes_expected = pd.read_csv(species_path / 'epitope_assignments.tsv', sep='\t')
  sources_expected = pd.read_csv(species_path / 'source_assignments.tsv', sep='\t')

  assert epitopes_expected.equals(epitope_assignments)
  assert sources_expected.equals(source_assignments)
