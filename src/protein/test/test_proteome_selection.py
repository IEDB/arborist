import unittest
import os
import shutil
import pandas as pd

from pathlib import Path
from protein_tree.select_proteome import ProteomeSelector

build_path = Path(__file__).parent / 'build'

class TestProteomeSelection(unittest.TestCase):
  def test_proteome_list_pull(self):
    selector = ProteomeSelector(taxon_id=12345, species_name='Gaunavirus GA1', group='virus', build_path=build_path)
    self.assertIsInstance(selector.proteome_list, pd.DataFrame)

  def test_gp_proteome_pull(self):
    selector = ProteomeSelector(taxon_id=694003, species_name='Betacoronavirus 1', group='virus', build_path=build_path)
    selector._get_gp_proteome_to_fasta('UP000007552', '31631')
    
    file_path = build_path / 'species' / '694003' / 'gp_proteome.fasta'
    self.assertTrue(os.path.getsize(file_path) > 0)

  def test_rest_proteome_pull(self):
    species_path = build_path / 'species' / '408688'
    ProteomeSelector.get_proteome_to_fasta('UP000000273', species_path)

    file_path = species_path / 'UP000000273.fasta'
    self.assertTrue(os.path.getsize(file_path) > 0)

  def test_taxonomy_protein_map_pull(self):
    ProteomeSelector.get_all_proteins(taxon_id=161600, species_path=build_path / 'species' / '161600')
    file_path = build_path / 'species' / '161600' / 'proteome.fasta'
    self.assertTrue(os.path.getsize(file_path) > 0)

  def tearDown(self):
    species_path = build_path / 'species'
    if species_path.exists() and species_path.is_dir():
      shutil.rmtree(species_path)