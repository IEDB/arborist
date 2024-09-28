import argparse
import re
import os
import csv
import requests
import gzip
import json
import polars as pl
import xml.etree.ElementTree as ET

from typing import Iterator
from Bio import SeqIO
from pathlib import Path
from pepmatch import Preprocessor, Matcher

from protein_tree.data_fetch import DataFetcher


class ProteomeSelector:
  def __init__(self, taxon_id: int, species_name: str, group: str, peptides: pl.DataFrame):
    self.taxon_id = taxon_id
    self.species_name = species_name
    self.group = group
    self.peptides = peptides
    self.species_path = build_path / 'species' / str(taxon_id)
    self.species_path.mkdir(parents=True, exist_ok=True)

    self._get_candidate_proteomes()
  
  def select_best_proteome(self):
    pass

  def _get_candidate_proteomes(self):
    url = f'https://rest.uniprot.org/proteomes/stream?format=json&query=taxonomy_id:{self.taxon_id}'
    try:
      r = requests.get(url)
      r.raise_for_status()
      proteome_list = self._parse_proteome_json(r.text)
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      proteome_list = self._get_candidate_proteomes()
    
    proteome_list = proteome_list.filter(pl.col('Proteome Type') != 'Excluded')
    proteome_list.write_csv(self.species_path / 'proteome-list.tsv', separator='\t')

    return proteome_list

  def _parse_proteome_json(self, json_text: str):
    data = json.loads(json_text)
    
    rows = []
    for proteome in data['results']:
      row = {
        'Species Taxon ID': self.taxon_id,
        'Proteome ID': proteome.get('id'),
        'Proteome Type': proteome.get('proteomeType'),
        'Proteome Taxon ID': proteome.get('taxonomy', {}).get('taxonId'),
        'Gene Count': proteome.get('geneCount'),
        'Protein Count': proteome.get('proteinCount'),
        'BUSCO Score': proteome.get('proteomeCompletenessReport', {}).get('buscoReport', {}).get('score'),
        'Annotation Score': proteome.get('annotationScore')
      }
      rows.append(row)

    return pl.DataFrame(rows)

def get_proteome(taxon_id:int):
  species_row = active_species.row(by_predicate=pl.col('Species ID') == taxon_id)
  species_name = species_row[2]
  group = species_row[4]
  active_taxa = [int(taxon_id) for taxon_id in species_row[3].split(', ')]
  peptides = data_fetcher.get_peptides_for_species(all_peptides, active_taxa)
  config = {
    'taxon_id': taxon_id,
    'species_name': species_name,
    'group': group,
    'peptides': peptides,
  }
  proteome_selector = ProteomeSelector(**config)
  proteome_selector.select_best_proteome()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-t', '--taxon_id', type=int, help='Taxon ID of the species to process.',
  )
  parser.add_argument(
    '-b', '--build_path', type=str, help='Path for all Arborist build files.',
    default=Path(__file__).parents[3] / 'build'
  )
  args = parser.parse_args()

  taxon_id = args.taxon_id
  build_path = Path(args.build_path)
  all_species = not bool(taxon_id)

  active_species = pl.read_csv(build_path / 'arborist' / 'active-species.tsv', separator='\t')
  
  data_fetcher = DataFetcher(build_path)
  all_peptides = data_fetcher.get_all_peptides()
  
  if all_species:
    for row in active_species.rows(named=True):
      get_proteome(row['Species ID'])
  else:
    get_proteome(taxon_id)
