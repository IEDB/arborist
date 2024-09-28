import argparse
import re
import requests
import gzip
import json
import polars as pl

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

    self.proteome_list = self._get_candidate_proteomes()
    self.num_proteomes = len(self.proteome_list) + 1
  
  def select(self):
    if self.proteome_list.is_empty():
      self._get_orphans()
    else:
      proteome_types = [
        'Reference and representative proteome',
        'Representative',
        'Reference proteome',
        'Other',
        'Redundant'
      ]

      for proteome_type in proteome_types:
        proteomes = self.proteome_list.filter(
          pl.col('Proteome Type').str.contains(proteome_type)
        )
        if not proteomes.is_empty():
          selected_proteomes = proteomes
          break

      if selected_proteomes.height > 1:
        proteome_id = self._proteome_tiebreak(selected_proteomes)
        self._remove_unselected_proteomes(proteome_id)
      else:
        proteome_id = selected_proteomes.item(0, 'Proteome ID')
        if not (self.species_path / f'{proteome_id}.fasta').exists():
          self._fetch_proteome_file(proteome_id)

      self._preprocess_proteome_if_needed(proteome_id)
      self._rename_proteome_files(proteome_id)
      # self._fetch_fragment_data(proteome_id)
      # self._fetch_synonym_data(proteome_id)
      

  def _get_orphans(self):
    url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&query=taxonomy_id:{taxon_id}&size=500'
    for batch in self._get_protein_batches(url):
      with open(self.species_path / 'proteome.fasta', 'a') as f:
        f.write(batch.text)

  def _get_protein_batches(self, batch_url: str):
    while batch_url:
      try:
        r = requests.get(batch_url)
        r.raise_for_status()
        yield r
        batch_url = self._get_next_link(r.headers)
      except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
        yield from self._get_protein_batches(batch_url)

  def _get_next_link(self, headers: dict):
    re_next_link = re.compile(r'<(.+)>; rel="next"') # regex to extract URL
    if 'Link' in headers:
      match = re_next_link.match(headers['Link'])
      if match:
        return match.group(1)

  def _proteome_tiebreak(self, selected_proteomes: pl.DataFrame):
    proteome_counts = {}
    peptide_seqs = self.peptides['Sequence'].to_list()
    if selected_proteomes.height > 20:
      proteome = selected_proteomes.sort('BUSCO Score').tail(1)
      proteome_id = proteome.item(0, 'Proteome ID')
      self._fetch_proteome_file(proteome_id)
      return proteome_id
    else:
      for proteome in selected_proteomes.rows(named=True):
        proteome_id = proteome['Proteome ID']
        
        if not (self.species_path / f'{proteome_id}.fasta').exists():
          self._fetch_proteome_file(proteome_id)

        match_count = self._get_match_count(peptide_seqs, proteome_id)
        proteome_counts[proteome_id] = match_count

    return max(proteome_counts.items(), key=lambda item: item[1])[0]

  def _get_match_count(self, peptide_seqs: list, proteome_id: str):
    Preprocessor(
      proteome = self.species_path / f'{proteome_id}.fasta',
      preprocessed_files_path = self.species_path,
    ).sql_proteome(k = 5)
    
    matches = Matcher(
      query = peptide_seqs,
      proteome_file = self.species_path / f'{proteome_id}.fasta', 
      max_mismatches = 0, 
      k = 5,
      preprocessed_files_path = self.species_path,
      output_format='dataframe'
    ).match()
    matches = pl.from_pandas(matches)
    matches = matches.filter(
      pl.col('Matched Sequence').is_not_null()
    ).unique('Query Sequence')

    return matches.height

  def _fetch_proteome_file(self, proteome_id: str):
    url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
    try:
      with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(self.species_path / f'{proteome_id}.fasta', 'w') as f:
          for chunk in r.iter_content(chunk_size=8192):
            if chunk:  # filter out keep-alive new chunks
              f.write(chunk.decode())
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      ProteomeSelector.get_proteome_to_fasta(proteome_id, species_path)  # recursive call on error

  def _remove_unselected_proteomes(self, proteome_id: str):
    for file in self.species_path.glob('*.fasta'):
      if file.name != f'{proteome_id}.fasta':
        file.unlink()
    for file in self.species_path.glob('*.db'):
      if file.name != f'{proteome_id}.db':
        file.unlink()
  
  def _preprocess_proteome_if_needed(self, proteome_id: str):
    if not (self.species_path / f'{proteome_id}.db').exists():
      Preprocessor(
        proteome = self.species_path / f'{proteome_id}.fasta',
        preprocessed_files_path = self.species_path,
      ).sql_proteome(k = 5)

  def _rename_proteome_files(self, proteome_id: str):
    (self.species_path / f'{proteome_id}.fasta').rename(self.species_path / 'proteome.fasta')
    (self.species_path / f'{proteome_id}.db').rename(self.species_path / 'proteome.db')

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

    proteome_list = pl.DataFrame(rows)
    if proteome_list.is_empty():
      proteome_list = pl.DataFrame([{
        'Species Taxon ID': None, 'Proteome ID': None, 'Proteome Type': None, 'Proteome Taxon ID': None, 
        'Gene Count': None, 'Protein Count': None, 'BUSCO Score': None, 'Annotation Score': None
      }]).head(0)

    return proteome_list

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
  proteome_selector.select()

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
