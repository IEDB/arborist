import argparse
import re
import requests
import gzip
import json
import polars as pl

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

    self.proteome_list = self._get_candidate_proteomes()
    self.num_proteomes = len(self.proteome_list) + 1
  
  def select(self):
    if self.proteome_list.is_empty():
      self._get_orphans()
      self._preprocess_proteome_if_needed()
      self._write_metadata('', self.taxon_id, 'Orphans', self.species_name)
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
          proteome_type = proteome_type
          break

      if selected_proteomes.height > 1:
        proteome_id, proteome_taxon, proteome_label = self._proteome_tiebreak(selected_proteomes)
        self._remove_unselected_proteomes(proteome_id)
      else:
        proteome_id = selected_proteomes.item(0, 'Proteome ID')
        proteome_taxon = selected_proteomes.item(0, 'Proteome Taxon ID')
        proteome_label = selected_proteomes.item(0, 'Proteome Label')
        if not (self.species_path / '{proteome_id}.fasta').exists():
          self._fetch_proteome_file(proteome_id)

      self._preprocess_proteome_if_needed(proteome_id)
      self._rename_proteome_files(proteome_id)
      self._fetch_fragment_data(proteome_id)
      self._fetch_synonym_data(proteome_id)
      self._fetch_gp_proteome(proteome_id, proteome_taxon)
      self._write_metadata(proteome_id, proteome_taxon, proteome_type, proteome_label)

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
      proteome_taxon = proteome.item(0, 'Proteome Taxon ID')
      self._fetch_proteome_file(proteome_id)
      return proteome_id, proteome_taxon
    else:
      for proteome in selected_proteomes.rows(named=True):
        proteome_id = proteome['Proteome ID']
        proteome_taxon = proteome['Proteome Taxon ID']
        proteome_label = proteome['Proteome Label']
        
        if not (self.species_path / f'{proteome_id}.fasta').exists():
          self._fetch_proteome_file(proteome_id)

        match_count = self._get_match_count(peptide_seqs, proteome_id)
        proteome_counts[(proteome_id, proteome_taxon, proteome_label)] = match_count

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
  
  def _preprocess_proteome_if_needed(self, proteome_id: str = 'proteome'):
    if not (self.species_path / f'{proteome_id}.db').exists():
      Preprocessor(
        proteome = self.species_path / f'{proteome_id}.fasta',
        preprocessed_files_path = self.species_path,
      ).sql_proteome(k = 5)

  def _rename_proteome_files(self, proteome_id: str):
    (self.species_path / f'{proteome_id}.fasta').rename(self.species_path / 'proteome.fasta')
    (self.species_path / f'{proteome_id}.db').rename(self.species_path / 'proteome.db')

  def _fetch_fragment_data(self, proteome_id: str):
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=json&query=proteome:{proteome_id}&fields=ft_chain,ft_peptide,ft_propep,ft_signal,ft_transit'
    
    try:
      r = requests.get(url)
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      self.get_fragment_data(proteome_id)

    r.raise_for_status()
    data = json.loads(r.text)

    fragment_map = {}
    for entry in data['results']: # protein entry loop
      uniprot_id = entry['primaryAccession']
      fragments = []
      features = entry.get('features', [])
      for feature in features: # fragment loop
        fragment_data = {
          'type': feature['type'],
          'start': feature['location']['start']['value'],
          'end': feature['location']['end']['value'],
          'description': feature.get('description', 'N/A'),
          'feature_id': feature.get('featureId', 'N/A')
        }
        fragments.append(fragment_data)

      if len(fragments) > 1:
        fragments = self._filter_fragment_data(fragments)

      fragment_map[uniprot_id] = fragments
    
    with open(self.species_path / 'fragment-data.json', 'w') as f:
      json.dump(fragment_map, f, indent=2)

  def _filter_fragment_data(self, fragments: list):
    if all(f['type'] == 'Chain' for f in fragments[:2]) and \
      fragments[0]['start'] == 1 and fragments[1]['start'] == 2 and \
      fragments[0]['end'] == fragments[1]['end']: # remove overlapping chains that differ by 1st residue
        fragments = [fragments[0]] + fragments[2:] if len(fragments) > 2 else [fragments[0]]
    return fragments 

  def _fetch_synonym_data(self, proteome_id: str):
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=json&query=proteome:{proteome_id}&fields=protein_name'

    try:
      r = requests.get(url)
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      self.get_synonym_data(proteome_id)

    r.raise_for_status()
    data = json.loads(r.text)

    synonym_map = {}
    for entry in data['results']:
      uniprot_id = entry['primaryAccession']
      alternative_names = [
        name['fullName']['value'] 
        for name in entry.get('proteinDescription', {}).get('alternativeNames', [])
        if 'fullName' in name and 'value' in name['fullName']
      ]
      if alternative_names:
        synonym_map[uniprot_id] = alternative_names

    with open(self.species_path / 'synonym-data.json', 'w') as f:
      json.dump(synonym_map, f, indent=2)

  def _fetch_gp_proteome(self, proteome_id: str, proteome_taxon: int):
    ftp_url = f'https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/'
    if self.group == 'archeobacterium':
      ftp_url += f'Archaea/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif self.group == 'bacterium':
      ftp_url += f'Bacteria/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif self.group in ['plant', 'vertebrate', 'other-eukaryote']:
      ftp_url += f'Eukaryota/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif self.group == 'virus':
      ftp_url += f'Viruses/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    else:
      return

    try:
      with requests.get(ftp_url, stream=True) as r:
        r.raise_for_status()
        with open(f'{self.species_path}/gp_proteome.fasta', 'wb') as f:
          f.write(gzip.open(r.raw, 'rb').read())
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon) # recursive call on error
    except requests.exceptions.HTTPError:
      return

  def _write_metadata(self, proteome_id: str, proteome_taxon: int, proteome_type: str, proteome_label: str):
    pl.DataFrame({
      'Species Taxon ID': self.taxon_id, 'Species Name': self.species_name, 'Proteome ID': proteome_id,
      'Proteome Taxon ID': proteome_taxon, 'Proteome Type': proteome_type, 'Proteome Label': proteome_label
    }).write_csv(self.species_path / 'species-data.tsv', separator='\t')

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
        'Proteome Label': proteome.get('taxonomy', {}).get('scientificName'),
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
  
  def to_tsv(self):
    regexes = {
      'protein_id': re.compile(r"\|([^|]*)\|"),     # between | and |
      'entry_name': re.compile(r"\|([^\|]*?)\s"),   # between | and space
      'protein_name': re.compile(r"\s(.+?)\sOS"),   # between space and space before OS
      'gene': re.compile(r"GN=([^\s]+?)(?=\s|$)"),  # between GN= and space or end of line
      'pe_level': re.compile(r"PE=(.+?)\s"),        # between PE= and space
    }

    proteins = list(SeqIO.parse(self.species_path / 'proteome.fasta', 'fasta'))
    gp_proteome_path = self.species_path / 'gp_proteome.fasta'
    if gp_proteome_path.exists():
      gp_ids = [str(protein.id.split('|')[1]) for protein in list(SeqIO.parse(gp_proteome_path, 'fasta'))]
    else:
      gp_ids = []

    proteome_data = []
    for protein in proteins:
      metadata = []
      for key in regexes:
        match = regexes[key].search(str(protein.description))

        if match:
          metadata.append(match.group(1))
        else:
          if key == 'protein_id':
            metadata.append(str(protein.id))
          elif key == 'pe_level':
            metadata.append(0)
          else:
            metadata.append('')
      
      gp = 1 if protein.id.split('|')[1] in gp_ids else 0
      metadata.append(gp)
      metadata.append(str(protein.seq))
      metadata.append(protein.id.split('|')[0])
      
      proteome_data.append(metadata)

    columns = [
      'Protein ID', 'Entry Name', 'Protein Name', 'Gene', 'Protein Existence Level', 
      'Gene Priority', 'Sequence', 'Database'
    ]
    proteome = pl.DataFrame(proteome_data, schema=columns, orient='row').with_columns(
      pl.when(pl.col('Protein ID').str.contains('-'))
      .then(pl.col('Protein ID').str.split('-').list.last())
      .otherwise(pl.lit('1')).alias('Isoform Count')
    ).select([
      'Database', 'Gene', 'Protein ID', 'Entry Name', 'Isoform Count', 'Protein Name', 
      'Protein Existence Level', 'Gene Priority', 'Sequence'
    ])
    proteome.write_csv(self.species_path / 'proteome.tsv', separator='\t')

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
  proteome_selector.to_tsv()


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
