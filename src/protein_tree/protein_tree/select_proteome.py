#!/usr/bin/env python3

import re
import os
import csv
import requests
import gzip
import json
import pandas as pd
import xml.etree.ElementTree as ET

from typing import Iterator
from Bio import SeqIO
from pathlib import Path
from pepmatch import Preprocessor, Matcher

from protein_tree.get_data import DataFetcher

class ProteomeSelector:
  def __init__(self, taxon_id, species_name, group, build_path = Path(__file__).parent.parent / 'build'):
    self.taxon_id = taxon_id
    self.group = group
    self.species_name = species_name
    self.species_path = build_path / 'species' / f'{taxon_id}'
    self.species_path.mkdir(parents=True, exist_ok=True)

    self.species_df = pd.read_csv(
      build_path / 'arborist' / 'active-species.tsv', sep='\t'
    )

    self.proteome_list = self._get_proteome_list() # get all candidate proteomes
    self.num_of_proteomes = len(self.proteome_list) + 1 # +1 for all proteins option

  def select_best_proteome(self, peptides_df: pd.DataFrame) -> list:
    """Select the best proteome to use for a species. Return the proteome ID,
    proteome taxon, and proteome type.
    
    Check UniProt for all candidate proteomes associated with that
    taxon. Then do the following checks:

    1. Are there any representative proteomes?
    2. Are there any reference proteomes?
    3. Are there any non-redudant proteomes?
    4. Are there any other proteomes?

    If yes to any of the above, check if there are multiples and do peptide
    search for tie breaks.

    If no to all of the above, then get every protein associated with
    the taxon ID using the get_all_proteins method.

    Args:
      peptides_df (pd.DataFrame): DataFrame of peptides for the species to use for tie breaks.
    """
    if self.proteome_list.empty:
      print('No proteomes found. Fetching orphan proteins.')
      self.get_all_proteins(self.taxon_id, self.species_path)
      proteome_id, proteome_taxon, proteome_type = '', self.taxon_id, 'Orphans'
      self._write_species_data(proteome_id, proteome_taxon, proteome_type)
      return [proteome_id, proteome_taxon, proteome_type]

    elif self.proteome_list['isRepresentativeProteome'].any():
      print('Found representative proteome(s).\n')
      proteome_type = 'Representative'
      self.proteome_list = self.proteome_list[self.proteome_list['isRepresentativeProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(peptides_df, True)
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)
    
    elif self.proteome_list['isReferenceProteome'].any():
      print('Found reference proteome(s).\n')
      proteome_type = 'Reference'
      self.proteome_list = self.proteome_list[self.proteome_list['isReferenceProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(peptides_df, True)
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)

    elif 'redundantTo' not in self.proteome_list.columns:
      print('Found other proteome(s).\n')
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(peptides_df, False)
    
    elif self.proteome_list['redundantTo'].isna().any():
      print('Found non-redundant proteome(s).\n')
      proteome_type = 'Non-redundant'
      self.proteome_list = self.proteome_list[self.proteome_list['redundantTo'].isna()]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(peptides_df, False)
    
    else:
      print('Found other proteome(s).\n')
      proteome_type = 'Other' # replace ID with redundant proteome ID
      self.proteome_list.loc[self.proteome_list['redundantTo'].notna(), 'upid'] = self.proteome_list['redundantTo']
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(peptides_df, False)

    self._remove_other_proteomes(proteome_id)
    self._write_species_data(proteome_id, proteome_taxon, proteome_type)

    return [proteome_id, proteome_taxon, proteome_type]

  def _write_species_data(self, proteome_id: str, proteome_taxon: str, proteome_type: str) -> None:
    """Write the proteome ID/taxon/type for a species to a CSV file for later use."""
    pd.DataFrame( # write proteome data to metrics file
      [[proteome_id, proteome_taxon, proteome_type, self.species_name]],
      columns=['Proteome ID', 'Proteome Taxon', 'Proteome Type', 'Species Name']
    ).to_csv(self.species_path / 'species-data.tsv', sep='\t', index=False)

  def proteome_to_tsv(self) -> None:
    """Write the proteome data for a species to a CSV file for later use."""
    
    if not (self.species_path / 'proteome.fasta').exists():
      return

    regexes = {
      'protein_id': re.compile(r"\|([^|]*)\|"),    # between | and |
      'entry_name': re.compile(r"\|([^\|]*?)\s"),   # between | and space
      'protein_name': re.compile(r"\s(.+?)\sOS"),  # between space and space before OS
      'gene': re.compile(r"GN=([^\s]+?)(?=\s|$)"), # between GN= and space or end of line
      'pe_level': re.compile(r"PE=(.+?)\s"),       # between PE= and space
    }

    proteins = list(SeqIO.parse(f'{self.species_path}/proteome.fasta', 'fasta'))
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
    
    columns = ['Protein ID', 'Entry Name', 'Protein Name', 'Gene', 'Protein Existence Level', 'Gene Priority', 'Sequence', 'Database']
    proteome = pd.DataFrame(proteome_data, columns=columns)
    
    proteome['Isoform Count'] = proteome['Protein ID'].apply(lambda x: x.split('-')[1] if '-' in x else '1')
    proteome = proteome[['Database', 'Gene', 'Protein ID', 'Entry Name', 'Isoform Count', 'Protein Name', 'Protein Existence Level', 'Gene Priority', 'Sequence']]
    proteome.to_csv(f'{self.species_path}/proteome.tsv', sep='\t', index=False)
  
  def _parse_proteome_xml(self, xml) -> pd.DataFrame:
    def str_to_bool(value):
      return value.lower() == 'true'

    root = ET.fromstring(xml)
    proteomes = root.findall('{http://uniprot.org/proteome}proteome')
    proteome_data = []
    for proteome in proteomes:
      excluded = True if proteome.find('{http://uniprot.org/proteome}excluded') is not None else False
      data = {
        'proteinCount': proteome.attrib['proteinCount'],
        'upid': proteome.find('{http://uniprot.org/proteome}upid').text,
        'taxonomy': proteome.find('{http://uniprot.org/proteome}taxonomy').text,
        'modified': proteome.find('{http://uniprot.org/proteome}modified').text,
        'isReferenceProteome': str_to_bool(proteome.find('{http://uniprot.org/proteome}isReferenceProteome').text),
        'isRepresentativeProteome': str_to_bool(proteome.find('{http://uniprot.org/proteome}isRepresentativeProteome').text),
        'excluded': excluded}
      for score_cols in proteome.findall('{http://uniprot.org/proteome}scores'):
        score_type = score_cols.attrib['name']
        for property in score_cols.findall('{http://uniprot.org/proteome}property'):
          prop_name = property.attrib['name']
          prop_value = property.attrib['value']
          data[f'{score_type}_{prop_name}'] = prop_value
      proteome_data.append(data)
    return pd.DataFrame(proteome_data)

  def _get_proteome_list(self) -> pd.DataFrame:
    """Get a list of proteomes for a species from the UniProt API.
    Check for proteome_type:1 first, which are the representative or
    reference proteomes.

    If there are no proteomes, return empty DataFrame.
    """
    # URL to get proteome list for a species - use proteome_type:1 first
    proteome_list_url = f'https://rest.uniprot.org/proteomes/stream?format=xml&query=taxonomy_id:{self.taxon_id}'
    try:
      r = requests.get(proteome_list_url)
      r.raise_for_status()
      proteome_list = self._parse_proteome_xml(r.text)
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      proteome_list = self._get_proteome_list()

    if 'excluded' in proteome_list.columns:
      proteome_list = proteome_list[proteome_list['excluded'] == False] # drop excluded
    proteome_list.to_csv(f'{self.species_path}/proteome-list.tsv', sep='\t', index=False)
    
    return proteome_list

  def _get_gp_proteome_to_fasta(self, proteome_id: str, proteome_taxon: str) -> None:
    """Write the gene priority proteome to a file. 
    This is only for representative and reference proteomes.
    Depending on the species group, the FTP URL will be different.

    Args:
      proteome_id (str): Proteome ID.
      proteome_taxon (str): Taxon ID for the proteome.
    """
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

  @staticmethod
  def get_proteome_to_fasta(proteome_id: str, species_path: Path) -> None:
    """Get the FASTA file for a proteome from UniProt API.
    Include all isoforms and do not compress the file.

    Args:
      proteome_id (str): UniProt Proteome ID.
      species_path (Path): Path to species directory to write the FASTA file.
    """
    url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
    species_path.mkdir(parents=True, exist_ok=True)
    try:
      with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(f'{species_path}/{proteome_id}.fasta', 'w') as f:
          for chunk in r.iter_content(chunk_size=8192):
            if chunk:  # filter out keep-alive new chunks
              f.write(chunk.decode())
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      ProteomeSelector.get_proteome_to_fasta(proteome_id, species_path)  # recursive call on error

  @staticmethod
  def get_all_proteins(taxon_id: int, species_path: Path) -> None:
    """Get every protein associated with a taxon ID on UniProt.
    Species on UniProt will have a proteome, but not every protein is
    stored within those proteomes. There is a way to get every protein
    using the taxonomy part of UniProt.

    Args:
      taxon_id (int): Taxon ID for the species.
      species_path (Path): Path to species directory.
    """
    def get_protein_batches(batch_url: str) -> Iterator[requests.Response]:
      while batch_url:
        try:
          r = requests.get(batch_url)
          r.raise_for_status()
          yield r
          batch_url = get_next_link(r.headers)
        except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
          yield from get_protein_batches(batch_url)

    def get_next_link(headers: dict) -> str:
      re_next_link = re.compile(r'<(.+)>; rel="next"') # regex to extract URL
      if 'Link' in headers:
        match = re_next_link.match(headers['Link'])
        if match:
          return match.group(1)

    # URL link to all proteins for a species - size = 500 proteins at a time
    url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&query=taxonomy_id:{taxon_id}&size=500'
    species_path.mkdir(parents=True, exist_ok=True)
    for batch in get_protein_batches(url):
      with open(f'{species_path}/proteome.fasta', 'a') as f:
        f.write(batch.text)

  def _get_proteome_with_most_matches(self, peptides_df: pd.DataFrame, rep_or_ref: bool) -> tuple:
    """Get the proteome ID and true taxon ID associated with
    that proteome with the most peptide matches in case there is a tie.
    We get the true taxon so we can extract the data from the FTP server
    if needed.

    Args:
      peptides_df (pd.DataFrame): DataFrame of peptides for the species.
    """
    if self.num_of_proteomes <= 2: # excludes all proteins option
      proteome_id = self.proteome_list['upid'].iloc[0]
      proteome_taxon = self.proteome_list['taxonomy'].iloc[0]
      ProteomeSelector.get_proteome_to_fasta(proteome_id, self.species_path)
      return proteome_id, proteome_taxon
    
    # use max BUSCO score if there are too many proteomes
    elif not rep_or_ref and self.num_of_proteomes >= 20 and 'busco_score' in self.proteome_list.columns:
      idx = self.proteome_list['busco_score'].idxmax()
      proteome_id = self.proteome_list['upid'].loc[idx]
      proteome_taxon = self.proteome_list['taxonomy'].loc[idx]
      ProteomeSelector.get_proteome_to_fasta(proteome_id, self.species_path)
      return proteome_id, proteome_taxon

    peptides_df = peptides_df[peptides_df['Sequence'].notna()] 
    peptides = peptides_df['Sequence'].tolist()

    if not peptides: # get proteome with highest protein count
      idx = self.proteome_list['proteinCount'].idxmax()
      proteome_id = self.proteome_list['upid'].loc[idx]
      ProteomeSelector.get_proteome_to_fasta(proteome_id, self.species_path)
      return proteome_id, self.taxon_id

    match_counts = {} # keep track of # of peptide matches for each proteome
    for proteome_id in list(self.proteome_list['upid']):
      ProteomeSelector.get_proteome_to_fasta(proteome_id, self.species_path)
      
      Preprocessor(
        proteome = f'{self.species_path}/{proteome_id}.fasta',
        preprocessed_files_path = f'{self.species_path}',
      ).sql_proteome(k = 5)

      matches_df = Matcher(
        query = peptides, 
        proteome_file = f'{self.species_path}/{proteome_id}.fasta', 
        max_mismatches = 0, 
        k = 5,
        preprocessed_files_path = f'{self.species_path}',
        output_format='dataframe'
      ).match()
      
      matches_df.drop_duplicates(subset=['Query Sequence'], inplace=True)
      
      try:
        match_counts[proteome_id] = matches_df['Matched Sequence'].dropna().count()
      except KeyError: # in case there are no matches
        match_counts[proteome_id] = 0 

    # select the proteome ID and proteome taxon with the most matches
    proteome_id = max(match_counts, key=match_counts.get)
    proteome_taxon = self.proteome_list[self.proteome_list['upid'] == proteome_id]['taxonomy'].iloc[0]

    return proteome_id, proteome_taxon
  
  def get_fragment_data(self, proteome_id) -> None:
    """Get the fragment data for a proteome from UniProt API which includes:
    chains, initiator methionine, peptides, propeptides, signal peptides, and transit peptides."""

    url = f'https://rest.uniprot.org/uniprotkb/stream?format=json&query=proteome:{proteome_id}&fields=ft_chain,ft_peptide,ft_propep,ft_signal,ft_transit'
    
    try:
      r = requests.get(url)
      r.raise_for_status()
      data = r.json()
      
      fragments = {}
      for entry in data["results"]:

        try: # some entries do not have features
          uniprot_id = entry["primaryAccession"]
          fragment_list = []
          for feature in entry["features"]:
            feature_type = feature["type"]
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            fragment_list.append(f"{feature_type}-{start}-{end}")

          if fragment_list:
            fragments[uniprot_id] = ", ".join(fragment_list)
        except KeyError:
          continue
      
      with open(f'{self.species_path}/fragment-data.json', 'w') as f:
        json.dump(fragments, f, indent=4)

    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      self.get_fragment_data(proteome_id)
  
  def get_synonym_data(self, proteome_id) -> None:
    """Get the synonym names of each protein for a proteome from UniProt API."""

    url = f'https://rest.uniprot.org/uniprotkb/stream?format=json&query=proteome:{proteome_id}&fields=protein_name'
    try:
      r = requests.get(url)
      r.raise_for_status()
      data = r.json()

      uniprot_alternative_names = {}
      for result in data["results"]:
        uniprot_id = result["primaryAccession"]
        alternative_names = [name["fullName"]["value"] for name in result["proteinDescription"].get("alternativeNames", []) if "fullName" in name and "value" in name["fullName"]]
        if alternative_names:
          uniprot_alternative_names[uniprot_id] = alternative_names

      if len(uniprot_alternative_names) != 0:
        with open(f'{self.species_path}/synonym-data.json', 'w') as f:
          json.dump(uniprot_alternative_names, f, indent=4)

    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout):
      self.get_synonym_data(proteome_id)

  def _remove_other_proteomes(self, proteome_id: str) -> None:
    """Remove the proteome FASTA files that are not the chosen proteome for that
    species. Also, remove the .db files and rename the chosen proteome to 
    "proteome.fasta".

    Args:
      proteome_id (str): Proteome ID of the chosen proteome.
    """
    proteome_list_to_remove = self.proteome_list[self.proteome_list['upid'] != proteome_id]
    for i in list(proteome_list_to_remove['upid']):
      (self.species_path / f'{i}.fasta').unlink(missing_ok=True)
      (self.species_path / f'{i}.db').unlink(missing_ok=True)

    Path(f'{self.species_path}/{proteome_id}.db').unlink(missing_ok=True)   
    os.rename(f'{self.species_path}/{proteome_id}.fasta', f'{self.species_path}/proteome.fasta')

def update_proteome(species_path: Path, taxon_id: int, data_path: Path) -> None:
  """Update the proteome for a species if there are new peptides.
  
  Args:
    species_path (Path): Path to species directory.
    taxon_id (int): Taxon ID for the species.
    data_path (Path): Path to species-data.tsv.
  """
  with open(data_path, newline='') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
      proteome_id = row['Proteome ID']
      proteome_taxon = row['Proteome Taxon']
      proteome_type = row['Proteome Type']
  
  if proteome_type == 'Orphans':
    ProteomeSelector.get_all_proteins(taxon_id, species_path)
  else:
    ProteomeSelector.get_proteome_to_fasta(proteome_id, species_path)
    os.rename(f'{species_path}/{proteome_id}.fasta', f'{species_path}/proteome.fasta')
  
  return [proteome_id, proteome_taxon, proteome_type]

def run(taxon_id: int, species_name: str, group: str, all_taxa: list, build_path: Path, all_peptides: pd.DataFrame, force: bool) -> list:
  """Run the proteome selection process for a species.
  
  Args:
    taxon_id (int): Taxon ID for the species.
    species_name (str): Name of the species.
    group (str): Group for the species (e.g. bacterium, virus, etc.).
    all_taxa (list): List of all children taxa for the species.
    build_path (Path): Path to build directory.
    all_peptides (pd.DataFrame): DataFrame of peptides for the species.
    force (bool): Force reselection of proteome.
  """
  species_path = build_path / 'species' / f'{taxon_id}'
  species_path.mkdir(parents=True, exist_ok=True)

  if not force and (data_path := species_path / 'species-data.tsv').exists():
    print(f'Updating proteome for species w/ taxon ID: {taxon_id}).')
    proteome_data = update_proteome(species_path, taxon_id, data_path)
    return proteome_data

  Fetcher = DataFetcher(build_path)
  peptides_df = Fetcher.get_peptides_for_species(all_peptides, all_taxa)

  print(f'Selecting best proteome for species w/ taxon ID: {taxon_id}).')
  
  Selector = ProteomeSelector(taxon_id, species_name, group, build_path)
  print(f'Number of candidate proteomes: {Selector.num_of_proteomes}')

  proteome_data = Selector.select_best_proteome(peptides_df)
  Selector.proteome_to_tsv()

  if proteome_data[2] != 'Orphans': # get fragment and synonym data for proteome
    Selector.get_fragment_data(proteome_data[0])
    Selector.get_synonym_data(proteome_data[0])

  # sanity check to make sure proteome.fasta is not empty
  if (species_path / 'proteome.fasta').stat().st_size == 0:
    proteome_data[0] = 'None'
    proteome_data[1] = taxon_id
    proteome_data[2] = 'Orphans'
    ProteomeSelector.get_all_proteins(taxon_id, species_path)
  
  print(f'Proteome ID: {proteome_data[0]}')
  print(f'Proteome taxon: {proteome_data[1]}')
  print(f'Proteome type: {proteome_data[2]}\n')

  return proteome_data

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser()
  
  parser.add_argument(
    '-b', '--build_path', 
    type=str,
    default=Path(__file__).parent.parent / 'build',
    help='Path to data directory.'
  )
  parser.add_argument(
    '-a', '--all_species', 
    action='store_true', 
    help='Build protein tree for all IEDB species.'
  )
  parser.add_argument(
    '-t', '--taxon_id', 
    type=int,
    help='Taxon ID for the species to pull data for.'
  )
  parser.add_argument(
    '-f', '--force',
    action='store_true',
    help='Force reselection of proteome.'
  )
  
  args = parser.parse_args()

  build_path = Path(args.build_path)
  all_species = args.all_species
  taxon_id = args.taxon_id
  force = args.force

  assert all_species or taxon_id, 'Please specify either --all_species or --taxon_id.'
  assert (all_species and not taxon_id) or (taxon_id and not all_species), 'Please specify either --all_species or --taxon_id.'

  # TODO: replace the data/active-species.tsv with updated arborist active-species somehow
  species_df = pd.read_csv(build_path / 'arborist' / 'active-species.tsv', sep='\t')
  valid_taxon_ids = species_df['Species ID'].tolist()

  all_peptides = DataFetcher(build_path).get_all_peptides()

  all_taxa_map = dict(zip( # map taxon ID to list of all children taxa
    species_df['Species ID'],
    species_df['Active Taxa']))
  taxon_to_species_map = dict(zip( # map taxon ID to species name
    species_df['Species ID'],
    species_df['Species Label']))
  
  if all_species: # run all species at once
    for taxon_id in valid_taxon_ids:
      species_name = taxon_to_species_map[taxon_id]
      group = species_df[species_df['Species ID'] == taxon_id]['Group'].iloc[0]
      all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]
      proteome_data = run(taxon_id, species_name, group, all_taxa, build_path, all_peptides, force)

  else: # one species at a time
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    species_name = taxon_to_species_map[taxon_id]
    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]
    group = species_df[species_df['Species ID'] == taxon_id]['Group'].iloc[0]
    proteome_data = run(taxon_id, species_name, group, all_taxa, build_path, all_peptides, force)