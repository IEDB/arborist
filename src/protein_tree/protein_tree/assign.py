#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import os
import glob
import pandas as pd

from ARC.classifier import SeqClassifier
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from pepmatch import Preprocessor, Matcher

from protein_tree.get_data import DataFetcher
from protein_tree.select_proteome import ProteomeSelector

class GeneAndProteinAssigner:
  def __init__(
    self,
    taxon_id,
    species_path,
    is_vertebrate,
    group,
    peptides_df,
    num_threads,
    build_path = Path(__file__).parent.parent / 'build',
  ):
    
    self.species_path = species_path
    self.taxon_id = taxon_id
    self.is_vertebrate = is_vertebrate
    self.build_path = build_path
    self.num_threads = num_threads

    # initialize dicts for assignments
    self.source_gene_assignment = {}
    self.source_protein_assignment = {}
    self.source_assignment_score = {}
    self.source_arc_assignment = {}
    self.peptide_protein_assignment = {}

    try: # read in proteome data
      self.proteome = pd.read_csv(f'{self.species_path}/proteome.tsv', sep='\t')
    except FileNotFoundError: # if proteome data doesn't exist, select best proteome
      print('No proteome file found. Selecting best proteome...')
      selector = ProteomeSelector(taxon_id, species_name, group, build_path)
      selector.select_best_proteome(peptides_df)
      selector.proteome_to_tsv()
      self.proteome = pd.read_csv(f'{self.species_path}/proteome.tsv', sep='\t')

    self.uniprot_id_to_gene_map = dict( # create UniProt ID -> gene map
      zip(
        self.proteome['Protein ID'], 
        self.proteome['Gene']
      )
    )
    # create UniProt ID -> protein name map
    self.uniprot_id_to_name_map = dict(
      zip(
        self.proteome['Protein ID'], 
        self.proteome['Protein Name']
      )
    )

  def assign(self, sources_df: pd.DataFrame, peptides_df: pd.DataFrame) -> None:
    """Overall function to assign genes and parent proteins to sources and peptides.

    Args:
      sources_df: DataFrame of peptide sources for a species.
      peptides_df: DataFrame of peptides for a species.
    """
    # assign None to all peptide sources and peptides to start
    for i, row in sources_df.iterrows():
      self.source_gene_assignment[row['Accession']] = None
      self.source_protein_assignment[row['Accession']] = None
      self.source_assignment_score[row['Accession']] = None
    for i, row in peptides_df.iterrows():
      self.peptide_protein_assignment[(row['Source Accession'], row['Sequence'])] = None

    self.source_to_peptides_map = self._create_source_to_peptides_map(peptides_df)

    sources_df['Length'] = sources_df['Sequence'].str.len() # add length column to sources
    self.source_length_map = dict( # create map of peptide sources to their length
      zip(
        sources_df['Accession'],
        sources_df['Length']
      )
    )
    print('Assigning peptide sources...')
    self._assign_sources(sources_df)
    print('Done assigning peptide sources.\n')

    self._assign_allergens()
    self._assign_manuals()

    print('Assigning peptides...')
    self._assign_peptides(peptides_df)
    print('Done.\n')

    # map peptide sources to their BLAST matches and ARC assignments
    sources_df.loc[:, 'Assigned Gene'] = sources_df['Accession'].map(self.source_gene_assignment)
    sources_df.loc[:, 'Assigned Gene'] = sources_df['Assigned Gene'].fillna(sources_df['Accession'].map(self.source_arc_assignment))
    sources_df.loc[:, 'Assigned Protein ID'] = sources_df['Accession'].map(self.source_protein_assignment)
    sources_df.loc[:, 'Assigned Protein Name'] = sources_df['Assigned Protein ID'].map(self.uniprot_id_to_name_map)
    sources_df.loc[:, 'Assignment Score'] = sources_df['Accession'].map(self.source_assignment_score)
    sources_df.loc[:, 'ARC Assignment'] = sources_df['Accession'].map(self.source_arc_assignment)

    # map peptide peptide sources to assignments above and then PEPMatch assignments
    peptides_df.loc[:, 'Parent Antigen ID'] = peptides_df['Source Accession'].map(self.source_protein_assignment)
    peptides_df.loc[:, 'Parent Antigen Name'] = peptides_df['Parent Antigen ID'].map(self.uniprot_id_to_name_map)
    peptides_df.loc[:, 'Parent Antigen Gene'] = peptides_df['Source Accession'].map(self.source_gene_assignment)
    peptides_df.loc[:, 'Parent Antigen Gene'] = peptides_df['Parent Antigen Gene'].fillna(peptides_df['Source Accession'].map(self.source_arc_assignment))
    peptides_df.set_index(['Source Accession', 'Sequence'], inplace=True)
    peptides_df.loc[:, 'Parent Antigen Gene Isoform ID'] = peptides_df.index.map(self.peptide_protein_assignment)
    peptides_df.loc[:, 'Parent Antigen Gene Isoform Name'] = peptides_df['Parent Antigen Gene Isoform ID'].map(self.uniprot_id_to_name_map)
    peptides_df.reset_index(inplace=True)
    peptides_df.loc[:, 'ARC Assignment'] = peptides_df['Source Accession'].map(self.source_arc_assignment)

    peptides_df.drop_duplicates(subset=['Source Accession', 'Sequence'], inplace=True) # drop duplicate peptides
    sources_df.drop(columns=['Sequence'], inplace=True) # drop sequence column for output

    self._remove_files()
    
    num_sources = len(sources_df['Accession'].drop_duplicates())
    num_peptides = len(peptides_df[['Source Accession', 'Sequence']].drop_duplicates())
    num_matched_sources = len(sources_df[sources_df['Assigned Protein ID'].notnull()])
    num_matched_peptides = len(peptides_df[peptides_df['Parent Antigen Gene Isoform ID'].notnull()])
  
    assigner_data = (
      num_sources,
      num_peptides,
      num_matched_sources,
      num_matched_peptides
    )

    return assigner_data, peptides_df, sources_df

  def _assign_sources(self, sources_df: pd.DataFrame) -> None:
    """Assign a gene to the peptide sources of a species.

    Run ARC for vertebrates to assign MHC/TCR/Ig to peptide sources first.
    Run BLAST for all other peptide sources to assign a gene and protein. 
    If there are ties, use PEPMatch to search the peptides within the protein 
    sequences andselect the protein with the most matches.

    Args:
      sources_df: DataFrame of peptide sources for a species.
    """    
    self._write_to_fasta(sources_df, 'sources') # write sources to FASTA file

    print('Running BLAST for peptide sources...')
    self._run_blast()

    if self.is_vertebrate:

      if sources_df['Sequence'].isna().all():
        return
      
      print('Running ARC for MHC/TCR/Ig assignments...')
      self._run_arc(sources_df)

  def _create_source_to_peptides_map(self, peptides_df: pd.DataFrame) -> dict:
    """Create a map from peptide sources to their peptides.
    
    Args:
      peptides_df: DataFrame of peptides for a species.
    """    
    source_to_peptides_map = {}
    for i, row in peptides_df.iterrows():
      if row['Source Accession'] in source_to_peptides_map.keys():
        source_to_peptides_map[row['Source Accession']].append(row['Sequence'])
      else:
        source_to_peptides_map[row['Source Accession']] = [row['Sequence']]
    
    return source_to_peptides_map 

  def _write_to_fasta(self, sources_df: pd.DataFrame, filename: str) -> None:
    """Write peptide sources to FASTA file.

    Args:
      sources_df: DataFrame of peptide sources for a species.
      filename: Name of the FASTA file to write to.    
    """          
    seq_records = [] # create seq records of sources with ID and sequence
    for i, row in sources_df.iterrows():

      # if there is no sequence, use empty string
      if pd.isnull(row['Sequence']):
        continue

      seq_records.append(
        SeqRecord(
          Seq(row['Sequence']),
          id=row['Accession'],
          description='')
      )
    with open(f'{self.species_path}/{filename}.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')

  def _preprocess_proteome_if_needed(self) -> None:
    """Preprocess the proteome if the preprocessed files don't exist."""
    if not os.path.exists(f'{self.species_path}/proteome.db'):
      gp_proteome = f'{self.species_path}/gp_proteome.fasta' if os.path.exists(f'{self.species_path}/gp_proteome.fasta') else ''
      Preprocessor(
        proteome = f'{self.species_path}/proteome.fasta',
        preprocessed_files_path = f'{self.species_path}',
        gene_priority_proteome=gp_proteome
      ).sql_proteome(k = 5)

  def _run_blast(self) -> None:
    """BLAST peptide sources against the selected proteome, then read in with
    pandas and assign column names. By default, blastp doesn't return header.

    Then, create a quality score based on % identity, alignment length, and
    query length. Select the best match for each peptide source and assign
    the gene symbol and UniProt ID to the source_gene_assignment and
    source_protein_assignment maps.

    If there are ties, use PEPMatch to search the peptides within the protein
    sequences and select the protein with the most matches.
    """
    # escape parentheses in species path
    species_path = str(self.species_path).replace('(', '\\(').replace(')', '\\)')

    # if proteome.fasta is empty then return
    if os.path.getsize(f'{self.species_path}/proteome.fasta') == 0:
      return
    
    os.system( # make BLAST database from proteome
      f'makeblastdb -in {species_path}/proteome.fasta '\
      f'-dbtype prot > /dev/null'
    )
    os.system( # run blastp
      f'blastp -query {species_path}/sources.fasta '\
      f'-db {species_path}/proteome.fasta '\
      f'-evalue 1  -num_threads {self.num_threads} -outfmt 10 '\
      f'-out {species_path}/blast_results.csv'
    )
    result_columns = [
      'Query', 'Target', 'Percentage Identity', 'Alignment Length', 
      'Mismatches', 'Gap Opens', 'Query Start', 'Query End', 
      'Target Start', 'Target End', 'e-Value', 'Bit Score'
    ]
    blast_results_df = pd.read_csv( # read in BLAST results
      f'{self.species_path}/blast_results.csv', names=result_columns
    )

    # extract the UniProt ID from the target column
    blast_results_df['Target'] = blast_results_df['Target'].str.split('|').str[1]

    # take out "-#" portion of the target UniProt ID because these are isoforms and 
    # won't be mapped properly to gene symbols
    blast_results_df['Target'] = blast_results_df['Target'].str.split('-').str[0]
    
    # map target UniProt IDs to gene symbols
    blast_results_df['Target Gene Symbol'] = blast_results_df['Target'].map(self.uniprot_id_to_gene_map)

    # create a quality score based on % identity, alignment length, and query length
    blast_results_df['Query Length'] = blast_results_df['Query'].astype(str).map(self.source_length_map)
    blast_results_df['Quality Score'] = blast_results_df['Percentage Identity'] * (blast_results_df['Alignment Length'] / blast_results_df['Query Length'])

    # join proteome metadata to select the best match for each source
    blast_results_df = blast_results_df.merge(
      self.proteome[['Protein ID', 'Protein Existence Level', 'Gene Priority']], 
      left_on='Target', 
      right_on='Protein ID', 
      how='left'
    )
    # sort by quality score (descending), gene priority (descending), and
    # protein existence level (ascending)
    blast_results_df.sort_values(
      by=['Quality Score', 'Gene Priority', 'Protein Existence Level'],
      ascending=[False, False, True],
      inplace=True
    )
    # after sorting, drop duplicates based on 'Query', keeping only the first (i.e., best) match.
    blast_results_df.drop_duplicates(subset='Query', keep='first', inplace=True)

    # assign gene symbols, protein ID, and score
    for i, row in blast_results_df.iterrows():
      self.source_gene_assignment[str(row['Query'])] = row['Target Gene Symbol']
      self.source_protein_assignment[str(row['Query'])] = row['Target']
      self.source_assignment_score[str(row['Query'])] = row['Quality Score']

  def _assign_peptides(self, peptides_df: pd.DataFrame) -> None:
    """Assign a parent protein to each peptide.
    
    Preprocess the proteome and then search all the peptides within
    the proteome using PEPMatch. Then, assign the parent protein
    to each peptide by selecting the best isoform of the assigned gene for
    its source.
    """
    self._preprocess_proteome_if_needed()

    # search all peptides within the proteome using PEPMatch
    all_peptides = peptides_df['Sequence'].unique().tolist()
    all_matches_df = Matcher(
      query = all_peptides,
      proteome_file = f'{self.species_path}/proteome.fasta',
      max_mismatches = 0, 
      k = 5, 
      preprocessed_files_path = f'{self.species_path}', 
      best_match=False, 
      output_format='dataframe',
      sequence_version=False
    ).match()
    
    # if no peptide sources were assigned, return
    if not self.source_gene_assignment or not self.source_protein_assignment:
      return

    # create dataframes of peptide source mappings so we can merge and perform operations
    peptide_source_map_df = pd.DataFrame({
      'Peptide': peptide,
      'Source': source
    } for source, peptides in self.source_to_peptides_map.items() for peptide in peptides)
    source_gene_assignment_df = pd.DataFrame({
      'Source': source, 
      'Gene': gene
    } for source, gene in self.source_gene_assignment.items())

    # merge the sources for each peptide
    merged_df = pd.merge( 
      all_matches_df, 
      peptide_source_map_df, 
      how='left', 
      left_on='Query Sequence', 
      right_on='Peptide'
    )

    # merge the protein ID for each source
    merged_df = pd.merge( 
      merged_df, 
      source_gene_assignment_df, 
      how='left', 
      on='Source',
      suffixes=('', '_assigned')
    )

    # isolate only those peptides that match to the assigned gene for their source
    merged_df = merged_df[merged_df['Gene'] == merged_df['Gene_assigned']]

    # now, get the isoform of the assigned gene with the best protein existence level
    best_isoform_indices = merged_df.groupby(['Gene', 'Query Sequence'])['Protein Existence Level'].idxmin()
    best_isoforms = merged_df.loc[best_isoform_indices, ['Gene', 'Query Sequence', 'Protein ID']].set_index(['Gene', 'Query Sequence'])['Protein ID']
    merged_df['Best Isoform ID'] = merged_df.set_index(['Gene', 'Query Sequence']).index.map(best_isoforms)
    
    # drop any unassigned peptides that couldn't be assigned an isoform
    merged_df = merged_df.dropna(subset=['Best Isoform ID'])

    # Update the protein assignment with the best isoform
    self.peptide_protein_assignment.update(
      dict(zip(
        zip(merged_df['Source'], merged_df['Query Sequence']),
        merged_df['Best Isoform ID']))
    )

  def _run_arc(self, sources_df: pd.DataFrame) -> None:
    """Run ARC to assign MHC/TCR/Ig to peptide sources.

    Args:
        sources_df: DataFrame of peptide sources for a species.
    """
    try:
      past_arc_results_df = pd.read_csv(f'{self.species_path}/ARC_results.tsv', sep='\t')
      past_arc_ids = set(past_arc_results_df['id'])
      arc_sources_df = sources_df[~sources_df['Accession'].isin(past_arc_ids)]
      past_results = True
    except FileNotFoundError:
      arc_sources_df = sources_df
      past_results = False

    if not arc_sources_df.empty:
      self._write_to_fasta(arc_sources_df, 'ARC_sources')

      if os.path.getsize(f'{self.species_path}/ARC_sources.fasta') != 0: # make sure file isn't empty
        SeqClassifier(
          outfile=f'{self.species_path}/ARC_temp_results.tsv',
          threads=self.num_threads,
        ).classify_seqfile(f'{self.species_path}/ARC_sources.fasta')

      new_arc_results_df = pd.read_csv(f'{self.species_path}/ARC_temp_results.tsv', sep='\t')

      if past_results:
        combined_results_df = pd.concat([past_arc_results_df, new_arc_results_df])
      else:
        combined_results_df = new_arc_results_df

      combined_results_df.to_csv(f'{self.species_path}/ARC_results.tsv', sep='\t', index=False)
      
      os.remove(f'{self.species_path}/ARC_temp_results.tsv')
      os.remove(f'{self.species_path}/ARC_sources.fasta')

    # Update source to ARC assignment dictionary
    arc_results_df = pd.read_csv(f'{self.species_path}/ARC_results.tsv', sep='\t')
    if not arc_results_df.dropna(subset=['class']).empty:
      self.source_arc_assignment = arc_results_df.set_index('id')['class'].to_dict()

  def _assign_allergens(self) -> None:
    """Get allergen data from allergen.org and then assign allergens to sources."""
    allergen_df = pd.read_csv(self.build_path / 'arborist' / 'allergens.tsv', sep='\t')
    allergen_map = allergen_df.set_index('AccProtein')['Name'].to_dict()

    for k, v in self.source_protein_assignment.items():
      if v in allergen_map.keys():
        self.uniprot_id_to_name_map[v] = allergen_map[v]

  def _assign_manuals(self) -> None:
    """Get manual assignments from manual-parents.tsv and then assign
    genes and proteins to sources.
    """
    manual_df = pd.read_csv(self.build_path / 'arborist' / 'manual-parents.tsv', sep='\t')

    manual_gene_map = manual_df.set_index('Accession')['Accession Gene'].to_dict()
    manual_protein_id_map = manual_df.set_index('Accession')['Parent Accession'].to_dict()
    manual_protein_name_map = manual_df.set_index('Accession')['Parent Name'].to_dict()
    
    for k, v in self.source_gene_assignment.items():
      if k in manual_gene_map.keys():
        self.source_gene_assignment[k] = manual_gene_map[k]
        self.source_protein_assignment[k] = manual_protein_id_map[k]
        self.source_assignment_score[k] = -1
        self.uniprot_id_to_name_map[k] = manual_protein_name_map[k]

  def _remove_files(self) -> None:
    """Delete all the files that were created when making the BLAST database."""
    for extension in ['pdb', 'phr', 'pin', 'pjs', 'pot', 'psq', 'ptf', 'pto']:
      try: # if DB wasn't created this will throw an error
        os.remove(glob.glob(f'{self.species_path}/*.{extension}')[0])
      except IndexError:
        pass 

    (self.species_path / 'blast_results.csv').unlink(missing_ok=True)
    (self.species_path / 'sources.fasta').unlink(missing_ok=True)

def run(taxon_id, species_name, group, all_taxa, build_path, all_peptides, all_sources, num_threads):
  species_path = build_path / 'species' / f'{taxon_id}' # directory to write species data

  peptides_df = DataFetcher(build_path).get_peptides_for_species(all_peptides, all_taxa)
  sources_df = DataFetcher(build_path).get_sources_for_species(all_sources, peptides_df['Source Accession'].tolist())

  if sources_df.empty or peptides_df.empty:
    return

  is_vertebrate = group == 'vertebrate'

  print(f'Assigning peptides and sources for species w/ taxon ID: {taxon_id}).')
  Assigner = GeneAndProteinAssigner(
    taxon_id,
    species_path,
    is_vertebrate,
    group,
    peptides_df,
    num_threads=num_threads,
    build_path=build_path,
  )

  assigner_data, peptide_assignments, source_assignments = Assigner.assign(sources_df, peptides_df)
  num_sources, num_peptides, num_matched_sources, num_matched_peptides = assigner_data

  peptide_assignments.to_csv(species_path / 'peptide-assignments.tsv', sep='\t', index=False)
  source_assignments.to_csv(species_path / 'source-assignments.tsv', sep='\t', index=False)

  assignment_data = (
    taxon_id,
    species_name,
    num_sources,
    num_peptides,
    num_matched_sources,
    num_matched_peptides
  )

  species_data = pd.read_csv(species_path / 'species-data.tsv', sep='\t')
  species_data.loc[:, 'Num Sources'] = num_sources
  species_data.loc[:, 'Num Peptides'] = num_peptides
  species_data.loc[:, '% Assigned Sources'] = round((num_matched_sources / num_sources)*100, 2)
  species_data.loc[:, '% Assigned Peptides'] = round((num_matched_peptides / num_peptides)*100, 2)
  species_data.to_csv(species_path / 'species-data.tsv', sep='\t', index=False)

  return assignment_data

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
    '-n', '--num_threads', 
    type=int,
    default=os.cpu_count() - 1,
    help='Number of threads to use.'
  )

  args = parser.parse_args()

  build_path = Path(args.build_path)
  all_species = args.all_species
  taxon_id = args.taxon_id
  num_threads = args.num_threads

  assert all_species or taxon_id, 'Please specify either --all_species or --taxon_id.'
  assert (all_species and not taxon_id) or (taxon_id and not all_species), 'Please specify either --all_species or --taxon_id.'

  # TODO: replace the data/active-species.tsv with updated arborist active-species somehow
  species_df = pd.read_csv(build_path / 'arborist' / 'active-species.tsv', sep='\t')
  valid_taxon_ids = species_df['Species ID'].tolist()

  all_peptides = DataFetcher(build_path).get_all_peptides()
  all_sources = DataFetcher(build_path).get_all_sources()

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
      assignment_data = run(
        taxon_id, species_name, group, all_taxa, build_path, all_peptides, all_sources, num_threads
      )

  else: # one species at a time
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    species_name = taxon_to_species_map[taxon_id]
    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]
    group = species_df[species_df['Species ID'] == taxon_id]['Group'].iloc[0]
    assignment_data = run(
      taxon_id, species_name, group, all_taxa, build_path, all_peptides, all_sources, num_threads
    )