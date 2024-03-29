#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

from get_data import DataFetcher
from select_proteome import ProteomeSelector
from protein_tree.assign import GeneAndProteinAssigner


def run_protein_tree(
  build_path: Path,
  taxon_id: int,
  species_name_map: dict,
  species_group_map: dict,
  epitopes_df: pd.DataFrame,
  sources_df: pd.DataFrame,
  update_proteome: bool,
  num_threads: int,
  force: bool
) -> None:
  """Run all steps for the protein tree.
  
  Args:
    build_path: Path to the build directory.
    taxon_id: Taxon ID for the species to run protein tree.
    species_name_map: Mapping of taxon ID to species name.
    species_group_map: Mapping of taxon ID to group. Groups:
      bacterium, archeobacterium, plant, vertebrate, virus, other-eukaryote, other.
    epitopes_df: Epitope data for the species.
    sources_df: Source antigen data for the species.
    update_proteome: Whether or not to update the proteome to be used for the species.
    num_threads: Number of threads to use for BLAST and ARC.
  """
  species_path = build_path / 'species' / f'{taxon_id}'
  species_name = species_name_map[taxon_id]
  group = species_group_map[taxon_id]

  print(f'Building tree for {species_name} (ID: {taxon_id})...\n')

  if epitopes_df.empty or sources_df.empty:
    print('No epitopes found or source antigens found for this species.')
    return
  
  if update_proteome or not (species_path / 'proteome.fasta').exists():
    update_proteome = True # if the file doesn't exist, update flag
    
    print('Getting the best proteome...')
    Selector = ProteomeSelector(taxon_id, species_path, build_path)
    print(f'Number of candidate proteomes: {Selector.num_of_proteomes}\n')

    proteome_data = Selector.select_best_proteome(epitopes_df)
    Selector.proteome_to_tsv()
    
    print('Got the best proteome:')
    print(f'Proteome ID: {proteome_data[0]}')
    print(f'Proteome taxon: {proteome_data[1]:.0f}')
    print(f'Proteome type: {proteome_data[2]}\n')
  
  else: 
    try: # check if epitopes and sources have changed since last run
      if force:
        raise FileNotFoundError # force run 

      previous_epitopes_df = pd.read_csv(
        species_path / 'epitope_assignments.tsv', sep='\t'
      )
      previous_sources_df = pd.read_csv(
        species_path / 'source_assignments.tsv', sep='\t'
      )
      
      previous_epitopes = set(previous_epitopes_df[['Sequence', 'Accession']].apply(
        lambda row: (row['Sequence'], row['Accession']), axis=1
      ).tolist())
      previous_sources = set(previous_sources_df['Accession'].apply(
        lambda accession: str(accession)
      ).tolist())

      epitopes = set(epitopes_df[['Sequence', 'Accession']].apply(
        lambda row: (row['Sequence'], row['Accession']), axis=1
      ).tolist())
      sources = set(sources_df['Accession'].apply(
        lambda accession: str(accession)
      ).tolist())

      if previous_epitopes == epitopes or previous_sources == sources:
        print('Sources, epitopes, and proteome have not changed since last run.\n')
        return
    
    except FileNotFoundError:
      pass

  # write raw data to files
  epitopes_df.to_csv(species_path / 'epitopes.tsv', sep='\t', index=False)
  sources_df.to_csv(species_path / 'sources.tsv', sep='\t', index=False)

  # boolean for vertebrate species - to be used for ARC assignments
  is_vertebrate = group == 'vertebrate'

  # assign genes to source antigens and parent proteins to epitopes
  Assigner = GeneAndProteinAssigner(
    taxon_id,
    species_path,
    is_vertebrate,
    num_threads,
    build_path=build_path,
    bin_path=build_path.parent.parent / 'bin'
  )
  assigner_data, epitope_assignments, source_assignments = Assigner.assign(sources_df, epitopes_df)

  # write assignments to files
  epitope_assignments.to_csv(
    species_path / 'epitope_assignments.tsv', sep='\t', index=False
  )
  source_assignments.to_csv(
    species_path / 'source_assignments.tsv', sep='\t', index=False
  )

  successful_source_assignment = (assigner_data[2] / assigner_data[0])*100
  successful_epitope_assignment = (assigner_data[3] / assigner_data[1])*100

  print(f'Number of sources: {assigner_data[0]}')
  print(f'Number of epitopes: {assigner_data[1]}')
  print(f'Successful source antigen assignments: {successful_source_assignment:.1f}%')
  print(f'Successful epitope assignments: {successful_epitope_assignment:.1f}%\n')

  print(f'Protein tree built for {species_name} (ID: {taxon_id}).\n\n')
  

def main():
  import argparse
  import multiprocessing

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
    help='Taxon ID for the species to run protein tree individually.'
  )
  parser.add_argument(
    '-d', '--update_data',
    action='store_true',
    help='Pull the epitope and source tables from the IEDB backend.'
  )
  parser.add_argument(
    '-p', '--update_proteome',
    action='store_true',
    help='Update the proteome(s) to be used for the species.'
  )
  parser.add_argument(
    '-s', '--update_species',
    action='store_true',
    help='Update the species table.'
  )
  parser.add_argument(
    '-n', '--num_threads',
    type=int,
    default=multiprocessing.cpu_count() - 2,
    help='Number of threads to use for BLAST and ARC.'
  )
  parser.add_argument(
    '-f', '--force',
    action='store_true',
    help='Force run the protein tree for a species, even if the data has not changed.'
  )

  args = parser.parse_args()

  build_path = Path(args.build_path)

  if args.update_species:
    print('Updating species table...')
    DataFetcher.update_species() # call update_species() static method
    print('Done.\n')

  species_df = pd.read_csv(build_path / 'arborist' / 'active-species.tsv', sep='\t')
  all_species_taxa = species_df['Species ID'].tolist()
  
  # key, taxa, species name, and group mapppings
  species_name_map = dict(
    zip(
      species_df['Species ID'],
      species_df['Species Label']
    )
  )
  all_taxa_map = dict(
    zip(
      species_df['Species ID'],
      species_df['Active Taxa']
    )
  )
  species_group_map = dict(
    zip(
      species_df['Species ID'],
      species_df['Group']
    )
  )

  Fetcher = DataFetcher(build_path)

  files_exist = (
    (build_path / 'iedb' / 'peptides.tsv').exists() and
    (build_path / 'iedb' / 'sources.tsv').exists() and
    (build_path / 'arborist' / 'allergens.tsv').exists()
  )
  if args.update_data or not files_exist:
    print('Getting all data...')
    Fetcher.get_all_data()
    print('All data written.')

  all_epitopes = Fetcher.get_all_epitopes()
  all_sources = Fetcher.get_all_sources()

  if args.all_species:
    for taxon_id in all_species_taxa:

      all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]

      epitopes_df = Fetcher.get_epitopes_for_species(all_epitopes, all_taxa)
      sources_df = Fetcher.get_sources_for_species(
        all_sources, epitopes_df['Accession'].tolist()
      )

      run_protein_tree(
        build_path, taxon_id, species_name_map, species_group_map, epitopes_df,
        sources_df, args.update_proteome, args.num_threads, args.force
      )

    print('All species complete.')

  else: # one species at a time

    taxon_id = args.taxon_id
    assert taxon_id in all_species_taxa, f'{taxon_id} is not a valid taxon ID.'

    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]

    epitopes_df = Fetcher.get_epitopes_for_species(all_epitopes, all_taxa)
    sources_df = Fetcher.get_sources_for_species(
      all_sources, epitopes_df['Accession'].tolist()
    )
    
    run_protein_tree(
      build_path, taxon_id, species_name_map, species_group_map, epitopes_df,
      sources_df, args.update_proteome, args.num_threads, args.force
    )

if __name__ == '__main__':
  main()
