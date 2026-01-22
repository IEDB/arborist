import argparse
import json
import polars as pl
from pathlib import Path

def create_species_map(allergens, active_species):
  iuis_to_clean = {}
  for s in allergens['Species'].unique().to_list():
    clean = s.lower().split('(')[0].strip()
    iuis_to_clean[s] = clean

  species_map = {}
  for row in active_species.iter_rows(named=True):
    taxon_id = row['Species ID']
    label = row['Species Label'].lower()
    for iuis_full, iuis_clean in iuis_to_clean.items():
      if iuis_clean in label or label.split('(')[0].strip() in iuis_clean:
        species_map[iuis_full] = taxon_id
        break

  manual_fixes = {
    'Bombyx mori': 7091,
    'Dermatophagoides pteronyssinus': 6956,
    'Brassica napus': 3708,
    'Prunus persica': 3760,
  }
  species_map.update(manual_fixes)
  species_map.pop('Exopalaemon modestus', None)

  return species_map

def main(build_path):
  allergens = pl.read_csv(build_path / 'arborist' / 'allergens.tsv', separator='\t')
  active_species = pl.read_csv(build_path / 'arborist' / 'active-species.tsv', separator='\t')

  species_map = create_species_map(allergens, active_species)
  map_output = build_path / 'arborist' / 'allergen-species-map.json'

  with open(map_output, 'w') as f:
    json.dump(species_map, f, indent=2, sort_keys=True)
  print(f"Created {map_output} with {len(species_map)} species mappings")

  allergens_with_taxa = allergens.with_columns(
    pl.col('Species').map_elements(lambda s: species_map.get(s), return_dtype=pl.Int64).alias('SpeciesID')
  )

  tsv_output = build_path / 'arborist' / 'allergens.tsv'
  allergens_with_taxa.write_csv(tsv_output, separator='\t')
  print(f"Updated {tsv_output} with taxon_id column")

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-b', '--build_path', type=str,
    default=Path(__file__).parents[3] / 'build'
  )
  args = parser.parse_args()
  main(Path(args.build_path))
