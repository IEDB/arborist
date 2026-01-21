import argparse
import json
import polars as pl
from pathlib import Path

def main():
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

  output_file = Path('src/protein/data/allergen_species_map.json')
  with open(output_file, 'w') as f:
    json.dump(species_map, f, indent=2, sort_keys=True)

  print(f"Saved {len(species_map)} species mappings to {output_file}")

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-b', '--build_path', type=str, help='Path for all Arborist build files.',
    default=Path(__file__).parents[3] / 'build'
  )
  args = parser.parse_args()

  build_path = Path(args.build_path)
  allergens = pl.read_csv(build_path / 'arborist' / 'allergens.csv')
  active_species = pl.read_csv(build_path / 'arborist' / 'active-species.tsv', separator='\t')
  
  main()
