#!/usr/bin/env python3

import os
import pandas as pd

def create_proteome_link(row):
  if 'Orphans' in row['Proteome Type']:
    return f"https://www.uniprot.org/uniprotkb?query=taxonomy_id:{row['Species Taxon ID']}"
  else:
    return f"https://www.uniprot.org/uniprotkb?query=proteome:{row['Proteome ID']}"

if __name__ == '__main__':
  df = pd.DataFrame()
  for path, _, files in os.walk('build/species'):
    for name in files:
      if name == 'species-data.tsv':
        species_data_df = pd.read_csv(os.path.join(path, name), sep='\t')
        species_data_df['Species Taxon ID'] = path.split('/')[-1]
        df = pd.concat([df, species_data_df])

  df.reset_index(inplace=True, drop=True)
  df.insert(0, 'Species Taxon ID', df.pop('Species Taxon ID')) # move taxon ID to first
  df['Proteome Link'] = df.apply(create_proteome_link, axis=1)
  df.sort_values(by=['Species Taxon ID'], inplace=True)
  df.to_csv('build/species-data.tsv', sep='\t', index=False)