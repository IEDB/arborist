import polars as pl
from pathlib import Path

species_map = {
  108: {'taxa': [108], 'sources': ['WP_012925058.1']},
  9593: {'taxa': [9593, 9595], 'sources': ['Q9BDE0.2', 'CAA42809.1']},
  9606: {'taxa': [9606], 'sources': ['P01116.1', 'AAH13572.1']},
  3052236: {'taxa': [2847087], 'sources': ['Q69422.1']}
}

prod_build_path = Path('../build/')
test_build_path = Path('build/')

# write active species file
active_species_df = pl.read_csv(prod_build_path / 'arborist' / 'active-species.tsv', separator='\t')
active_species_df.filter(
  pl.col('Species ID').is_in(species_map.keys())
).write_csv(test_build_path / 'arborist' / 'active-species.tsv', separator='\t')

sources = []
for species in species_map.values():
  sources += species['sources']

# write peptide and peptide source files
peptide_df = pl.read_csv(prod_build_path / 'iedb' / 'peptide.tsv', separator='\t')
peptide_source_df = pl.read_csv(prod_build_path / 'iedb' / 'peptide_source.tsv', separator='\t')

peptide_df.filter(
  pl.col('Source Accession').is_in(sources)
).write_csv(test_build_path / 'iedb' / 'peptide.tsv', separator='\t')

peptide_source_df.filter(
  pl.col('Source Accession').is_in(sources)
).write_csv(test_build_path / 'iedb' / 'peptide_source.tsv', separator='\t')
