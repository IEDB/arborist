import polars as pl
from pathlib import Path


def calculate_source_diffs():
  previous_sources = pl.read_csv(proteins_path / 'previous' / 'source-parents.tsv', separator='\t')
  latest_sources = pl.read_csv(proteins_path / 'latest' / 'source-parents.tsv', separator='\t')

  previous_sources = previous_sources.filter(pl.col('Status') == 'Active')

  previous_counts = previous_sources.group_by(['Species ID', 'Species Label']).agg(
    pl.count('Source ID').alias('Previous Count')
  )
  latest_counts = latest_sources.group_by(['Species ID', 'Species Label']).agg(
    pl.count('Source ID').alias('Latest Count')
  )

  diffs = latest_counts.join(
    previous_counts, on='Species ID', how='left', coalesce=True
  ).fill_null(0).select(
    pl.col('Species ID', 'Species Label', 'Previous Count', 'Latest Count')
  ).with_columns(
    (pl.col('Latest Count') - pl.col('Previous Count')).alias('Diff')
  )
  
  total = diffs.select(pl.sum('Previous Count'), pl.sum('Latest Count'), pl.sum('Diff')).with_columns(
    pl.lit('Total').alias('Species Label'), pl.lit(0).cast(pl.Int64).alias('Species ID')
  ).select(pl.col('Species ID', 'Species Label', 'Previous Count', 'Latest Count', 'Diff'))

  diffs = pl.concat([diffs, total]).sort('Previous Count', descending=True)
  diffs.write_csv(proteins_path / 'diffs' / 'source-diffs.tsv', separator='\t')


def calculate_parent_diffs():
  previous_parents = pl.read_csv(proteins_path / 'previous' / 'parent-proteins.tsv', separator='\t')
  latest_parents = pl.read_csv(proteins_path / 'latest' / 'parent-proteins.tsv', separator='\t')
  previous_sources = pl.read_csv(proteins_path / 'previous' / 'source-parents.tsv', separator='\t')
  latest_sources = pl.read_csv(proteins_path / 'latest' / 'source-parents.tsv', separator='\t')


  previous_counts = previous_parents.join(
    previous_sources, how='left', left_on='Accession', right_on='Parent Protein Accession', coalesce=True
  ).unique(subset='Accession').group_by(['Species ID', 'Species Label']).agg(
    pl.count('Accession').alias('Previous Count')
  )

  latest_counts = latest_parents.join(
    latest_sources, how='left', left_on='Accession', right_on='Parent Protein Accession', coalesce=True
  ).unique(subset='Accession').group_by(['Species ID', 'Species Label']).agg(
    pl.count('Accession').alias('Latest Count')
  )

  diffs = latest_counts.join(
    previous_counts, on='Species ID', how='left', coalesce=True
  ).fill_null(0).select(
    pl.col('Species ID', 'Species Label', 'Previous Count', 'Latest Count')
  ).with_columns(
    (pl.col('Latest Count') - pl.col('Previous Count')).alias('Diff')
  )

  total = diffs.select(pl.sum('Previous Count'), pl.sum('Latest Count'), pl.sum('Diff')).with_columns(
    pl.lit('Total').alias('Species Label'), pl.lit(0).cast(pl.Int64).alias('Species ID')
  ).select(pl.col('Species ID', 'Species Label', 'Previous Count', 'Latest Count', 'Diff'))

  diffs = pl.concat([diffs, total]).sort('Previous Count', descending=True)
  diffs.write_csv(proteins_path / 'diffs' / 'parent-diffs.tsv', separator='\t')


def calculate_epitope_diffs():
  previous_epitopes = pl.read_csv(proteins_path / 'previous' / 'epitope-mappings.tsv', separator='\t')
  latest_epitopes = pl.read_csv(proteins_path / 'latest' / 'epitope-mappings.tsv', separator='\t')
  previous_sources = pl.read_csv(proteins_path / 'previous' / 'source-parents.tsv', separator='\t')
  latest_sources = pl.read_csv(proteins_path / 'latest' / 'source-parents.tsv', separator='\t')

  previous_counts = previous_epitopes.join(
    previous_sources, how='left', left_on='source_accession', right_on='Accession', coalesce=True
  ).unique(subset='epitope_id').group_by(['Species ID', 'Species Label']).agg(
    pl.count('epitope_id').alias('Previous Count')
  )

  latest_counts = latest_epitopes.join(
    latest_sources, how='left', left_on='source_accession', right_on='Accession', coalesce=True
  ).unique(subset='epitope_id').group_by(['Species ID', 'Species Label']).agg(
    pl.count('epitope_id').alias('Latest Count')
  )

  diffs = latest_counts.join(
    previous_counts, on='Species ID', how='left', coalesce=True
  ).fill_null(0).select(
    pl.col('Species ID', 'Species Label', 'Previous Count', 'Latest Count')
  ).with_columns(
    (pl.col('Latest Count') - pl.col('Previous Count')).alias('Diff')
  )

  total = diffs.select(pl.sum('Previous Count'), pl.sum('Latest Count'), pl.sum('Diff')).with_columns(
    pl.lit('Total').alias('Species Label'), pl.lit(0).cast(pl.Int64).alias('Species ID')
  ).select(pl.col('Species ID', 'Species Label', 'Previous Count', 'Latest Count', 'Diff'))

  diffs = pl.concat([diffs, total]).sort('Previous Count', descending=True)
  diffs.write_csv(proteins_path / 'diffs' / 'epitope-diffs.tsv', separator='\t')


if __name__ == '__main__':

  proteins_path = Path(__file__).parents[2] / 'build' / 'proteins'
  (proteins_path / 'diffs').mkdir(parents=True, exist_ok=True)

  calculate_source_diffs()
  calculate_parent_diffs()
  calculate_epitope_diffs()