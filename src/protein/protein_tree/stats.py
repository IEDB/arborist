#!/usr/bin/env python3
"""Capture week-to-week parent-assignment stats into one master sheet.

Runs at the END of `make protein`. Reads only build/arborist/all-peptide-assignments.tsv
(every per-species metric is derivable there: it carries Species Taxon ID, Source
Accession, ARC Assignment, Epitope ID, Assigned Protein ID and the exact-match
position). Never touches assignment logic.

Metrics per species (and a TOTAL = column sum):
  source_antigens_total     unique Source Accession
  source_antigens_assigned  unique Source Accession with a protein hit OR an ARC class
  epitopes_total            unique Epitope ID
  epitopes_assigned         unique Epitope ID with a parent (protein id OR ARC) -- "has a home"
  epitopes_located_exact    unique Epitope ID found EXACTLY in the parent/isoform (has a
                            match position)

The master sheet (repo dir, git-tracked, survives `make clean`) is wide: rows are
(Species ID, Species Name, Metric), one column per ISO week. Re-running the same week
overwrites that week's column (idempotent); species entering/leaving are handled by an
outer pivot.
"""

import argparse
import sys
import datetime
from pathlib import Path

import polars as pl

REPO_ROOT = Path(__file__).parents[3]

METRICS = [
  'source_antigens_total',
  'source_antigens_assigned',
  'epitopes_total',
  'epitopes_assigned',
  'epitopes_located_exact',
]
ID_COLUMNS = ['Species ID', 'Species Name', 'Metric']

# A source/epitope counts as "assigned" if ARC classified it as a receptor/MHC node.
_ARC_RECEPTOR = r'TCR|BCR|MHC'


def iso_week(day: datetime.date = None) -> str:
  day = day or datetime.date.today()
  year, week, _ = day.isocalendar()
  return f'{year}-W{week:02d}'


def _nonempty(column: str) -> pl.Expr:
  return pl.col(column).is_not_null() & (pl.col(column) != '')


def compute_week_stats(assignments: pl.DataFrame) -> pl.DataFrame:
  """Per-species metrics as a long frame [Species ID, Species Name, Metric, value] + TOTAL."""
  arc_receptor = pl.col('ARC Assignment').fill_null('').str.contains(_ARC_RECEPTOR)
  source_assigned = _nonempty('Source Assigned Protein ID') | arc_receptor
  epitope_assigned = _nonempty('Assigned Protein ID') | arc_receptor
  epitope_located = _nonempty('Assigned Protein Starting Position')

  per_species = assignments.group_by(['Species Taxon ID', 'Species Name']).agg(
    pl.col('Source Accession').n_unique().alias('source_antigens_total'),
    pl.col('Source Accession').filter(source_assigned).n_unique().alias('source_antigens_assigned'),
    pl.col('Epitope ID').n_unique().alias('epitopes_total'),
    pl.col('Epitope ID').filter(epitope_assigned).n_unique().alias('epitopes_assigned'),
    pl.col('Epitope ID').filter(epitope_located).n_unique().alias('epitopes_located_exact'),
  ).rename({'Species Taxon ID': 'Species ID'})

  long = per_species.unpivot(
    index=['Species ID', 'Species Name'], on=METRICS,
    variable_name='Metric', value_name='value',
  ).with_columns(
    pl.col('Species ID').cast(pl.String), pl.col('value').cast(pl.Int64)
  )

  totals = pl.DataFrame({
    'Species ID': ['TOTAL'] * len(METRICS),
    'Species Name': [''] * len(METRICS),
    'Metric': METRICS,
    'value': [int(per_species[m].sum()) for m in METRICS],
  }, schema={'Species ID': pl.String, 'Species Name': pl.String, 'Metric': pl.String, 'value': pl.Int64})

  return pl.concat([long.select(['Species ID', 'Species Name', 'Metric', 'value']), totals])


def _week_columns(sheet: pl.DataFrame) -> list:
  return [c for c in sheet.columns if c not in ID_COLUMNS]


def _sort_sheet(wide: pl.DataFrame) -> pl.DataFrame:
  metric_rank = {m: i for i, m in enumerate(METRICS)}
  weeks = sorted(_week_columns(wide))
  return wide.with_columns(
    pl.when(pl.col('Species ID') == 'TOTAL').then(0).otherwise(1).alias('_grp'),
    pl.col('Species ID').cast(pl.Int64, strict=False).alias('_sid'),
    pl.col('Metric').replace_strict(metric_rank, default=99, return_dtype=pl.Int64).alias('_mr'),
  ).sort(['_grp', '_sid', '_mr']).drop('_grp', '_sid', '_mr').select(ID_COLUMNS + weeks)


def upsert_stats(existing: pl.DataFrame | None, week: str, week_long: pl.DataFrame) -> pl.DataFrame:
  """Insert/replace `week`'s column into the wide master sheet."""
  this_week = week_long.with_columns(pl.lit(week).alias('week')).select(
    ['Species ID', 'Species Name', 'Metric', 'week', 'value']
  )

  if existing is not None and existing.height > 0:
    old_long = existing.unpivot(
      index=ID_COLUMNS, on=_week_columns(existing),
      variable_name='week', value_name='value',
    ).with_columns(pl.col('value').cast(pl.Int64, strict=False))
    old_long = old_long.filter(pl.col('week') != week)
    combined = pl.concat([old_long, this_week])
  else:
    combined = this_week

  wide = combined.pivot(on='week', index=ID_COLUMNS, values='value')
  return _sort_sheet(wide)


def load_sheet(path: Path) -> pl.DataFrame | None:
  if not path.exists():
    return None
  return pl.read_csv(path, separator='\t', infer_schema_length=0)


def run(args) -> None:
  build_path = Path(args.build_path)
  stats_file = Path(args.stats_file)

  assignments = pl.read_csv(
    build_path / 'arborist' / 'all-peptide-assignments.tsv', separator='\t', infer_schema_length=0
  )
  week_long = compute_week_stats(assignments)
  sheet = upsert_stats(load_sheet(stats_file), args.week, week_long)
  sheet.write_csv(stats_file, separator='\t')

  total = {m: v for _, m, v in week_long.filter(pl.col('Species ID') == 'TOTAL')
           .select('Species ID', 'Metric', 'value').iter_rows()}
  print(f"[{args.week}] epitopes {total.get('epitopes_assigned', 0)}/{total.get('epitopes_total', 0)} "
        f"assigned, sources {total.get('source_antigens_assigned', 0)}/"
        f"{total.get('source_antigens_total', 0)} -> {stats_file}")


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('-b', '--build_path', type=str, default=str(REPO_ROOT / 'build'),
                      help='Path for all Arborist build files.')
  parser.add_argument('-o', '--stats_file', type=str,
                      default=str(REPO_ROOT / 'protein-tree-weekly-stats.tsv'),
                      help='Master weekly stats sheet (wide, one column per ISO week).')
  parser.add_argument('-w', '--week', type=str, default=iso_week(),
                      help='ISO week label for this run (default: current ISO week).')
  args = parser.parse_args()

  # Observability must never break the weekly run: log and exit 0 on any error.
  try:
    run(args)
  except Exception as error:  # noqa: BLE001 - deliberately non-breaking
    print(f'stats: non-fatal error, weekly stats not updated: {error}', file=sys.stderr)
  return 0


if __name__ == '__main__':
  sys.exit(main())
