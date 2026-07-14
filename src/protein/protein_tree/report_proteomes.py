#!/usr/bin/env python3
"""Report the proteome selected for each active species.

Runs AFTER `make proteome`. It only scans the filesystem produced by
select_proteome.py; it never triggers or changes proteome selection. The goal is
to surface, in one place, which active species got a proteome and which came back
EMPTY or MISSING (a silent cause of dropped epitopes downstream).

Status per active species:
  SELECTED - proteome.fasta has content and a real proteome was chosen
  ORPHANS  - proteome.fasta has content but was built from orphan proteins
  EMPTY    - a species dir exists but proteome.fasta is missing or 0 bytes
             (Reason: no-candidates if the candidate list was empty, else fetch-empty)
  MISSING  - no species dir at all for an active species
"""

import argparse
import sys
from pathlib import Path

import polars as pl

STATUS_SELECTED = 'SELECTED'
STATUS_ORPHANS = 'ORPHANS'
STATUS_EMPTY = 'EMPTY'
STATUS_MISSING = 'MISSING'

REPORT_COLUMNS = [
  'Species ID', 'Species Name', 'Group', 'Status', 'Proteome ID', 'Proteome Type',
  'Proteome Label', 'Candidate Count', 'FASTA Bytes', 'Reason',
]
REPORT_SCHEMA = {
  'Species ID': pl.String, 'Species Name': pl.String, 'Group': pl.String,
  'Status': pl.String, 'Proteome ID': pl.String, 'Proteome Type': pl.String,
  'Proteome Label': pl.String, 'Candidate Count': pl.Int64, 'FASTA Bytes': pl.Int64,
  'Reason': pl.String,
}


def _read_species_metadata(species_path: Path):
  """Return the single species-data.tsv row as a dict, or None."""
  metadata_file = species_path / 'species-data.tsv'
  if not metadata_file.exists():
    return None
  df = pl.read_csv(metadata_file, separator='\t', infer_schema_length=0)
  if df.height == 0:
    return None
  return df.row(0, named=True)


def _candidate_count(species_path: Path):
  """Number of candidate proteomes considered, or None if the list is absent."""
  list_file = species_path / 'proteome-list.tsv'
  if not list_file.exists():
    return None
  try:
    return pl.read_csv(list_file, separator='\t', infer_schema_length=0).height
  except pl.exceptions.NoDataError:
    return 0


def classify_species(species_path: Path):
  """Classify one species dir into a report row (values only).

  Returns a dict keyed by REPORT_COLUMNS minus Species ID/Name/Group (which come
  from active-species.tsv).
  """
  if not species_path.exists():
    return {
      'Status': STATUS_MISSING, 'Proteome ID': '', 'Proteome Type': '',
      'Proteome Label': '', 'Candidate Count': None, 'FASTA Bytes': 0,
      'Reason': 'no species dir',
    }

  fasta = species_path / 'proteome.fasta'
  fasta_bytes = fasta.stat().st_size if fasta.exists() else 0
  candidate_count = _candidate_count(species_path)
  metadata = _read_species_metadata(species_path)
  proteome_id = str((metadata or {}).get('Proteome ID', '') or '')
  proteome_type = str((metadata or {}).get('Proteome Type', '') or '')
  proteome_label = str((metadata or {}).get('Proteome Label', '') or '')

  if fasta_bytes > 0:
    if proteome_type.strip().lower() == 'orphans':
      status, reason = STATUS_ORPHANS, 'ok'
    else:
      status = STATUS_SELECTED
      reason = 'ok' if metadata is not None else 'no-metadata'
  else:
    status = STATUS_EMPTY
    reason = 'no-candidates' if not candidate_count else 'fetch-empty'

  return {
    'Status': status, 'Proteome ID': proteome_id, 'Proteome Type': proteome_type,
    'Proteome Label': proteome_label, 'Candidate Count': candidate_count,
    'FASTA Bytes': fasta_bytes, 'Reason': reason,
  }


def build_report(build_path: Path) -> pl.DataFrame:
  """Build the proteome-selection report for every active species."""
  build_path = Path(build_path)
  active_species = pl.read_csv(
    build_path / 'arborist' / 'active-species.tsv', separator='\t', infer_schema_length=0
  )
  rows = []
  for species in active_species.rows(named=True):
    taxon_id = str(species['Species ID'])
    species_path = build_path / 'species' / taxon_id
    row = {
      'Species ID': taxon_id,
      'Species Name': species.get('Species Label', '') or '',
      'Group': species.get('Group', '') or '',
    }
    row.update(classify_species(species_path))
    rows.append(row)

  return pl.DataFrame(rows, schema=REPORT_SCHEMA).select(REPORT_COLUMNS)


def summarize(report: pl.DataFrame) -> str:
  """One-line stderr summary of the status counts."""
  counts = dict(
    report.group_by('Status').len().iter_rows()
  )
  return (
    f"{counts.get(STATUS_SELECTED, 0)} selected, {counts.get(STATUS_ORPHANS, 0)} orphans, "
    f"{counts.get(STATUS_EMPTY, 0)} EMPTY, {counts.get(STATUS_MISSING, 0)} missing"
  )


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    '-b', '--build_path', type=str, help='Path for all Arborist build files.',
    default=Path(__file__).parents[3] / 'build',
  )
  args = parser.parse_args()
  build_path = Path(args.build_path)

  report = build_report(build_path)
  output = build_path / 'arborist' / 'proteome-selection-report.tsv'
  report.write_csv(output, separator='\t')

  print(f'Proteome selection: {summarize(report)}', file=sys.stderr)
  print(f'Wrote {output}', file=sys.stderr)


if __name__ == '__main__':
  main()
