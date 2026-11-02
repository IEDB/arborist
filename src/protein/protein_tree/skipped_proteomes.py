"""Per-species record of UniProt proteome fetches that failed after all retries.

A proteome that cannot be fetched must never crash the build: select_proteome
drops that candidate (or, when nothing usable remains, skips the taxon) and
records it here instead. The record lives next to the species' other outputs
(build/species/<id>/fetch-skipped.tsv) so that

  * the next build can see it and retry (assign.check_for_proteome and
    select_proteome.proteome_complete treat a species with a record as
    unfinished, so a skipped proteome is never cached as a final answer), and
  * the end-of-run summary and the weekly digest email can list every skip.
"""

import csv
from pathlib import Path

SKIP_FILE = 'fetch-skipped.tsv'
FIELDS = ['Species ID', 'Species Name', 'Proteome ID', 'Outcome', 'Error']

# Other candidates were usable; the best remaining one was selected.
OUTCOME_CANDIDATE_DROPPED = 'candidate dropped'
# No candidate was usable; the taxon gets an empty proteome and is skipped.
OUTCOME_TAXON_SKIPPED = 'taxon skipped'


def _clean(value) -> str:
  """One-line, tab-free text so a record is always a single TSV row."""
  return ' '.join(str(value).split())


def has_skips(species_path: Path) -> bool:
  return (Path(species_path) / SKIP_FILE).exists()


def clear_skips(species_path: Path) -> None:
  path = Path(species_path) / SKIP_FILE
  if path.exists():
    path.unlink()


def record_skip(species_path: Path, species_id, species_name: str, proteome_id: str,
                outcome: str, error) -> None:
  path = Path(species_path) / SKIP_FILE
  is_new = not path.exists()
  with open(path, 'a', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    if is_new:
      writer.writerow(FIELDS)
    writer.writerow([
      _clean(species_id), _clean(species_name), _clean(proteome_id), outcome, _clean(error),
    ])


def load_skips(build_path: Path) -> list:
  """Every recorded skip under build/species/, ordered by species then record."""
  rows = []
  species_dir = Path(build_path) / 'species'
  if not species_dir.is_dir():
    return rows
  files = sorted(species_dir.glob(f'*/{SKIP_FILE}'), key=lambda p: (not p.parent.name.isdigit(), p.parent.name.zfill(12)))
  for path in files:
    with open(path, newline='') as f:
      rows.extend(csv.DictReader(f, delimiter='\t'))
  return rows


def format_skips(skips: list) -> str:
  """Plain-text list for the end-of-run log; empty string when nothing was skipped."""
  if not skips:
    return ''
  lines = [f'{len(skips)} UniProt proteome fetch(es) failed and were skipped; they are retried on the next build:']
  for row in skips:
    lines.append(
      f"  taxon {row['Species ID']} ({row['Species Name']}): proteome {row['Proteome ID'] or '-'} "
      f"[{row['Outcome']}] {row['Error']}"
    )
  return '\n'.join(lines)
