#!/usr/bin/env python3
"""Build k=3 pepmatch indexes for the species that need them.

Arborist itself only ever uses k=5 (select_proteome.py `_preprocess_proteome`).
pepmatch-curation additionally serves 3-mer indexes for a fixed short list of
large, heavily queried species. Rather than teach the build about a second k,
this runs after a refresh and before a promote, over that list only.

Idempotent: an index newer than its proteome.fasta is left alone unless --force.
"""

import argparse
import sys
from pathlib import Path

from pepmatch import Preprocessor

# Fixed list, from what pepmatch-curation actually serves today.
THREE_MER_SPECIES = [
  9606,   # Homo sapiens
  10090,  # Mus musculus
  10116,  # Rattus norvegicus
  9612,   # Canis lupus
  9796,   # Equus caballus
  9823,   # Sus scrofa
  9913,   # Bos taurus
  9986,   # Oryctolagus cuniculus
]

K = 3


def index_path(species_path: Path) -> Path:
  return species_path / f'proteome_{K}mers.pepidx'


def needs_build(species_path: Path, force: bool = False) -> bool:
  """Stale or missing index, given a proteome that actually has content."""
  fasta = species_path / 'proteome.fasta'
  if not fasta.exists() or fasta.stat().st_size == 0:
    return False
  index = index_path(species_path)
  if force or not index.exists():
    return True
  return index.stat().st_mtime < fasta.stat().st_mtime


def build(species_path: Path) -> Path:
  Preprocessor(
    proteome=species_path / 'proteome.fasta',
    preprocessed_files_path=species_path,
  ).preprocess(k=K)
  return index_path(species_path)


def main(argv=None) -> int:
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-b', '--build_path', type=Path, required=True)
  parser.add_argument('-f', '--force', action='store_true',
                      help='rebuild even when the index is newer than the FASTA')
  parser.add_argument('-t', '--taxon', type=int, action='append',
                      help='override the species list (repeatable)')
  args = parser.parse_args(argv)

  species_ids = args.taxon or THREE_MER_SPECIES
  built, skipped, missing = 0, 0, []

  for taxon_id in species_ids:
    species_path = args.build_path / 'species' / str(taxon_id)
    if not (species_path / 'proteome.fasta').exists():
      missing.append(taxon_id)
      print(f'{taxon_id}: no proteome.fasta, skipping', file=sys.stderr)
      continue
    if not needs_build(species_path, args.force):
      skipped += 1
      print(f'{taxon_id}: {K}-mer index current', file=sys.stderr)
      continue
    print(f'{taxon_id}: building {K}-mer index', file=sys.stderr)
    print(f'{taxon_id}: wrote {build(species_path)}', file=sys.stderr)
    built += 1

  print(f'{K}-mers: {built} built, {skipped} current, {len(missing)} missing '
        f'of {len(species_ids)}', file=sys.stderr)
  # A missing proteome for one of these species is a real problem: they are the
  # ones prod queries most.
  return 1 if missing else 0


if __name__ == '__main__':
  sys.exit(main())
