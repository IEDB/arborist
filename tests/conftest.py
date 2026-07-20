import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

# Make the editable `protein_tree` package importable without a full install,
# so the offline tests run on any box (asahidake included, no bin/, no MySQL).
SRC_PROTEIN = Path(__file__).parent.parent / 'src' / 'protein'
if str(SRC_PROTEIN) not in sys.path:
  sys.path.insert(0, str(SRC_PROTEIN))


@pytest.fixture
def silo(tmp_path, monkeypatch):
  """A hermetic, offline build tree + wired `assign` module globals.

  Gives each test its own build/ under tmp_path with header-only arborist
  inputs, and points assign.build_path at it. No network, no MySQL, no bins.
  """
  from protein_tree import assign

  build = tmp_path / 'build'
  arborist = build / 'arborist'
  arborist.mkdir(parents=True)
  (arborist / 'manual-parents.tsv').write_text(
    'Organism ID\tOrganism Name\tAccession\tName\t'
    'Parent Database\tParent Accession\tParent Name\tGene\n'
  )
  (arborist / 'manual-synonyms.tsv').write_text('Accession\tSynonyms\n')
  # A numeric placeholder row so SpeciesID infers as int (matches production);
  # SpeciesID 0 never matches a real taxon, so create_allergen_fasta skips it.
  (arborist / 'allergens.tsv').write_text('SpeciesID\tName\tSequence\n0\t\t\n')

  monkeypatch.setattr(assign, 'build_path', build, raising=False)

  def species_path(taxon_id):
    path = build / 'species' / str(taxon_id)
    path.mkdir(parents=True, exist_ok=True)
    return path

  return SimpleNamespace(
    assign=assign, build=build, arborist=arborist, species_path=species_path,
  )
