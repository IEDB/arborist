"""Siloed, offline tests for AssignmentHandler.generate_proteome_tsv.

Parses a UniProt-style proteome.fasta into proteome.tsv. Pure Biopython + regex
-- no bins, no network. Also covers the defensive early-returns (empty and
non-FASTA payloads) that must never crash a run.
"""

import polars as pl

TAXON = 424242

PROTEOME_FASTA = (
  '>sp|P0001|GENEA_TEST Alpha protein OS=Testus specius OX=424242 GN=geneA PE=1 SV=1\n'
  'MSIINFEKLGGKLVAALKAGG\n'
  '>tr|Q0002|GENEB_TEST Beta protein OS=Testus specius OX=424242 GN=geneB PE=4 SV=1\n'
  'MKLVAALKAWWWWSIINFEKL\n'
)


def _handler(silo, taxon=TAXON):
  empty = pl.DataFrame()
  return silo.assign.AssignmentHandler(
    taxon_id=taxon, species_name='Testus specius', group='bacterium',
    peptides=empty, sources=empty, num_threads=1,
  )


def test_generate_proteome_tsv_parses_metadata(silo):
  species_path = silo.species_path(TAXON)
  (species_path / 'proteome.fasta').write_text(PROTEOME_FASTA)

  _handler(silo).generate_proteome_tsv()

  tsv = species_path / 'proteome.tsv'
  assert tsv.exists()
  df = pl.read_csv(tsv, separator='\t', infer_schema_length=0)

  p = df.filter(pl.col('Protein ID') == 'P0001').row(0, named=True)
  assert p['Database'] == 'sp'
  assert p['Gene'] == 'geneA'
  assert p['Entry Name'] == 'GENEA_TEST'
  assert p['Protein Name'] == 'Alpha protein'
  assert p['Protein Existence Level'] == '1'
  assert p['Sequence'] == 'MSIINFEKLGGKLVAALKAGG'

  q = df.filter(pl.col('Protein ID') == 'Q0002').row(0, named=True)
  assert q['Database'] == 'tr'
  assert q['Protein Existence Level'] == '4'


def test_generate_proteome_tsv_skips_empty_proteome(silo):
  species_path = silo.species_path(TAXON)
  (species_path / 'proteome.fasta').touch()  # zero bytes

  _handler(silo).generate_proteome_tsv()  # must not raise
  assert not (species_path / 'proteome.tsv').exists()


def test_generate_proteome_tsv_survives_non_fasta_payload(silo):
  species_path = silo.species_path(TAXON)
  # e.g. a UniProt stream-error body that slipped through instead of FASTA.
  (species_path / 'proteome.fasta').write_text(
    'Error encountered when streaming data. Please try again later.\n'
  )

  _handler(silo).generate_proteome_tsv()  # must not raise
  assert not (species_path / 'proteome.tsv').exists()
