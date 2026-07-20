"""Siloed, offline test for PeptideProcessor (epitope -> parent protein).

This is the IEDB-facing output. It runs fully offline: peptide search uses
pepmatch (a pip library, no bins/network), and with no allergens for the taxon
the allergen-blast branch is never taken. No MySQL, no bin/, no network.
"""

import polars as pl

TAXON = 999999
SPECIES_NAME = 'Testus specius'

PROTEOME_FASTA = (
  '>sp|P0001|GENEA_TEST Alpha protein OS=Testus specius OX=999999 GN=geneA PE=1 SV=1\n'
  'MSIINFEKLGGKLVAALKAGGRRRR\n'
  '>tr|Q0002|GENEB_TEST Beta protein OS=Testus specius OX=999999 GN=geneB PE=4 SV=1\n'
  'MKLVAALKAWWWWSIINFEKLDDDD\n'
)

PROTEOME_TSV = (
  'Database\tGene\tProtein ID\tEntry Name\tIsoform Count\tProtein Name\t'
  'Protein Existence Level\tSequence\n'
  'sp\tgeneA\tP0001\tGENEA_TEST\t1\tAlpha protein\t1\tMSIINFEKLGGKLVAALKAGGRRRR\n'
  'tr\tgeneB\tQ0002\tGENEB_TEST\t1\tBeta protein\t4\tMKLVAALKAWWWWSIINFEKLDDDD\n'
)


def _peptides():
  # Shape mirrors data_fetcher peptides already left-joined onto the source
  # assignments (what AssignmentHandler.process_species feeds PeptideProcessor).
  return pl.DataFrame({
    'Epitope ID': [900001],
    'Sequence': ['SIINFEKL'],
    'Starting Position': [1],
    'Ending Position': [8],
    'Source Accession': ['SRC1'],
    'Organism ID': [TAXON],
    'Organism Name': [SPECIES_NAME],
    'Source Alignment Score': [100.0],
    'Source Assigned Gene': [None],
    'Source Assigned Protein ID': ['P0001'],
    'Source Assigned Protein Name': ['Alpha protein'],
    'ARC Assignment': [None],
  }, schema_overrides={
    'Source Assigned Gene': pl.String, 'ARC Assignment': pl.String,
  })


def test_process_assigns_epitope_to_reviewed_protein_offline(silo):
  species_path = silo.species_path(TAXON)
  (species_path / 'proteome.fasta').write_text(PROTEOME_FASTA)
  (species_path / 'proteome.tsv').write_text(PROTEOME_TSV)

  processor = silo.assign.PeptideProcessor(
    taxon_id=TAXON, species_name=SPECIES_NAME, peptides=_peptides(),
    source_assignments=None, species_path=species_path,
  )
  processor.process()

  out = species_path / 'peptide-assignments.tsv'
  assert out.exists()
  df = pl.read_csv(out, separator='\t', infer_schema_length=0)
  row = df.filter(pl.col('Epitope Sequence') == 'SIINFEKL')
  assert row.height == 1
  # SIINFEKL occurs in both proteins; the reviewed (sp) P0001 must win.
  assert row.item(0, 'Assigned Protein ID') == 'P0001'
  assert row.item(0, 'Assigned Protein Entry Name') == 'GENEA_TEST'
  assert row.item(0, 'Assigned Protein Review Status') == 'true'
  # SIINFEKL is at residues 2-9 of MSIINFEKL... (1-indexed pepmatch positions).
  assert row.item(0, 'Assigned Protein Starting Position') == '2'
  assert row.item(0, 'Assigned Protein Ending Position') == '9'


def test_search_peptides_finds_all_occurrences_offline(silo):
  species_path = silo.species_path(TAXON)
  (species_path / 'proteome.fasta').write_text(PROTEOME_FASTA)

  processor = silo.assign.PeptideProcessor(
    taxon_id=TAXON, species_name=SPECIES_NAME, peptides=_peptides(),
    source_assignments=None, species_path=species_path,
  )
  processor.preprocess_proteome()
  processor.search_peptides()

  matches = pl.read_csv(species_path / 'peptide-matches.tsv', separator='\t')
  hits = matches.filter(pl.col('Query Sequence') == 'SIINFEKL')
  assert set(hits['Protein ID'].to_list()) == {'P0001', 'Q0002'}
