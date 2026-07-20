"""Siloed, offline tests for SourceProcessor (source antigen -> parent protein).

The only external dependency in this path is the alignment step
(blastp/mmseqs2), which writes alignments.csv / alignments.tsv. We freeze a
tiny alignment table and exercise the pure-polars logic directly, so these run
with no bins, no MySQL, and no network -- on asahidake included.
"""

import polars as pl

TAXON = 108  # bacterium -> blastp path (not in the mmseqs2 list), no ARC step

PROTEOME_TSV = (
  'Database\tGene\tProtein ID\tEntry Name\tIsoform Count\tProtein Name\t'
  'Protein Existence Level\tSequence\n'
  'sp\tgeneA\tP0001\tGENEA_TEST\t1\tAlpha protein\t1\tMSIINFEKLGGKLVAALKAGGRRRR\n'
  'tr\tgeneB\tQ0002\tGENEB_TEST\t1\tBeta protein\t4\tMKLVAALKAWWWWSIINFEKLDDDD\n'
)

# blastp -outfmt 10 output: comma-separated, no header. SRC1 hits both proteins;
# the reviewed P0001 wins on Score. SRC2 has no hit -> must survive as unassigned.
ALIGNMENTS_CSV = (
  'SRC1,sp|P0001|GENEA_TEST,100.0,25,0,0,1,25,1,25,1e-30,60.0\n'
  'SRC1,tr|Q0002|GENEB_TEST,60.0,20,5,0,1,20,1,20,1e-03,30.0\n'
)


def _sources():
  return pl.DataFrame({
    'Source Accession': ['SRC1', 'SRC2'],
    'Sequence': ['MSIINFEKLGGKLVAALKAGGRRRR', 'MKKKKKKKKKK'],
    'Length': [25, 11],
  })


def _processor(silo, alignments=ALIGNMENTS_CSV):
  species_path = silo.species_path(TAXON)
  (species_path / 'proteome.tsv').write_text(PROTEOME_TSV)
  if alignments is not None:
    (species_path / 'alignments.csv').write_text(alignments)
  return silo.assign.SourceProcessor(
    taxon_id=TAXON, group='bacterium', sources=_sources(),
    num_threads=1, species_path=species_path,
  )


def test_select_top_proteins_ranks_by_score_then_review(silo):
  top = _processor(silo).select_top_proteins()
  src1 = top.filter(pl.col('Query') == 'SRC1')
  assert src1.height == 1
  assert src1.item(0, 'Subject') == 'P0001'  # higher Score, reviewed


def test_select_top_proteins_keeps_sources_without_a_hit(silo):
  top = _processor(silo).select_top_proteins()
  assert set(top['Query'].to_list()) == {'SRC1', 'SRC2'}
  assert top.filter(pl.col('Query') == 'SRC2').item(0, 'Subject') is None


def test_process_offline_assigns_reviewed_protein(silo, monkeypatch):
  proc = _processor(silo)
  # blast is external: no-op it, leaving the frozen alignments.csv in place.
  monkeypatch.setattr(proc, 'blast_sources', lambda: None)
  result = proc.process()
  src1 = result.filter(pl.col('Source Accession') == 'SRC1')
  assert src1.item(0, 'Source Assigned Protein ID') == 'P0001'
  assert src1.item(0, 'Source Assigned Gene') == 'geneA'
  assert src1.item(0, 'ARC Assignment') is None  # bacterium: no ARC classification


def test_combine_arc_data_builds_class_string(silo):
  proc = _processor(silo)
  top_protein_data = pl.DataFrame({'Source Accession': ['SRC1'], 'keep': [1]})
  arc_df = pl.DataFrame({
    'id': ['SRC1'], 'class': ['abc'], 'chain_type': ['TRA'], 'calc_mhc_allele': ['HLA'],
  })
  out = proc.combine_arc_data(top_protein_data, arc_df)
  assert out.filter(pl.col('Source Accession') == 'SRC1').item(0, 'ARC Assignment') == 'abc_TRA_HLA'
