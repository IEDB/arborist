import polars as pl
from pathlib import Path


class DataFetcher:
  def __init__(self, build_path: Path = Path(__file__).parent.parent / 'build') -> None:
    self.build_path = build_path
  
  def get_all_peptides(self) -> pl.DataFrame:
    return pl.read_csv(self.build_path / 'iedb' / 'peptide.tsv', separator='\t')
  
  def get_all_sources(self) -> pl.DataFrame:
    return pl.read_csv(self.build_path / 'iedb' / 'peptide_source.tsv', separator='\t')

  def get_peptides_for_species(self, all_peptides: pl.DataFrame, all_taxa: list) -> pl.DataFrame:
    return all_peptides.filter(pl.col('Organism ID').is_in(all_taxa))

  def get_sources_for_species(self, all_sources: pl.DataFrame, accessions: list) -> pl.DataFrame:
    return all_sources.filter(pl.col('Source Accession').is_in(accessions))