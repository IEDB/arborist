import pandas as pd
from pathlib import Path


class DataFetcher:
  def __init__(self, build_path: Path = Path(__file__).parent.parent / 'build') -> None:
    self.build_path = build_path
  
  def get_all_peptides(self) -> pd.DataFrame:
    return pd.read_csv(self.build_path / 'iedb' / 'peptide.tsv', sep='\t')
  
  def get_all_sources(self) -> pd.DataFrame:
    return pd.read_csv(self.build_path / 'iedb' / 'peptide_source.tsv', sep='\t')

  def get_peptides_for_species(self, all_peptides: pd.DataFrame, all_taxa: list) -> pd.DataFrame:
    """Get peptides for a species using a list of children taxa."""
    return all_peptides[all_peptides['Organism ID'].isin(all_taxa)]

  def get_sources_for_species(self, all_sources: pd.DataFrame, accessions: list) -> pd.DataFrame:
    """Get sources for a species using a list of accessions."""
    return all_sources[all_sources['Source Accession'].isin(accessions)]