#!/usr/bin/env python3

import requests
import pandas as pd
from pathlib import Path


class DataFetcher:
  def __init__(self, build_path: Path = Path(__file__).parent.parent / 'build') -> None:
    self.build_path = build_path
  
  def get_all_peptides(self) -> pd.DataFrame:
    """Get all peptides from the written file."""
    return pd.read_csv(self.build_path / 'iedb' / 'peptide.tsv', sep='\t')
  
  def get_all_sources(self) -> pd.DataFrame:
    """Get all peptide sources from the written file."""
    return pd.read_csv(self.build_path / 'iedb' / 'peptide_source.tsv', sep='\t')

  def get_peptides_for_species(
    self, all_peptides: pd.DataFrame, all_taxa: list
  ) -> pd.DataFrame:
    """Get peptides from the written file only for a specific species.
    
    Args:
      all_peptides: list of all peptides from the backend. 
      all_taxa: list of all active children taxa for a species.
    """
    return all_peptides[all_peptides['Organism ID'].isin(all_taxa)]

  def get_sources_for_species(
    self, all_sources: pd.DataFrame, accessions: list
  ) -> pd.DataFrame:
    """Get peptide sources from the written file only for a specific species.
    
    Args:
      all_sources: list of all peptide sources from the backend.
      all_taxa: list of all active children taxa for a species.
    """
    return all_sources[all_sources['Source Accession'].isin(accessions)]

