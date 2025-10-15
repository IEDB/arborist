# Introduction

## What is Arborist?

Arborist is a software suite responsible for building and maintaining the biological taxonomy and other ontological trees for the Immune Epitope Database (IEDB). The name "Arborist" reflects its primary function: cultivating and managing tree-like hierarchical data structures that are essential for organizing the immunological data.

The project automates the complex process of integrating data from multiple standard biological databases to produce consistent and up-to-date ontologies for internal IEDB use.

## The Problem Arborist Solves

The IEDB contains millions of epitope records that need to be categorized by the organism they come from, the disease they are associated with, and the source proteins they are derived from. To make this data useful, it must be organized hierarchically. For example, a user should be able to find all epitopes associated with Type 1 Diabetes or all antigens associated with a _Mycobacterium tuberculosis_.

Arborist addresses this by creating four core trees:

1. **The Organism Tree**: A comprehensive taxonomy of organisms relevant to the IEDB, based on the NCBI Taxonomy but curated and pruned to meet the IEDB's specific needs.
2. **The Disease Tree**: A disease ontology built by combining terms from the Disease Ontology (DOID) and the IEDB's internal ontology (ONTIE).
3. **The Molecule (Protein) Tree**: A hierarchy that assigns IEDB source antigens and epitopes to specific genes and proteins from UniProt, organized under their respective species in the organism tree.

## How It Works

The main components of its workflow are:

* **Data Fetching**: Scripts pull data from an IEDB MySQL database, public resources (like UniProt and ontology sources), and curated internal files (which we call the Source of Truth or SoT).
* **Ontology Building**: The system uses `make` and the `ROBOT` tool to process, filter, and merge ontology files (`.owl`).
* **Tree Construction**: A series of Python scripts process the raw data to build parent-child relationships, assign metadata, and generate the final tree structures.
* **Web Interface**: The project uses **Nanobot** to serve the resulting data through a local web interface, allowing curators and developers to browse the trees and underlying data tables.