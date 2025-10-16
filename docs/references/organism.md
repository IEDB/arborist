# Organism Tree Reference

## Overview

The Organism Tree is the taxonomic backbone of the IEDB, providing a hierarchical classification of all organisms relevant to the IEDB. It combines data from the NCBI Taxonomy with IEDB-specific taxa and manual curation to create a pruned taxonomy for organizing epitope data.

## Structure and Components

### Three-Tier Architecture

The organism tree system consists of three related but distinct trees:

1. **NCBI Taxonomy**: The complete NCBI Taxonomy converted to LDTab format
2. **Organism Tree**: A curated subset focused on immunologically relevant organisms
3. **Subspecies Tree**: An extended version that includes all descendant taxa below species level

### Taxonomic Levels

Each taxon in the tree is assigned to one of three levels:

- **Upper**: High-level taxonomic ranks (kingdom, phylum, class, order, family, genus, etc.)
- **Species**: Taxa at the species rank, which serve as anchor points for proteome selection and data organization
- **Lower**: Subspecies, strains, and "no rank" taxa

### Data Sources

#### 1. NCBI Taxonomy
The primary source for taxonomic structure and scientific names. Downloaded from NCBI and converted to a local SQLite database using `ncbitaxon2ldtab.py`.

**Included Data:**
- Scientific names
- Taxonomic ranks
- Parent-child relationships
- Synonyms (common names, GenBank names, etc.)
- Alternative IDs for merged taxa
- PubMed citations

#### 2. IEDB Taxa
Custom taxa created by IEDB curators for specific immunological contexts, stored in the `iedb_taxa` table.

**Properties:**
- Unique IEDB IDs (format: `iedb-taxon:XXXXXXXX`)
- One or two parent taxa
- Custom labels and synonyms
- Rank and level assignments

#### 3. NCBI Include List
A curated list of NCBI taxa (`ncbi_include`) that should be explicitly included in the organism tree.

#### 4. Organism Core
The primary curation interface (`organism_core.tsv`) where curators define:
- Which taxa appear in the organism tree
- Custom labels (often including common names)
- Parent relationships
- Special nodes (marked with `use_other` flag)

## Build Process

### Step 1: Convert NCBI Taxonomy
```bash
make ncbitaxon
```

Downloads `taxdmp.zip` from NCBI and converts it to an LDTab SQLite database:
- Processes `names.dmp` for labels and synonyms
- Processes `nodes.dmp` for hierarchy and ranks
- Processes `merged.dmp` for alternative IDs
- Processes `citations.dmp` for literature references

### Step 2: Assign Species
```bash
python src/organism/assign_species.py
```

The `assign_species.py` script creates the foundation by:

1. **Loading organism_core.tsv**: Reads manually curated taxa
2. **Adding IEDB taxa**: Includes custom IEDB-defined organisms
3. **Processing NCBI includes**: Adds explicitly included NCBI taxa
4. **Identifying active taxa**: Adds organisms with epitope data
5. **Assigning species**: For each lower taxon, finds its species ancestor
6. **Assigning parents**: Determines parent taxa within the tree

**Output**: `organism-tree.tsv` with columns:
- `curie`: Unique identifier (NCBITaxon:XXXXX or iedb-taxon:XXXXX)
- `label`: Display name
- `label_source`: Origin of the label
- `rank`: Taxonomic rank
- `level`: Upper, species, or lower
- `epitope_count`: Number of associated epitopes (if any)
- `parent`: Primary parent taxon
- `parent_label`: Parent's display name
- `species`: Associated species (for lower taxa)
- `species_label`: Species display name

### Step 3: Build Organism Tree
```bash
python src/organism/build_organism_tree.py
```

Creates the `organism_tree` table in the LDTab database:

1. **Creates base structure**: Sets up annotation properties
2. **Inserts core triples**: Adds taxa from organism-tree.tsv as OWL classes
3. **Copies additional data**: Pulls relevant triples from NCBI Taxonomy
4. **Adds metadata**:
   - Taxon IDs (`ONTIE:0003615`)
   - NCBI browser URLs (`ONTIE:0003616`)
   - Taxonomic ranks (`ONTIE:0003617`)
   - Usage flags (`ONTIE:0003618`)
   - Label sources
   - IEDB synonyms
   - Epitope counts

### Step 4: Build Subspecies Tree
```bash
python src/organism/build_subspecies_tree.py
```

Extends the organism tree by:

1. **Copying organism_tree**: Starts with all curated taxa
2. **Adding descendants**: For each species, recursively adds all NCBI descendants
3. **Annotating species relationships**: Links lower taxa to their species using `ONTIE:0003619`

## Special Features

### "Other" Nodes

Taxa marked with `use_other = TRUE` in organism_core get special treatment:
- A sibling node `iedb-taxon:XXXXXX-other` is created
- Used for taxa that need a catch-all category for unclassified members
- Example: "Other viruses" under a viral family

### Species Assignment

Lower taxa are automatically assigned to their closest species ancestor:
- Enables proteome selection at the species level
- Groups strains and variants under their parent species
- Essential for protein tree construction

### Label Customization

Labels can be customized with common names:
```
Scientific name: Mycobacterium tuberculosis
Label with common name: Mycobacterium tuberculosis (tubercle bacillus)
Label source: NCBI Taxonomy scientific name (GenBank common name)
```

## Maintenance and Curation

### Validation and Checking

The `check_organism_core.py` script validates entries:

```bash
python src/organism/check_organism_core.py build/arborist/nanobot.db
```

**Checks performed:**
- Path to root exists
- Parent taxa exist in the tree
- Level assignments are appropriate
- Labels match NCBI (or document differences)
- Label sources are accurate

**Auto-update mode:**
```bash
python src/organism/check_organism_core.py -u build/arborist/nanobot.db
```
Automatically applies suggested fixes for info-level messages.

### Active Species

The `get_active_species.py` script identifies species with epitope data:

```bash
python src/organism/get_active_species.py
```

**Output** (`active-species.tsv`):
- Species Key: Unique identifier for proteome selection
- Species ID: NCBI taxon ID
- Species Label: Display name
- Active Taxa: Comma-separated list of descendants with epitopes
- Group: Biological category (virus, bacterium, vertebrate, etc.)
- Epitope Count: Total epitopes for this species and descendants

## File Formats

### organism_core.tsv

The source of truth for manual curation:

| Column | Description | Required |
|--------|-------------|----------|
| curie | Unique identifier | Yes |
| label | Display name | Yes |
| label_source | Origin of label | Yes |
| rank | Taxonomic rank | No |
| level | Upper/species/lower | Yes |
| parent | Parent taxon CURIE | Yes* |
| parent_label | Parent display name | No |
| use_other | Create "other" node | No |
| iedb_synonyms | Additional synonyms | No |

*Not required for root node (NCBITaxon:1)

### organism-tree.tsv

Intermediate file generated by `assign_species.py`:

Contains all taxa that will appear in the final organism tree, with:
- Complete parent-child relationships
- Species assignments for lower taxa
- Epitope counts for active taxa
- Source table tracking (organism_core, iedb_taxa, ncbi_include, count)

## Web Interface

### Tree Navigation

Browse the organism tree at `/arborist/organism_tree/NCBITaxon:1`:
- Hierarchical display with expand/collapse
- Search functionality with typeahead
- Links to NCBI Taxonomy Browser
- Quick access to related views

### Compare Trees

Compare different tree versions at `/arborist/ncbitaxon organism_tree/[TAXON_ID]`:
- Side-by-side comparison
- Highlights differences
- Useful for validation during curation

### Table View

View and edit organism_core at `/arborist/organism_core`:
- Sortable columns
- Filtering by any field
- Form view for adding/editing entries
- Validation messages

## Integration with Other Trees/Data

### Protein Tree

The organism tree serves as the foundation for the protein tree:
- Each species node gets a corresponding "species protein" node
- Proteins are organized under their species
- Gene and protein assignments depend on species identification

### Disease Tree

Disease associations use organism taxonomy:
- Link diseases to causative organisms
- Enable queries like "all diseases caused by Mycobacterium"

### Epitope Data

Epitopes are linked to organisms via:
- Direct organism ID on epitope records
- Source antigen organism
- Host organism for immune responses

## Command Reference

```bash
# Full organism tree build
make organism

# Individual steps
make ncbitaxon                           # Download and convert NCBI Taxonomy
make build/organism-tree.tsv             # Assign species and parents
make build/organism-tree.built           # Build organism_tree table
make build/subspecies-tree.built         # Build subspecies_tree table
make build/organism-tree.owl             # Export organism tree to OWL
make build/subspecies-tree.owl           # Export subspecies tree to OWL

# Validation and checking
python src/organism/check_organism_core.py build/arborist/nanobot.db
python src/organism/check_organism_core.py -u build/arborist/nanobot.db

# Active species
python src/organism/get_active_species.py \
  build/arborist/nanobot.db \
  build/iedb/peptide-count.tsv \
  build/arborist/active-species.tsv

# Render organism core tree
python src/organism/render_organism_core.py \
  build/arborist/organism_core.tsv \
  build/arborist/organism_core.html

# Sort organism core
python src/organism/sort_organism_core.py \
  build/arborist/organism_core.tsv
```