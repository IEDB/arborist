# Protein Tree Reference

## Overview

The Protein Tree assigns IEDB source antigens and epitopes to their corresponding genes and proteins from UniProt. It builds upon the Organism Tree by selecting the best reference proteome for each active species and then mapping epitopes and source antigens to specific proteins using sequence similarity and structural analysis.

## Purpose and Goals

The protein tree solves several critical challenges:

1. **Standardization**: Maps diverse source antigen annotations to canonical UniProt proteins
2. **Hierarchy**: Organizes proteins under their species in a browsable tree structure
3. **Traceability**: Links epitopes to their source proteins with position information
4. **Gene-level organization**: Groups proteins by their genes for better biological context
5. **Fragment handling**: Identifies and represents protein fragments (chains, peptides, signal peptides)

## Architecture

### Two Types of Protein Trees

1. Protein Tree (without Gene Layer)
- Direct assignment: Species → Protein
2. Protein Tree (with Gene Layer)  
- Gene-mediated assignment: Species → Gene → Protein

### Tree Structure

```
BFO:0000040 (material entity)
└── PR:000000001 (protein)
    ├── iedb-protein:9606 (Homo sapiens protein)
    │   ├── iedb-protein:9606-HLA-A (Gene: HLA-A)
    │   │   ├── UP:P01892 (HLA class I histocompatibility antigen...)
    │   │   │   ├── UP:P01892-PRO_0000338644 (mature protein (25-365))
    │   │   │   └── UP:P01892-PRO_0000338645 (mature protein (25-308))
    │   │   └── UP:P30443 (HLA class I histocompatibility antigen...)
    │   ├── iedb-protein:9606-tcr-ig (T Cell Receptor chain)
    │   │   ├── iedb-protein:9606-tcr-ig-TRAV (Gene: TRAV)
    │   │   └── iedb-protein:9606-tcr-ig-TRBV (Gene: TRBV)
    │   └── iedb-protein:9606-other (Other Homo sapiens protein)
    └── iedb-protein:10090 (Mus musculus protein)
```

## Data Sources

### 1. IEDB Data

**Peptide Data** (`build/iedb/peptide.tsv`):
- Epitope sequences
- Source antigen accessions
- Position information (start/end)
- Organism IDs

**Source Data** (`build/iedb/peptide_source.tsv`):
- Source antigen accessions
- Database origins (UniProt, NCBI, etc.)
- Names and aliases
- Protein sequences
- Organism information

### 2. UniProt Proteomes

Selected for each active species via `ProteomeSelector`:

**Selection Criteria** (in priority order):
1. Reference and representative proteome
2. Representative proteome
3. Reference proteome
4. Other proteome types
5. Redundant proteomes

**Tiebreaking** (when multiple proteomes of same type exist):
1. BUSCO completeness score (if available and >20 candidates)
2. Epitope coverage (peptide matching across candidates)
3. Protein count (as a last resort)

**Downloaded Files**:
- `proteome.fasta`: Complete proteome sequences
- `gp_proteome.fasta`: Gene-prioritized reference proteome (if available)
- `proteome.tsv`: Parsed metadata

### 3. Curated Data

**Manual Assignments** (`build/arborist/manual-parents.tsv`):
- Source accessions that require manual protein assignment
- Overrides automatic BLAST/matching results

**Manual Synonyms** (`build/arborist/manual-synonyms.tsv`):
- Additional protein name synonyms
- Common names not in UniProt

**Allergen Data** (`build/arborist/allergens.json`):
- IUIS official allergen nomenclature
- Mapping to UniProt IDs
- Retrieved from WHO/IUIS Allergen Nomenclature database
- Manual overrides of parent protein results

## Build Process

### Step 1: Identify Active Species

```bash
make build/arborist/active-species.tsv
```

Uses `get_active_species.py` to create a list of species with epitope data, including:
- Species taxonomic information
- List of all active descendant taxa
- Biological group (virus, bacterium, vertebrate, etc.)
- Total epitope count

### Step 2: Select Proteomes

```bash
make proteome
# Or for a single species:
python src/protein/protein_tree/select_proteome.py -t 9606
```

**Process** (`ProteomeSelector` class):

1. **Query UniProt**: Fetch all proteomes for the species taxon ID
2. **Filter**: Remove excluded proteomes
3. **Categorize**: Group by proteome type
4. **Select**:
   - If single best type: use it
   - If multiple candidates: run tiebreaker
5. **Download**: Fetch FASTA and metadata
6. **Process**: Parse to TSV format
7. **Fetch Additional Data**:
   - Fragment information (chains, peptides, signal sequences)
   - Protein name synonyms
   - Gene-priority proteome (for reference proteomes)

**Special Cases**:
- **Orphan species**: No proteome available, fetch all proteins for taxon
- **Allergen species**: Always fetch as orphans for completeness
- **Replaced taxa**: Use mapping for recently merged NCBI taxa. This is an ongoing issue with the NCBI taxonomy changing over time, so there are manual overrides.

**Outputs** (per species in `build/species/{taxon_id}/`):
```
proteome-list.tsv      # All available proteomes
proteome.fasta         # Selected proteome FASTA file
proteome.tsv           # Parsed proteome metadata (like gene symbols)
gp_proteome.fasta      # Gene priority proteome (if available)
fragment-data.json     # Protein fragment information
synonym-data.json      # Additional protein synonyms
species-data.tsv       # Metadata on proteome selected
```

### Step 3: Assign Source Antigens

```bash
make protein
# Or for a single species (n = # of CPU cores to use):
python src/protein/protein_tree/assign.py -t 9606 -n 8
```

**Process** (`SourceProcessor` class):

#### 3a. Sequence Alignment

Uses BLAST (default) or MMseqs2 (for human, mouse, rat):

```bash
# BLAST approach
makeblastdb -in proteome.fasta -dbtype prot
blastp -query sources.fasta -db proteome.fasta -outfmt 10
```

```bash
# MMseqs2 approach (faster, more sensitive)
mmseqs easy-search sources.fasta proteome.fasta alignments.tsv tmp/
```

**Scoring**: `(Alignment Length × % Identity) / Query Length`

#### 3b. Top Protein Selection

For each source antigen:
1. Sort alignments by score and SwissProt review status
2. Select highest-scoring protein
3. Apply manual overrides from `manual-parents.tsv`

#### 3c. ARC Classification (Vertebrates Only)

For vertebrate species, run ARC (Antigen Receptor Classifier):

```python
SeqClassifier().classify_seqfile('sources.fasta')
```

**Identifies**:
- T Cell Receptors (TCR)
- B Cell Receptors / Immunoglobulins (BCR)
- Major Histocompatibility Complex (MHC) class I and II

**ARC Output**: Classification and chain type for immune receptors

#### 3d. Combine Results

Merge alignment data with:
- Proteome metadata (gene, protein name, review status)
- ARC classification (if applicable)
- Allergen mappings (if applicable)
- Manual synonyms

**Output**: `source-data.tsv` with columns:
```
Source ID, Source Accession, Database, Name, Aliases, 
Assigned Protein Synonyms, Taxon ID, Taxon Name, 
Species ID, Species Label, Proteome ID, Proteome Label,
Protein Strategy, Parent IRI, Parent Protein Database,
Parent Protein Accession, Parent Sequence Length, Sequence,
Parent Protein Gene
```

### Step 4: Assign Epitopes

Continues in `AssignmentHandler` → `PeptideProcessor`:

#### 4a. Preprocess Proteome with PEPMatch

```python
Preprocessor(proteome='proteome.fasta', k=5).sql_proteome()
```

#### 4b. Search Peptides with PEPMatch

```python
Matcher(query=peptides, k=5, max_mismatches=0).match()
```

**Two-stage matching**:

1. **With gene information**: 
   - Match peptide sequence AND source gene
   - Prefer proteins from the same gene as source antigen assigned

2. **Without gene information**:
   - Match peptide sequence to any protein
   - Prefer source antigen protein assignment if available

**Selection criteria** (in order):
- SwissProt reviewed status
- Exact source match (same protein as source)
- Gene priority (from gene-prioritized proteome)
- Protein existence level (1-5, lower is better)
- Protein ID (lexicographic sort for consistency)

#### 4c. Add Metadata

Enrich assignments with:
- Protein sequences and lengths
- Start and end positions
- Review status (SwissProt vs TrEMBL)
- Fragment information
- Synonyms (UniProt + manual)
- Species information

**Output**: `peptide-assignments.tsv` with columns:
```
Species Taxon ID, Species Name, Organism ID, Organism Name,
Source Accession, Source Alignment Score, Source Assigned Gene,
Source Assigned Protein ID, Source Assigned Protein Name,
ARC Assignment, Epitope ID, Epitope Sequence,
Source Starting Position, Source Ending Position,
Assigned Protein ID, Assigned Protein Name, 
Assigned Protein Entry Name, Assigned Protein Review Status,
Assigned Protein Starting Position, Assigned Protein Ending Position,
Assigned Protein Sequence, Assigned Protein Length,
Assigned Protein Fragments, Assigned Protein Synonyms
```

### Step 5: Combine Data

```bash
make build/arborist/all-peptide-assignments.tsv
```

Combines assignments from all species into unified tables:
- `all-peptide-assignments.tsv`: All epitope-to-protein mappings
- `all-source-data.tsv`: All source antigen data
- `all-species-data.tsv`: All proteome selection metadata

### Step 6: Build Tree Structure

```bash
make build/arborist/protein-tree.owl
make build/arborist/protein-tree-with-gene.owl
```

**Process** (`build_tree()` function):

#### 6a. Create Base Structure

Copy organism tree and transform:
- Replace `NCBITaxon:` with `iedb-protein:`
- Replace `iedb-taxon:` with `iedb-protein:`
- Add top-level `PR:000000001` (protein)
- Re-parent to protein instead of root
- Remove lower-level taxa (keep upper and species only)

#### 6b. Add Protein Nodes

**Three categories**:

1. **Normal Proteins**:
   - Group by gene (if gene layer enabled)
   - Create `UP:{accession}` nodes
   - Parent: `iedb-protein:{taxon_id}` or gene node

2. **Antigen Receptor Proteins** (from ARC):
   - Create receptor category nodes (TCR, BCR, MHC-I, MHC-II)
   - Group by gene within receptor category
   - Parent: `iedb-protein:{taxon_id}-{receptor_type}`

3. **Other Proteins** (no assignment):
   - Create "other" category nodes
   - Parent: `iedb-protein:{taxon_id}-other`

#### 6c. Add Protein Fragments

For proteins with multiple structural features:
- Chains (mature proteins)
- Peptides
- Propeptides  
- Signal peptides
- Transit peptides

Create child nodes: `UP:{accession}-{fragment_id}`

**Fragment criteria**:
- At least 2 fragments required
- Skip fragments spanning entire protein
- Include position information (start-end)

#### 6d. Add Metadata

For each protein node:
- **ONTIE:0003673**: Canonical status (SwissProt reviewed)
- **ONTIE:0003622**: Synonyms (exact synonym)
- **ONTIE:0003623**: Accession
- **ONTIE:0003624**: UniProt URL
- **ONTIE:0003625**: Source database
- **ONTIE:0003627**: Fragment start position
- **ONTIE:0003628**: Fragment end position
- **ONTIE:0003620**: Fragment label
- **ONTIE:0003674**: Gene symbol

**Outputs**:
- `protein_tree_old` table: Without gene layer
- `protein_tree_new` table: With gene layer

## ImmunomeBrowser

### Step 7: Create Immunome Files

```bash
make build/arborist/epitope-mappings.tsv
# Runs src/protein/protein_tree/immunomebrowser.py
```

Generates files for the IEDB Immunome Browser:

#### 7a. Source Parents (`source-parents.tsv`)

Links source antigens to their parent proteins:

```
Source ID, Accession, Database, Name, Aliases, Synonyms,
Taxon ID, Taxon Name, Species ID, Species Label,
Proteome ID, Proteome Label, Protein Strategy,
Parent IRI, Parent Protein Database, Parent Protein Accession,
Parent Sequence Length, Sequence, Parent Protein Gene
```

#### 7b. Parent Proteins (`parent-proteins.tsv`)

Unique list of all parent proteins:

```
Accession, Database, Name, Title, 
Proteome ID, Proteome Label, Sequence
```

#### 7c. Epitope Mappings (`epitope-mappings.tsv`)

Detailed epitope-to-protein mappings with alignment information:

```
epitope_id, epitope_seq, epitope_start, epitope_end,
source_accession, parent_accession, parent_seq,
parent_start, parent_end, identity_alignment,
similarity_alignment, gaps_source_alignment,
gaps_parent_alignment, all_gaps, source_alignment,
parent_alignment, parent_alignment_modified
```

**Three mapping types**:

1. **Exact matches**: Epitope positions known from PEPMatch
2. **Linear peptides**: BLAST alignment of unpositioned linear epitopes
3. **Discontinuous epitopes**: Residue-by-residue matching

**Alignment strings**:
- `source_alignment`: Epitope sequence with gaps
- `parent_alignment`: Protein sequence with gaps
- `parent_alignment_modified`: Mismatches in lowercase

## File Formats

### Active Species (`active-species.tsv`)

| Column | Description |
|--------|-------------|
| Species Key | Unique key for file naming |
| Species ID | NCBI Taxon ID |
| Species Label | Scientific name |
| Active Taxa | Comma-separated descendant taxa with epitopes |
| Group | Biological category |
| Epitope Count | Total epitopes for species |

### Proteome Metadata (`proteome.tsv`)

| Column | Description |
|--------|-------------|
| Database | sp (SwissProt) or tr (TrEMBL) |
| Gene | Gene symbol |
| Protein ID | UniProt accession |
| Entry Name | UniProt entry name |
| Isoform Count | Isoform number (or "1") |
| Protein Name | Full protein name |
| Protein Existence Level | 1-5 evidence level |
| Gene Priority | 1 if in gene-priority proteome |
| Sequence | Amino acid sequence |

### Peptide Assignments (`peptide-assignments.tsv`)

Primary output linking epitopes to proteins. See Step 4d above for full column list.

**Key relationships**:
- `Epitope ID` → `Epitope Sequence`
- `Source Accession` → `Source Assigned Protein ID`
- `Epitope Sequence` + `Assigned Protein ID` → `Assigned Protein Starting Position`

## Parameters

### Command-Line Arguments

**select_proteome.py**:
```bash
-t, --taxon_id      Species taxon ID to process
-b, --build_path    Path for build files (default: build/)
```

**assign.py**:
```bash
-t, --taxon_id      Species taxon ID to process
-b, --build_path    Path for build files (default: build/)
-n, --num_threads   Number of threads for BLAST/MMseqs2 (default: 1)
```

**immunomebrowser.py**:
```bash
-n, --num_threads   Number of threads for BLAST (default: 1)
```

### Algorithm Parameters

**BLAST**:
- E-value: 1.0 (permissive to catch distant matches)
- Output format: CSV with custom fields
- Threads: Configurable via `-n`

**MMseqs2**:
- Sensitivity: 7.0 (high sensitivity)
- Used for human, mouse, rat (taxon IDs: 9606, 10090, 10116)

**PEPMatch**:
- k-mer size: 5
- Max mismatches: 0 (exact matches only)
- Best match: False (return all matches)

**ARC**:
- Uses default HMM profiles for TCR, BCR, MHC
- BLAST and HMMER for classification
- Outputs class, chain type, and calculated MHC allele

## Tools and Dependencies

### Core Tools

**External Binaries** (auto-installed to `bin/`):
- `blastp` / `makeblastdb`
- `mmseqs`
- `hmmscan`

### Python Packages

**From `requirements.txt`**:
- `biopython>=1.78`
- `polars>=1.0`
- `pepmatch>=0.9.4`
- `ARC`

## Command Reference

```bash
# Full protein tree build
make protein

# Individual steps
make build/arborist/active-species.tsv           # Identify active species
make proteome                                    # Select proteomes (all species)
make build/arborist/all-peptide-assignments.tsv  # Assign all
make build/arborist/protein-tree.owl             # Build tree structure
make build/arborist/protein-tree-with-gene.owl   # Build tree structure with gene layer
make build/arborist/epitope-mappings.tsv         # Generate ImmunomeBrowser files

# Single species operations
python src/protein/protein_tree/select_proteome.py -t 9606 -b build/
python src/protein/protein_tree/assign.py -t 9606 -n 8 -b build/

# Generate new mappings (this is comparing previous build to newest build)
python src/util/generate_new_mappings.py \
  build/proteins/latest/epitope-mappings.tsv \
  build/proteins/previous/epitope-mappings.tsv \
  build/proteins/diffs/new-mappings.tsv
```