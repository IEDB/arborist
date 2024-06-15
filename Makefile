# Arborist Makefile
# James A. Overton <james@overton.ca>
#
# Arborist builds "trees" that drive various IEDB finders.
# It merges IEDB data with community ontologies and other databases.
# The build follows these steps:
#
# 1. Fetch IEDB data: requires MySQL connection
# 2. Build NCBI Taxonomy
# 3. Build Organism Tree and determine active species
# 4. Select Proteomes for each active species
# 5. Assign Proteins
# 6. Build Protein Tree
# 7. TODO Build Molecule Tree
# 8. TODO Build Assay Tree
# 9. TODO Build Disease Tree
#
# TODO: test suite
# TODO: geolocation tree, MHC tree
# TODO: merged SoT tree
# TODO: test data, symlink?


### Configuration
#
# These are standard options to make Make sane:
# <http://clarkgrubb.com/makefile-style-guide#toc2>

MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help
.DELETE_ON_ERROR:
.PRECIOUS:
.SUFFIXES:

export PATH := $(shell pwd)/bin:$(PATH)


### Main Tasks
#
# The main tasks for running Arborist.

.PHONY: help
help:
	@echo "Arborist: build trees for the IEDB"
	@echo ""
	@echo "TASKS"
	@echo "  deps        install dependencies"
	@echo "  iedb        load IEDB data"
	@echo "  ncbitaxon   build the NCBI Taxonomy"
	@echo "  organism    build the organism and subspecies trees"
	@echo "  proteome    select proteomes"
	@echo "  protein     build the protein tree"
	@echo "  molecule    build the molecule tree"
	@echo "  leidos      copy files for Leidos"
	@echo "  all         build all trees"
	@echo "  serve       run the web interface on localhost:3000"
	@echo "  clean       remove all build files"
	@echo "  clobber     remove all generated files"
	@echo "  help        print this message"

# Dependencies are added to this list below.
.PHONY: deps
deps:

.PHONY: all
all: deps iedb ncbitaxon organism protein molecule disease leidos

.PHONY: leidos
leidos: build/organisms/latest/ build/proteins/latest/

.PHONY: serve
serve: src/util/serve.py
	python3 $<

.PHONY: clean
clean:
	rm -rf build

.PHONY: clobber
clobber:
	chmod +w -R cache/
	rm -rf bin/ build/ cache/ current/

bin/ build/ cache/ current/:
	mkdir -p $@


### Install Dependencies
#
# For each software dependency
# we use Make's `ifeq` conditional and `command -v <name>`
# to check if the dependency is already on the PATH
# (including `bin/`).
# If not, then we define a Make task to install it to `bin/`,
# and add that dependency to `deps`.

# Require SQLite
ifeq ($(shell command -v sqlite3),)
$(error 'Please install SQLite 3')
endif

# Require MySQL or MariaDB
ifeq ($(shell command -v mysql),)
$(error "Please install 'mariadb' from MariaDB")
endif

# Require Python
ifeq ($(shell command -v python3),)
$(error 'Please install Python 3, so we can run various scripts')
endif

# Require Java
ifeq ($(shell command -v java),)
$(error 'Please install Java, so we can run ROBOT and LDTab')
endif

# Install ROBOT if not already present
ifeq ($(shell command -v robot),)
bin/robot.jar: | bin/
	curl -L -o $@ 'https://github.com/ontodev/robot/releases/download/v1.9.4/robot.jar'
bin/robot: bin/robot.jar
	curl -L -o $@ 'https://raw.githubusercontent.com/ontodev/robot/master/bin/robot'
	chmod +x $@
deps: bin/robot
endif

# Install LDTab if not already present
ifeq ($(shell command -v ldtab),)
bin/ldtab.jar: | bin/
	curl -L -o $@ 'https://github.com/ontodev/ldtab.clj/releases/download/v2023-12-21/ldtab.jar'
bin/ldtab: bin/ldtab.jar
	echo '#!/bin/sh' > $@
	echo 'java -jar "$$(dirname $$0)/ldtab.jar" "$$@"' >> $@
	chmod +x $@
deps: bin/ldtab
endif

# Install Nanobot if not already present
ifeq ($(shell command -v nanobot),)
bin/nanobot: | bin/
	curl -L -k -o $@ 'https://github.com/ontodev/nanobot.rs/releases/download/v2023-10-26/nanobot-x86_64-unknown-linux-musl'
	chmod +x $@
deps: bin/nanobot
endif

# Install valve-export script if not already present
ifeq ($(shell command -v valve-export),)
bin/valve-export: | bin/
	curl -L -o $@ 'https://github.com/ontodev/valve.rs/raw/main/scripts/export.py'
	chmod +x $@
deps: bin/valve-export
endif

# Install QSV if not already present
ifeq ($(shell command -v qsv),)
QSV_VERSION := 0.118.0
bin/qsv: | bin/ build/
	curl -L -k -o build/qsv.zip 'https://github.com/jqnatividad/qsv/releases/download/$(QSV_VERSION)/qsv-$(QSV_VERSION)-x86_64-unknown-linux-musl.zip'
	cd build && unzip qsv.zip qsv
	mv build/qsv $@
deps: bin/qsv
endif

# Install BLAST if not already present
ifeq ($(shell command -v blastp),)
BLAST_VERSION := 2.15.0
build/ncbi-blast.tar.gz: | build/
	curl -L -o $@ 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(BLAST_VERSION)/ncbi-blast-$(BLAST_VERSION)+-x64-linux.tar.gz'
bin/blastp bin/makeblastdb: build/ncbi-blast.tar.gz | bin/
	cd build/ && tar -zxvf $(notdir $<) ncbi-blast-$(BLAST_VERSION)+/$@
	mv build/ncbi-blast-$(BLAST_VERSION)+/$@ $@
deps: bin/blastp bin/makeblastdb
endif

# Install HMMER if not already present
ifeq ($(shell command -v hmmscan),)
HMMER_VERSION := 3.4
build/hmmer-$(HMMER_VERSION).tar.gz: | build/
	curl -L -o $@ 'http://eddylab.org/software/hmmer/hmmer-$(HMMER_VERSION).tar.gz'
build/hmmer-$(HMMER_VERSION): build/hmmer-$(HMMER_VERSION).tar.gz
	cd build/ && tar xvf $(notdir $<)
bin/hmmscan: build/hmmer-$(HMMER_VERSION) | bin/
	cd $< && ./configure --prefix $(shell pwd)/$< && make install
	cp $</$@ $@
deps: bin/hmmscan
endif

# Install MMseqs2 if not already present
ifeq ($(shell command -v mmseqs),)
MMSEQS_VERSION := 15-6f452
build/mmseqs-linux-sse41.tar.gz: | build/
	curl -L -o $@ 'https://github.com/soedinglab/MMseqs2/releases/download/$(MMSEQS_VERSION)/mmseqs-linux-sse41.tar.gz'
build/mmseqs: build/mmseqs-linux-sse41.tar.gz
	cd build/ && tar xvf $(notdir $<) && mv mmseqs mmseqs-linux-sse41
bin/mmseqs: build/mmseqs | bin/
	cp build/mmseqs-linux-sse41/bin/mmseqs $@
	rm -rf build/mmseqs-linux-sse41
deps: bin/mmseqs
endif

### 1. Fetch IEDB Data
#
# Copy some tables from the IEDB,
# validate some of them,
# and load them into a SQLite database.

build/iedb/:
	mkdir -p $@

build/iedb/nanobot.toml: src/iedb/nanobot.toml | build/iedb/
	cp $< $@

# Update the cached IEDB data if current/iedb/ is not yet set.
# Requires MySQL/MariaDB connection details.
current/iedb/%:
	src/iedb/update-cache "$$IEDB_MYSQL_DATABASE"

# Unzip an IEDB table from the cache into the build directory,
build/iedb/%.tsv: current/iedb/%.tsv.gz | build/iedb/
	zcat $< > $@

# Load IEDB tables into SQLite using Nanobot.
# For some tables we just create a header row, with no data rows,
# then remove those files after Nanobot has initialized.
# TODO: This task will be much simpler with planned VALVE features.
build/iedb/nanobot.db: build/iedb/nanobot.toml
build/iedb/nanobot.db: build/iedb/ncbi_include.tsv
build/iedb/nanobot.db: build/iedb/iedb_taxa.tsv
build/iedb/nanobot.db: current/iedb/source.tsv.gz
build/iedb/nanobot.db: current/iedb/object.tsv.gz
build/iedb/nanobot.db: current/iedb/epitope.tsv.gz
	rm -f $@
	zcat current/iedb/source.tsv.gz | head -n1 > build/iedb/source.tsv || exit 0
	zcat current/iedb/object.tsv.gz | head -n1 > build/iedb/object.tsv || exit 0
	zcat current/iedb/epitope.tsv.gz | head -n1 > build/iedb/epitope.tsv || exit 0
	# create peptide.tsv with just headers
	qsv search --select table 'peptide\b' src/iedb/column.tsv \
	| qsv select column \
	| qsv behead \
	| tr '\n' '\t' \
	| sed 's/	$$//' \
	> build/iedb/peptide.tsv
	# create peptide_source.tsv with just headers
	qsv search --select table 'peptide_source\b' src/iedb/column.tsv \
	| qsv select column \
	| qsv behead \
	| tr '\n' '\t' \
	| sed 's/	$$//' \
	> build/iedb/peptide_source.tsv
	cd build/iedb/ && nanobot init
	rm -f build/iedb/source.tsv
	rm -f build/iedb/object.tsv
	rm -f build/iedb/epitope.tsv
	rm -f build/iedb/peptide.tsv
	rm -f build/iedb/peptide_source.tsv

# Load a table in SQLite, without VALVE validation.
build/iedb/%.built: build/iedb/%.tsv | build/iedb/nanobot.db
	sqlite3 $| "DELETE FROM '$*'"
	src/util/tsv2sqlite $| $<
	touch $@

# Extract tables of peptides and sources from SQLite tables.
build/iedb/peptide.tsv: src/iedb/peptide.sql build/iedb/epitope.built build/iedb/object.built | build/iedb/nanobot.db
	src/util/sqlite2tsv $| $< $@

build/iedb/peptide_source.tsv: src/iedb/peptide_source.sql build/iedb/source.built | build/iedb/nanobot.db
	src/util/sqlite2tsv $| $< $@

build/iedb/structure.tsv: build/iedb/peptide.tsv
	zcat current/iedb/structure.tsv.gz | cut -f1,2,13,15 > $@
	src/util/map_structure_ids.py $@ $<

.PHONY: iedb
iedb: build/iedb/peptide.built build/iedb/peptide_source.built build/iedb/structure.tsv

### 2. Build NCBI Taxonomy
#
# Set up a Nanobot instance.
# Fetch the NCBI Taxonomy `taxdmp.zip` file
# and use it to populate a table.

build/arborist/:
	mkdir -p $@

build/arborist/nanobot.toml: src/arborist/nanobot.toml | build/arborist/
	cp $< $@

# Initialize a Nanobot database.
# Create an empty organism-tree.tsv.
build/arborist/nanobot.db: build/arborist/nanobot.toml src/arborist/*.tsv
	rm -f $@ $(dir $@)*.built
	cd build/arborist/ && nanobot init

cache/ncbitaxon/:
	mkdir $@

TAXDMP_VERSION := $(shell date +"%Y-%m-01")

# Fetch the taxdmp.zip for this month.
cache/ncbitaxon/taxdmp_$(TAXDMP_VERSION).zip: | cache/ncbitaxon/ current/
	curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_$(TAXDMP_VERSION).zip

current/taxdmp.zip: cache/ncbitaxon/taxdmp_$(TAXDMP_VERSION).zip
	rm -f $@
	cd current/ && ln -s ../$< taxdmp.zip

build/arborist/ncbitaxon.built: src/organism/ncbitaxon2ldtab.py current/taxdmp.zip | build/arborist/nanobot.db
	sqlite3 $| "DROP TABLE IF EXISTS ncbitaxon"
	python3 $^ $|
	sqlite3 $| "CREATE INDEX idx_ncbitaxon_subject ON ncbitaxon(subject)"
	sqlite3 $| "CREATE INDEX idx_ncbitaxon_predicate ON ncbitaxon(predicate)"
	sqlite3 $| "CREATE INDEX idx_ncbitaxon_object ON ncbitaxon(object)"
	sqlite3 $| "ANALYZE ncbitaxon"
	touch $@

.PHONY: ncbitaxon
ncbitaxon: build/arborist/ncbitaxon.built


### 3. Build Organism Tree
#
# Build the organism tree from organism_core, IEDB taxa, and active NCBI taxo.
# Determine all the active species.
# Build the subspecies tree by adding all descendants of active species
# from the full NCBI Taxonomy.

# TODO: Check that this is a reasonable way to count.
build/arborist/peptide-count.tsv: src/organism/peptide-count.sql build/iedb/peptide.built | build/iedb/nanobot.db
	src/util/sqlite2tsv $| $< $@

# Render the organism_core as HTML.
build/arborist/organism_core.html: src/organism/render_organism_core.py src/organism/organism_core.tsv | build/arborist/
	python3 $^ $@

# Build a new organism tree.
build/arborist/organism-tree.tsv: src/organism/assign_species.py build/arborist/ncbitaxon.built src/organism/organism_core.tsv build/iedb/ncbi_include.tsv build/iedb/iedb_taxa.tsv build/arborist/peptide-count.tsv | build/arborist/nanobot.db
	python3 $< $| $(filter %.tsv, $^) $@
	qsv sort $@ --output $@

# Convert the organism tree to an LDTab table in SQLite.
build/arborist/organism-tree.built: src/organism/build_organism_tree.py build/arborist/organism-tree.tsv | build/arborist/nanobot.db
	sqlite3 $| "DROP TABLE IF EXISTS organism_tree"
	python3 $< $| $(filter %.tsv, $^)
	touch $@

build/arborist/subspecies-tree.built: src/organism/build_subspecies_tree.py build/arborist/organism-tree.built | build/arborist/nanobot.db
	sqlite3 $| "DROP TABLE IF EXISTS subspecies_tree"
	python3 $< $|
	touch $@

# Export sorted LDTab TSV files.
build/arborist/%-tree-ldtab.tsv: build/arborist/%-tree.built | build/arborist/nanobot.db
	rm -f $@
	ldtab export $| --table $*_tree $@
	mv $@ $@-tmp.tsv
	qsv sort $@-tmp.tsv --output $@
	rm $@-tmp.tsv

build/arborist/%-tree.ttl: build/arborist/%-tree.built | build/arborist/nanobot.db
	rm -f $@
	ldtab export $| --table $*_tree $@

build/arborist/%-tree.owl: build/arborist/%-tree.ttl src/organism/predicates.ttl
	robot merge --input $< --input $(word 2,$^) \
	annotate \
	--ontology-iri https://ontology.iedb.org/ontology/$(notdir $@) \
	--output $@

build/arborist/active-species.tsv: src/organism/get_active_species.py build/arborist/organism-tree.built build/arborist/peptide-count.tsv | build/arborist/nanobot.db
	python3 $< $| $(filter %.tsv, $^) $@

.PHONY: organism
organism: build/arborist/organism_core.html
organism: build/arborist/organism-tree-ldtab.tsv build/arborist/subspecies-tree-ldtab.tsv
organism: build/arborist/organism-tree.owl build/arborist/subspecies-tree.owl
organism: build/arborist/active-species.tsv build/arborist/proteome.tsv
	make reload

build/organisms/latest/: build/arborist/subspecies-tree.owl
	rm -rf $@
	mkdir -p $@
	cp $^ $@
	chmod 644 $@*

### 4. Select Proteomes for each active species
#
# Create a directory for any active species.
# Use the active_taxa column to get all descendant taxa.
# Copy peptides and sources into the directory.
# Select a proteome for that species,
# fetching FASTA and XML annotations.

# TODO: Use previously selected proteomes or force refresh.
build/arborist/proteome.tsv: build/arborist/active-species.tsv src/proteome/proteome.tsv
	qsv join --left 'Species ID' $< 'Species ID' $(word 2,$^) \
	| qsv select 1-6,12- --output $@

build/species/%/:
	mkdir -p $@

# Get active taxa list as a regular expression pattern of alternates,
# for use by `qsv select`, e.g. 1053|11057|11059|11060|...
build/species/%/taxa.txt: build/arborist/active-species.tsv | build/species/%/
	awk 'BEGIN {FS="\t"} $$2==$* {print $$4}' $< | \
	sed 's/, /|/g' > $@

build/species/%/epitopes.tsv: build/iedb/peptide.tsv build/species/%/taxa.txt
	qsv search --select 'Organism ID' `cat build/species/$*/taxa.txt` $< --output $@

build/species/%/sources.tsv: build/iedb/peptide_source.tsv build/species/%/taxa.txt
	qsv search --select 'Organism ID' `cat build/species/$*/taxa.txt` $< --output $@

.PRECIOUS: build/species/%/epitopes.tsv build/species/%/sources.tsv

build/arborist/proteomes.built: build/arborist/proteome.tsv
	python3 src/protein_tree/protein_tree/select_proteome.py -b build/ -a

# Remove epitope counts from proteome table.
build/arborist/proteome_after.tsv: build/arborist/proteome.tsv
	qsv select 1-5,7- $< --output $@

# Compare new proteomes to src/proteome.
build/arborist/proteome.html: src/proteome/proteome.tsv build/arborist/proteome_after.tsv
	daff --output $@ $^

.PHONY: proteome
proteome: build/arborist/proteomes.built


### 5. Build Protein Tree

build/arborist/allergens.csv: | build/
	curl -L -o $@ 'http://www.allergen.org/csv.php?table=joint'

build/arborist/allergens.tsv: src/util/csv2tsv.py build/arborist/allergens.csv
	python3 $^ $@

build/arborist/manual-parents.tsv: build/arborist/allergens.tsv
	# wget --no-check-certificate 'https://docs.google.com/spreadsheets/d/1VUDYmmnQURRnuqyVxZGyF8JCgAIiooaKCmi3_mf03o8/export?format=tsv&gid=2087231134' -O $@
	cp src/protein/data/manual-parents.tsv $@

# build/arborist/protein-tree.assigned: build/arborist/allergens.tsv build/arborist/manual-parents.tsv
# 	python src/protein/protein_tree/assign.py -n 14
# 	touch $@

build/arborist/mro.owl: | build/arborist/
	curl -L -o $@ 'http://purl.obolibrary.org/obo/mro/2023-12-21/mro.owl'

build/arborist/mro.built: build/arborist/mro.owl | build/arborist/nanobot.db
	$(eval DB := build/arborist/nanobot.db)
	$(eval TABLE := mro)
	sqlite3 $(DB) 'DROP TABLE IF EXISTS $(TABLE)'
	ldtab init $(DB) --table $(TABLE)
	ldtab import $(DB) $< --table $(TABLE)
	sqlite3 $(DB) 'CREATE INDEX idx_$(TABLE)_subject ON $(TABLE)(subject)'
	sqlite3 $(DB) 'CREATE INDEX idx_$(TABLE)_predicate ON $(TABLE)(predicate)'
	sqlite3 $(DB) 'CREATE INDEX idx_$(TABLE)_object ON $(TABLE)(object)'
	sqlite3 $(DB) 'ANALYZE $(TABLE)'
	touch $@

.PHONY: mro
mro: build/arborist/mro.built

build/arborist/protein-tree.built: src/protein/protein_tree/build.py build/arborist/all-peptide-assignments.tsv build/arborist/organism-tree.built
	$(eval DB := build/arborist/nanobot.db)
	sqlite3 $(DB) 'DROP TABLE IF EXISTS protein_tree_old'
	sqlite3 $(DB) 'DROP TABLE IF EXISTS protein_tree_new'
	python $<
	sqlite3 $(DB) 'CREATE INDEX idx_protein_tree_old_subject ON protein_tree_old(subject)'
	sqlite3 $(DB) 'CREATE INDEX idx_protein_tree_old_predicate ON protein_tree_old(predicate)'
	sqlite3 $(DB) 'CREATE INDEX idx_protein_tree_old_object ON protein_tree_old(object)'
	sqlite3 $(DB) 'ANALYZE protein_tree_old'
	sqlite3 $(DB) 'CREATE INDEX idx_protein_tree_new_subject ON protein_tree_new(subject)'
	sqlite3 $(DB) 'CREATE INDEX idx_protein_tree_new_predicate ON protein_tree_new(predicate)'
	sqlite3 $(DB) 'CREATE INDEX idx_protein_tree_new_object ON protein_tree_new(object)'
	sqlite3 $(DB) 'ANALYZE protein_tree_new'
	touch $@

build/arborist/protein-tree.ttl: build/arborist/protein-tree.built | build/arborist/nanobot.db
	rm -f $@
	ldtab export $| $@ --table protein_tree_old

build/arborist/protein-tree.owl: build/arborist/protein-tree.ttl
	robot convert -i $< -o $@

build/arborist/epitope-mappings_new.tsv: build/arborist/epitope-mappings.tsv
	qsv slice --start -10 $< --output $@

.PHONY: protein
protein: build/arborist/protein-tree.owl


### 6. TODO Build Molecule Tree

build/arborist/molecule-tree.owl: nonpeptide-tree-20240305.owl build/arborist/protein-tree.owl
	robot remove \
	--input $< \
	--term BFO:0000023 \
	--select "self descendants" \
	merge \
	--input $(word 2,$^) \
	annotate \
	--ontology-iri https://ontology.iedb.org/ontology/molecule-tree.owl \
	--version-iri https://ontology.iedb.org/ontology/$(shell date +%Y-%m-%d)/molecule-tree.owl \
	--output $@

build/arborist/molecule-tree.built: build/arborist/molecule-tree.owl
	$(eval DB := build/arborist/nanobot.db)
	$(eval TABLE := molecule_tree)
	sqlite3 $(DB) 'DROP TABLE IF EXISTS $(TABLE)'
	ldtab init $(DB) --table $(TABLE)
	ldtab import $(DB) $< --table $(TABLE)
	sqlite3 $(DB) 'CREATE INDEX idx_$(TABLE)_subject ON $(TABLE)(subject)'
	sqlite3 $(DB) 'CREATE INDEX idx_$(TABLE)_predicate ON $(TABLE)(predicate)'
	sqlite3 $(DB) 'CREATE INDEX idx_$(TABLE)_object ON $(TABLE)(object)'
	sqlite3 $(DB) 'ANALYZE $(TABLE)'
	touch $@

.PHONY: molecule
molecule: build/arborist/molecule-tree.built

build/proteins/latest/: build/disease-tree.owl build/arborist/molecule-tree.owl build/arborist/parent-proteins.tsv build/arborist/source-parents.tsv build/arborist/epitope-mappings.tsv build/arborist/epitope-mappings_new.tsv
	rm -rf $@
	mkdir -p $@
	cp $^ $@
	chmod 644 $@*


### 7. TODO Build Assay Tree


### 8. Build Disease Tree

# Copy disease-tree.owl from production.
# TODO: Build proper disease tree
build/disease-tree.owl: | build/
	curl -L -k -o $@ 'https://10.0.7.92/proteins/latest/disease-tree.owl'

build/disease-tree.tsv: build/disease-tree.owl | build/arborist/nanobot.db
	sqlite3 $| 'DROP TABLE IF EXISTS disease_tree'
	ldtab init $| --table disease_tree
	ldtab import $| $< --table disease_tree
	sqlite3 $| 'CREATE INDEX idx_disease_tree_subject ON disease_tree(subject)'
	sqlite3 $| 'CREATE INDEX idx_disease_tree_predicate ON disease_tree(predicate)'
	sqlite3 $| 'CREATE INDEX idx_disease_tree_object ON disease_tree(object)'
	sqlite3 $| 'ANALYZE disease_tree'
	ldtab export $| $@ --table disease_tree

.PHONY: disease
disease: build/disease-tree.tsv

# TODO: geolocation tree, MHC tree
# TODO: merged SoT tree
# TODO: test data, symlink?


### Nanobot Actions
#
# Editing operations for the Nanobot Arborist database.

.PHONY: save
save: | build/arborist/nanobot.db
	valve-export data $| src/arborist/ table column datatype
	valve-export data $| src/organism/ organism_core
	valve-export data $| build/arborist/ $$(grep build src/arborist/table.tsv | cut -f1 | tr '\n' ' ')
	python3 src/organism/sort_organism_core.py src/organism/organism_core.tsv

DROPTABLES := proteomes active_species organism_core organism_tree_tsv prefix column datatype table message history
.PHONY: reload
reload: src/organism/check_organism_core.py | build/arborist/nanobot.db
	sqlite3 $| $(foreach DT,$(DROPTABLES),"DROP VIEW IF EXISTS '$(DT)_text_view'" "DROP VIEW IF EXISTS '$(DT)_view'" "DROP TABLE IF EXISTS '$(DT)_conflict'" "DROP TABLE IF EXISTS '$(DT)'")
	cd $(dir $|) && nanobot init
	-python3 $< $|


### Comparisons
#
# These tasks build other trees for comparison.

# Load an existing organism-tree.owl
build/arborist/organism-tree-old.built: build/arborist/organism-tree-old.owl | build/arborist/nanobot.db
	sqlite3 $| 'DROP TABLE IF EXISTS organism_tree_old'
	ldtab init $| --table organism_tree_old
	ldtab import $| $(word 2,$^) -t organism_tree_old
	sqlite3 $| 'ANALYZE organism_tree_old'
	touch $@

# Load an existing subspecies-tree.owl
build/arborist/subspecies-tree-old.built: build/arboirst/subspecies-tree-old.owl | build/arborist/nanobot.db
	sqlite3 $| 'DROP TABLE IF EXISTS subspecies_tree_old'
	ldtab init $| --table subspecies_tree_old
	ldtab import $| $< --table subspecies_tree_old
	sqlite3 $| 'ANALYZE subspecies_tree_old'
	touch $@

# Load an existing molecule-tree.owl
build/arborist/molecule-tree-old.built: molecule-tree-20240317.owl | build/arborist/nanobot.db
	$(eval TABLE := molecule_tree_old)
	sqlite3 $| 'DROP TABLE IF EXISTS $(TABLE)'
	ldtab init $| --table $(TABLE)
	ldtab import $| $< --table $(TABLE)
	sqlite3 $| 'CREATE INDEX idx_$(TABLE)_subject ON $(TABLE)(subject)'
	sqlite3 $| 'CREATE INDEX idx_$(TABLE)_predicate ON $(TABLE)(predicate)'
	sqlite3 $| 'CREATE INDEX idx_$(TABLE)_object ON $(TABLE)(object)'
	sqlite3 $| 'ANALYZE $(TABLE)'
	touch $@

