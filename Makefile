### Arborist Makefile
#
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
# 7. Build Molecule Tree
# 8. Build Assay Tree
# 9. Build Disease Tree
#
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
.DEFAULT_GOAL := usage
.DELETE_ON_ERROR:
.SUFFIXES:

export PATH := $(shell pwd)/bin:$(PATH)


### Main Tasks
#
# The main tasks for running Arborist.

.PHONY: usage
usage:
	@echo "Arborist: build trees for the IEDB"
	@echo ""
	@echo "TASKS"
	@echo "  deps     install dependencies"
	@echo "  iedb     load IEDB data"
	@echo "  all      build all trees"
	@echo "  serve    run web interface on localhost:3000"
	@echo "  clean    remove all build files"
	@echo "  clobber  remove all generated files"
	@echo "  usage    print this message"

# Dependencies are added to this list below.
.PHONY: deps
deps:

.PHONY: all
all: build/organism-tree.owl build/subspecies-tree.owl build/active-species.tsv

.PHONY: serve
serve: build/iedb/nanobot.db
	cd build/iedb/ && nanobot serve

.PHONY: clean
clean:
	rm -rf build

.PHONY: clobber
clobber:
	rm -rf bin build cache

bin/ build/ cache/:
	mkdir -p $@


### Install Dependencies
#
# For each software dependency
# we use Make's `ifeq` macro and `command -v <name>`
# to check if the dependency is already on the PATH
# (including `bin/`).
# If not, then we define a Make task to install it to `bin/`,
# and add that dependency to `deps`.

# Require SQLite
ifeq (, $(shell command -v sqlite3))
$(error 'Please install SQLite')
endif

# Require Java
ifeq (, $(shell command -v java))
$(error 'Please install Java, so we can run ROBOT and LDTab')
endif

# Install ROBOT if not already present
ifeq (, $(shell command -v robot))
bin/robot.jar: | bin/
	curl -L -o $@ 'https://github.com/ontodev/robot/releases/download/v1.9.4/robot.jar'
bin/robot: bin/robot.jar
	curl -L -o $@ 'https://raw.githubusercontent.com/ontodev/robot/master/bin/robot'
	chmod +x $@
deps: bin/robot
endif

# Install LDTab if not already present
ifeq (, $(shell command -v ldtab))
bin/ldtab.jar: | bin/
	curl -L -o $@ 'https://github.com/ontodev/ldtab.clj/releases/download/v2023-08-19/ldtab.jar'
bin/ldtab: bin/ldtab.jar
	echo '#!/bin/sh' > $@
	echo 'java -jar build/ldtab.jar' >> $@
	chmod +x $@
deps: bin/ldtab
endif

# Install Nanobot if not already present
ifeq (, $(shell command -v nanobot))
bin/nanobot: | bin/
	curl -L -k -o $@ 'https://github.com/ontodev/nanobot.rs/releases/download/v2023-10-26/nanobot-x86_64-unknown-linux-musl'
	chmod +x $@
deps: bin/nanobot
endif

# Install valve-export script if not already present
ifeq (, $(shell command -v valve-export))
bin/valve-export: | bin/
	curl -L -o $@ 'https://github.com/ontodev/valve.rs/raw/main/scripts/export.py'
	chmod +x $@
deps: valve-export
endif

# Install valve-export script if not already present
ifeq (, $(shell command -v qsv))
bin/qsv: | bin/ build/
	curl -L -k -o build/qsv.zip 'https://github.com/jqnatividad/qsv/releases/download/0.118.0/qsv-0.118.0-x86_64-unknown-linux-musl.zip'
	cd build && unzip qsv.zip qsv
	mv build/qsv $@
deps: bin/qsv
endif

# Install BLAST if not already present
ifeq (, $(shell command -v blastp))
BLAST_VERSION := 2.15.0
build/ncbi-blast.tar.gz: | build/
	curl -L -o $@ 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(BLAST_VERSION)/ncbi-blast-$(BLAST_VERSION)+-x64-linux.tar.gz'
bin/blastp bin/makeblastdb: build/ncbi-blast.tar.gz | bin/
	cd build/ && tar -zxvf $(notdir $<) ncbi-blast-$(BLAST_VERSION)+/$@
	mv build/ncbi-blast-$(BLAST_VERSION)+/$@ $@
deps: bin/blastp bin/makeblastdb
endif

# Install hmmer if not already present
ifeq (, $(shell command -v hmmscan))
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


### 1. Fetch IEDB Data
#
# Copy some tables from the IEDB,
# then validate and load them into a SQLite database.
# Requires MySQL/MariaDB client and connection details.

build/iedb/:
	mkdir -p $@

build/iedb/nanobot.toml: src/iedb/nanobot.toml | build/iedb/
	cp $< $@

build/iedb/table.tsv: src/iedb/table.tsv | build/iedb/
	cp $< $@

build/iedb/column.tsv: src/iedb/column.tsv | build/iedb/
	cp $< $@

build/iedb/datatype.tsv: src/iedb/datatype.tsv | build/iedb/
	cp $< $@

.PRECIOUS: current/iedb/%
current/iedb/%:
	@echo "Updating cached IEDB data..."
	src/iedb/update-cache

# TODO: This task will be much simpler with planned VALVE features.
build/iedb/nanobot.db: build/iedb/nanobot.toml
build/iedb/nanobot.db: build/iedb/datatype.tsv
build/iedb/nanobot.db: build/iedb/table.tsv
build/iedb/nanobot.db: build/iedb/column.tsv
build/iedb/nanobot.db: current/iedb/
build/iedb/nanobot.db: current/iedb/ncbi_include.tsv.gz
build/iedb/nanobot.db: current/iedb/iedb_taxa.tsv.gz
build/iedb/nanobot.db: current/iedb/source.tsv.gz
build/iedb/nanobot.db: current/iedb/object.tsv.gz
build/iedb/nanobot.db: current/iedb/epitope.tsv.gz
	rm -f $@
	zcat current/iedb/ncbi_include.tsv.gz > build/iedb/ncbi_include.tsv
	zcat current/iedb/iedb_taxa.tsv.gz > build/iedb/iedb_taxa.tsv
	zcat current/iedb/source.tsv.gz | head -n1 > build/iedb/source.tsv || exit 0
	zcat current/iedb/object.tsv.gz | head -n1 > build/iedb/object.tsv || exit 0
	zcat current/iedb/epitope.tsv.gz | head -n1 > build/iedb/epitope.tsv || exit 0
	echo 'Accession	Name	Sequence	Organism ID' > build/iedb/peptide_source.tsv
	echo 'Sequence	Source Name	Accession	Organism ID' > build/iedb/peptide.tsv
	cd build/iedb/ && nanobot init
	zcat current/iedb/source.tsv.gz > build/iedb/source.tsv
	zcat current/iedb/object.tsv.gz > build/iedb/object.tsv
	zcat current/iedb/epitope.tsv.gz > build/iedb/epitope.tsv
	src/util/tsv2sqlite $@ build/iedb/source.tsv
	src/util/tsv2sqlite $@ build/iedb/object.tsv
	src/util/tsv2sqlite $@ build/iedb/epitope.tsv

build/iedb/peptide.tsv: build/iedb/nanobot.db src/iedb/peptide.sql
	src/util/sqlite2tsv $^ $@
	src/util/tsv2sqlite $< $@

build/iedb/peptide_source.tsv: build/iedb/nanobot.db src/iedb/peptide_source.sql
	src/util/sqlite2tsv $^ $@
	src/util/tsv2sqlite $< $@

.PHONY: iedb
iedb: build/iedb/peptide.tsv build/iedb/peptide_source.tsv


# TODO: Replace with local queries

# Tax ID -> epitope count
build/counts.tsv: src/get-counts.sql | build/
	$(MIRROR_QUERY) < $< > $@

build/ncbi_include.tsv: src/ncbi-active-taxa.sql | build/
	$(MIRROR_QUERY) < $< > $@

build/counts_full.tsv: src/combine_taxids.py build/counts.tsv build/ncbi_include.tsv
	python3 $^ $@

# Custom IEDB taxa
build/iedb_taxa.tsv: src/get-iedb-taxa.sql | build/
	$(MIRROR_QUERY) < $< > $@


### TODO: Find a home

build/nanobot.db: | build/nanobot
	rm -f $@ build/*.built
	echo 'curie	label	label_source	iedb_synonyms	rank	level	epitope_count	parent	parent_label	parent2	parent2_label	species	species_label	source_table	use_other' > build/organism-tree.tsv
	nanobot init

.PHONY: save
save: build/nanobot.db
	valve-export data build/nanobot.db src/schema/ table column datatype
	valve-export data build/nanobot.db src/ organism_core
	valve-export data build/nanobot.db build/ $$(grep build src/schema/table.tsv | cut -f1 | tr '\n' ' ')
	python3 src/sort_organism_core.py src/organism_core.tsv

DROPTABLES := epitope_object proteomes active_species organism_core organism_tree_tsv iedb_taxa prefix column datatype table message history
.PHONY: reload
reload: src/check_organism_core.py | build/nanobot.db
	sqlite3 build/nanobot.db $(foreach DT,$(DROPTABLES),"DROP VIEW IF EXISTS '$(DT)_view'" "DROP TABLE IF EXISTS '$(DT)_conflict'" "DROP TABLE IF EXISTS '$(DT)'")
	nanobot init
	-python3 $< build/nanobot.db


### 3. Build NCBI Taxonomy

build/taxdmp.zip: | build/
	# curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2023-11-01.zip

build/ncbitaxon.built: build/taxdmp.zip | build/nanobot.db
	sqlite3 build/nanobot.db "DROP TABLE IF EXISTS ncbitaxon"
	python3 src/ncbitaxon2ldtab.py $< build/nanobot.db
	sqlite3 build/nanobot.db "CREATE INDEX idx_ncbitaxon_subject ON ncbitaxon(subject)"
	sqlite3 build/nanobot.db "CREATE INDEX idx_ncbitaxon_predicate ON ncbitaxon(predicate)"
	sqlite3 build/nanobot.db "CREATE INDEX idx_ncbitaxon_object ON ncbitaxon(object)"
	sqlite3 build/nanobot.db "ANALYZE ncbitaxon"
	touch $@


### 3. Build Organism Tree and determine active species

build/organism_core.html: src/build_organism_core.py src/organism_core.tsv
	python3 $^ $@

build/organism-tree.tsv: build/ncbitaxon.built src/assign_species.py src/organism_core.tsv build/iedb_taxa.tsv build/counts_full.tsv
	python3 $(word 2, $^) build/nanobot.db $(filter %.tsv, $^) $@
	qsv sort $@ --output $@

# Build a new organism tree
build/organism-tree.built: build/ncbitaxon.built src/build_organism_tree.py build/organism-tree.tsv build/organism_core.html
	sqlite3 build/nanobot.db "DROP TABLE IF EXISTS organism_tree"
	python3 $(word 2, $^) build/nanobot.db $(filter %.tsv, $^)
	sqlite3 build/nanobot.db "ANALYZE organism_tree"
	touch $@

build/subspecies-tree.built: build/organism-tree.built src/build_subspecies_tree.py
	sqlite3 build/nanobot.db "DROP TABLE IF EXISTS subspecies_tree"
	python3 $(word 2, $^) build/nanobot.db
	sqlite3 build/nanobot.db "ANALYZE subspecies_tree"
	touch $@

build/organism-tree.ttl: build/organism-tree.built | build/ldtab.jar
	rm -f $@
	ldtab export build/nanobot.db --table organism_tree $@

build/subspecies-tree.ttl: build/subspecies-tree.built | build/ldtab.jar
	rm -f $@
	ldtab export build/nanobot.db --table subspecies_tree $@

build/organism-tree.sorted.tsv: build/organism-tree.built | build/ldtab.jar
	rm -f build/organism-tree.ldtab.tsv
	ldtab export build/nanobot.db --table organism_tree build/organism-tree.ldtab.tsv
	cut -f4- build/organism-tree.ldtab.tsv | sort > $@

build/%-tree.owl: build/%-tree.ttl src/predicates.ttl | build/robot.jar
	robot merge --input $< --input $(word 2,$^) \
	annotate \
	--ontology-iri https://ontology.iedb.org/ontology/$(notdir $@) \
	--output $@

build/active-species.tsv: src/get_active_species.py build/organism-tree.built build/counts_full.tsv
	python3 $< build/nanobot.db build/counts_full.tsv $@

build/proteomes.tsv: build/active-species.tsv src/selected_proteomes.tsv
	qsv join --left 'Species ID' $< 'Species ID' $(word 2,$^) \
	| qsv select 1-6,12- --output $@


### 4. Select Proteomes for each active species

build/allergens.csv: | build/
	curl -L -o $@ 'http://www.allergen.org/csv.php?table=joint'

build/allergens.tsv: build/allergens.csv
	qsv input $< --output $@

build/%/:
	mkdir -p $@

build/%/taxa.txt: build/active-species.tsv | build/%/
	awk 'BEGIN {FS="\t"} $$2==$* {print $$4}' build/active-species.tsv | \
	sed 's/, /|/g' > $@

build/%/epitopes.tsv: build/epitope_object.tsv build/%/taxa.txt
	qsv search --select 'Organism ID' `cat build/$*/taxa.txt` $< --output $@

build/%/sources.tsv: build/sources.tsv build/%/taxa.txt
	qsv search --select 'Organism ID' `cat build/$*/taxa.txt` $< --output $@

build/%/proteome.tsv: build/%/epitopes.tsv build/%/sources.tsv
	protein_tree/protein_tree/select_proteome.py -t $*


### 5. Assign Proteins

build/%/epitope_assignments.tsv: build/%/epitopes.tsv build/%/sources.tsv build/%/proteome.tsv
	protein_tree/protein_tree/run.py -t $*

.PHONY: dengue
dengue: build/12637/epitope_assignments.tsv

### 6. Build Protein Tree
### 7. Build Molecule Tree


### Comparisons
#
# These tasks build other trees for comparison.

# Load an existing organism-tree.owl
build/organism-tree-old.built: build/ncbitaxon.built build/organism-tree-old.owl | build/ldtab.jar
	sqlite3 build/nanobot.db "DROP TABLE IF EXISTS organism_tree_old"
	sed s/statement/organism_tree_old/g src/statement.sql | sqlite3 build/nanobot.db
	ldtab import build/nanobot.db $(word 2,$^) -t organism_tree_old
	sqlite3 build/nanobot.db "ANALYZE organism_tree_old"
	touch $@

# Load an existing subspecies-tree.owl
build/subspecies-tree-old.built: build/ncbitaxon.built build/subspecies-tree-old.owl | build/ldtab.jar
	sqlite3 build/nanobot.db "DROP TABLE IF EXISTS subspecies_tree_old"
	sed s/statement/subspecies_tree_old/g src/statement.sql | sqlite3 build/nanobot.db
	ldtab import build/nanobot.db $(word 2,$^) -t subspecies_tree_old
	sqlite3 build/nanobot.db "ANALYZE subspecies_tree_old"
	touch $@


