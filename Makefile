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

export PATH := bin:$(PATH)


### Main Tasks
#
# The main tasks for running Arborist.

.PHONY: usage
usage:
	@echo "Arborist: build trees for the IEDB"
	@echo ""
	@echo "TASKS"
	@echo "  deps     install dependencies"
	@echo "  all      build all trees"
	@echo "  clean    remove all build files"
	@echo "  clobber  remove all cached files"
	@echo "  usage    print this message"

.PHONY: all
all: build/organism-tree.owl build/subspecies-tree.owl build/active-species.tsv

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

.PHONY: deps
deps:

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

# TODO: This task will be much simpler with planned VALVE features.
build/iedb/nanobot.db: src/tsv2sqlite
build/iedb/nanobot.db: build/iedb/nanobot.toml
build/iedb/nanobot.db: build/iedb/datatype.tsv
build/iedb/nanobot.db: build/iedb/table.tsv
build/iedb/nanobot.db: build/iedb/column.tsv
build/iedb/nanobot.db: cache/iedb/ncbi_include.tsv
build/iedb/nanobot.db: cache/iedb/iedb_taxa.tsv
build/iedb/nanobot.db: cache/iedb/source.tsv
build/iedb/nanobot.db: cache/iedb/object.tsv
build/iedb/nanobot.db: cache/iedb/epitope.tsv
build/iedb/nanobot.db: | build/nanobot
	rm -f $@
	head -n1 cache/iedb/source.tsv > build/iedb/source.tsv
	head -n1 cache/iedb/object.tsv > build/iedb/object.tsv
	head -n1 cache/iedb/epitope.tsv > build/iedb/epitope.tsv
	echo 'Accession	Name	Sequence	Organism ID' > build/iedb/peptide_source.tsv
	echo 'Sequence	Source Name	Accession	Organism ID' > build/iedb/peptide.tsv
	cd build/iedb/ && ../../$(NANOBOT) init
	cp cache/iedb/source.tsv build/iedb/source.tsv
	cp cache/iedb/object.tsv build/iedb/object.tsv
	cp cache/iedb/epitope.tsv build/iedb/epitope.tsv
	$< $@ build/iedb/source.tsv
	$< $@ build/iedb/object.tsv
	$< $@ build/iedb/epitope.tsv

build/iedb/peptide.tsv: src/sqlite2tsv build/iedb/nanobot.db src/peptide.sql | build/qsv
	$^ $@
	src/tsv2sqlite build/iedb/nanobot.db $@

build/iedb/peptide_source.tsv: src/sqlite2tsv build/iedb/nanobot.db src/peptide_source.sql | build/qsv
	$^ $@
	src/tsv2sqlite build/iedb/nanobot.db $@

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

$(DB): | build/nanobot
	rm -f $@ build/*.built
	echo 'curie	label	label_source	iedb_synonyms	rank	level	epitope_count	parent	parent_label	parent2	parent2_label	species	species_label	source_table	use_other' > build/organism-tree.tsv
	$(NANOBOT) init

.PHONY: save
save: $(EXPORT) $(DB)
	python3 $(EXPORT) data $(DB) src/schema/ table column datatype
	python3 $(EXPORT) data $(DB) src/ organism_core
	python3 $(EXPORT) data $(DB) build/ $$(grep build src/schema/table.tsv | cut -f1 | tr '\n' ' ')
	python3 src/sort_organism_core.py src/organism_core.tsv

DROPTABLES := epitope_object proteomes active_species organism_core organism_tree_tsv iedb_taxa prefix column datatype table message history
.PHONY: reload
reload: src/check_organism_core.py | $(DB)
	sqlite3 $(DB) $(foreach DT,$(DROPTABLES),"DROP VIEW IF EXISTS '$(DT)_view'" "DROP TABLE IF EXISTS '$(DT)_conflict'" "DROP TABLE IF EXISTS '$(DT)'")
	$(NANOBOT) init
	-python3 $< $(DB)


### 3. Build NCBI Taxonomy

build/taxdmp.zip: | build/
	# curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2023-11-01.zip

build/ncbitaxon.built: build/taxdmp.zip | $(DB)
	sqlite3 $(DB) "DROP TABLE IF EXISTS ncbitaxon"
	python3 src/ncbitaxon2ldtab.py $< $(DB)
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_subject ON ncbitaxon(subject)"
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_predicate ON ncbitaxon(predicate)"
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_object ON ncbitaxon(object)"
	sqlite3 $(DB) "ANALYZE ncbitaxon"
	touch $@


### 3. Build Organism Tree and determine active species

build/organism_core.html: src/build_organism_core.py src/organism_core.tsv
	python3 $^ $@

build/organism-tree.tsv: build/ncbitaxon.built src/assign_species.py src/organism_core.tsv build/iedb_taxa.tsv build/counts_full.tsv
	python3 $(word 2, $^) $(DB) $(filter %.tsv, $^) $@
	${QSV} sort $@ --output $@

# Build a new organism tree
build/organism-tree.built: build/ncbitaxon.built src/build_organism_tree.py build/organism-tree.tsv build/organism_core.html
	sqlite3 $(DB) "DROP TABLE IF EXISTS organism_tree"
	python3 $(word 2, $^) $(DB) $(filter %.tsv, $^)
	sqlite3 $(DB) "ANALYZE organism_tree"
	touch $@

build/subspecies-tree.built: build/organism-tree.built src/build_subspecies_tree.py
	sqlite3 $(DB) "DROP TABLE IF EXISTS subspecies_tree"
	python3 $(word 2, $^) $(DB)
	sqlite3 $(DB) "ANALYZE subspecies_tree"
	touch $@

build/organism-tree.ttl: build/organism-tree.built | build/ldtab.jar
	rm -f $@
	$(LDTAB) export $(DB) --table organism_tree $@

build/subspecies-tree.ttl: build/subspecies-tree.built | build/ldtab.jar
	rm -f $@
	$(LDTAB) export $(DB) --table subspecies_tree $@

build/organism-tree.sorted.tsv: build/organism-tree.built | build/ldtab.jar
	rm -f build/organism-tree.ldtab.tsv
	$(LDTAB) export $(DB) --table organism_tree build/organism-tree.ldtab.tsv
	cut -f4- build/organism-tree.ldtab.tsv | sort > $@

build/%-tree.owl: build/%-tree.ttl src/predicates.ttl | build/robot.jar
	$(ROBOT) merge --input $< --input $(word 2,$^) \
	annotate \
	--ontology-iri https://ontology.iedb.org/ontology/$(notdir $@) \
	--output $@

build/active-species.tsv: src/get_active_species.py build/organism-tree.built build/counts_full.tsv
	python3 $< $(DB) build/counts_full.tsv $@

build/proteomes.tsv: build/active-species.tsv src/selected_proteomes.tsv | $(QSV)
	$(QSV) join --left 'Species ID' $< 'Species ID' $(word 2,$^) \
	| $(QSV) select 1-6,12- --output $@


### 4. Select Proteomes for each active species

build/allergens.csv: | build/
	curl -L -o $@ 'http://www.allergen.org/csv.php?table=joint'

build/allergens.tsv: build/allergens.csv
	$(QSV) input $< --output $@

build/%/:
	mkdir -p $@

build/%/taxa.txt: build/active-species.tsv | build/%/
	awk 'BEGIN {FS="\t"} $$2==$* {print $$4}' build/active-species.tsv | \
	sed 's/, /|/g' > $@

build/%/epitopes.tsv: build/epitope_object.tsv build/%/taxa.txt
	$(QSV) search --select 'Organism ID' `cat build/$*/taxa.txt` $< --output $@

build/%/sources.tsv: build/sources.tsv build/%/taxa.txt
	$(QSV) search --select 'Organism ID' `cat build/$*/taxa.txt` $< --output $@

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
	sqlite3 $(DB) "DROP TABLE IF EXISTS organism_tree_old"
	sed s/statement/organism_tree_old/g src/statement.sql | sqlite3 $(DB)
	$(LDTAB) import $(DB) $(word 2,$^) -t organism_tree_old
	sqlite3 $(DB) "ANALYZE organism_tree_old"
	touch $@

# Load an existing subspecies-tree.owl
build/subspecies-tree-old.built: build/ncbitaxon.built build/subspecies-tree-old.owl | build/ldtab.jar
	sqlite3 $(DB) "DROP TABLE IF EXISTS subspecies_tree_old"
	sed s/statement/subspecies_tree_old/g src/statement.sql | sqlite3 $(DB)
	$(LDTAB) import $(DB) $(word 2,$^) -t subspecies_tree_old
	sqlite3 $(DB) "ANALYZE subspecies_tree_old"
	touch $@


