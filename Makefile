### Configuration
#
# These are standard options to make Make sane:
# <http://clarkgrubb.com/makefile-style-guide#toc2>

MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.SUFFIXES:

all: build/organism-tree.owl build/subspecies-tree.owl build/active-species.tsv

clean:
	rm -rf build


### Set Up

build:
	mkdir -p $@

ROBOT := java -jar build/robot.jar --prefix "iedb-taxon: https://ontology.iedb.org/taxon/" --prefix "ONTIE: https://ontology.iedb.org/ONTIE_"
LDTAB := java -jar build/ldtab.jar
NANOBOT := build/nanobot
QSV := build/qsv
EXPORT := build/export.py
DB := build/nanobot.db

build/robot.jar: | build
	curl -L -o $@ "https://github.com/ontodev/robot/releases/download/v1.9.4/robot.jar"

build/ldtab.jar: | build
	curl -L -o $@ "https://github.com/ontodev/ldtab.clj/releases/download/v2023-08-19/ldtab.jar"

build/nanobot: | build
	curl -L -k -o $@ "https://github.com/ontodev/nanobot.rs/releases/download/v2023-10-26/nanobot-x86_64-unknown-linux-musl"
	chmod +x $@

# Download qsv binary for Linux ARM64
build/qsv: | build/
	curl -L -k -o build/qsv.zip "https://github.com/jqnatividad/qsv/releases/download/0.112.0/qsv-0.112.0-x86_64-unknown-linux-musl.zip"
	cd build && unzip qsv.zip qsv

build/export.py: | build/
	curl -L -o $@ "https://github.com/ontodev/valve.rs/raw/main/scripts/export.py"


### Nanobot Database

$(DB): | build build/nanobot
	rm -f $@ build/*.built
	$(NANOBOT) init

.PHONY: save
save: $(EXPORT) $(DB)
	python3 $(EXPORT) data $(DB) src/schema/ table column datatype
	python3 $(EXPORT) data $(DB) src/ organism_core
	python3 $(EXPORT) data $(DB) build/ $$(grep build src/schema/table.tsv | cut -f1 | tr '\n' ' ')
	python3 src/sort_organism_core.py src/organism_core.tsv

DROPTABLES := proteomes active_species organism_core organism_tree_tsv iedb_taxa prefix column datatype table message
.PHONY: reload
reload: src/check_organism_core.py | $(DB)
	sqlite3 $(DB) $(foreach DT,$(DROPTABLES),"DROP VIEW IF EXISTS '$(DT)_view'" "DROP TABLE IF EXISTS '$(DT)_conflict'" "DROP TABLE IF EXISTS '$(DT)'")
	$(NANOBOT) init
	python3 $< $(DB)


### IEDB Data

# Tax ID -> epitope count
build/counts.tsv: src/get-counts.sql | build
	$(MIRROR_QUERY) < $< > $@

build/ncbi_include.tsv: src/ncbi-active-taxa.sql | build
	$(MIRROR_QUERY) < $< > $@

build/counts_full.tsv: src/combine_taxids.py build/counts.tsv build/ncbi_include.tsv
	python3 $^ $@

# Custom IEDB taxa
build/iedb_taxa.tsv: src/get-iedb-taxa.sql | build
	$(MIRROR_QUERY) < $< > $@


### NCBI Taxonomy

build/taxdmp.zip: | build
	# curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2023-10-01.zip

build/ncbitaxon.built: build/taxdmp.zip | $(DB)
	sqlite3 $(DB) "DROP TABLE IF EXISTS ncbitaxon"
	python3 src/ncbitaxon2ldtab.py $< $(DB)
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_subject ON ncbitaxon(subject)"
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_predicate ON ncbitaxon(predicate)"
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_object ON ncbitaxon(object)"
	sqlite3 $(DB) "ANALYZE ncbitaxon"
	touch $@


### Old Organism Tree

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


### Organism Tree

build/organism_core.html: src/build_organism_core.py src/organism_core.tsv
	python3 $^ $@

build/organism-tree.tsv: build/ncbitaxon.built src/assign_species.py src/organism_core.tsv build/iedb_taxa.tsv build/counts_full.tsv
	python3 $(word 2, $^) $(DB) $(filter %.tsv, $^) $@

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
