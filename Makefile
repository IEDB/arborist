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

# all: build/organism-tree.owl build/subspecies-tree.owl build/active-species.tsv
all: build/organism-tree.owl build/organism-tree-old.built

clean:
	rm -rf build


### Set Up

build:
	mkdir -p $@

ROBOT := java -jar bin/robot.jar --prefix "iedb-taxon: https://ontology.iedb.org/taxon/" --prefix "ONTIE: https://ontology.iedb.org/ONTIE_"
NANOBOT := /home/knocean/nanobot.rs/target/release/nanobot
EXPORT := build/export.py
LDTAB := java -jar bin/ldtab.jar
DB := build/nanobot.db

build/export.py: | build/
	curl -L -o $@ "https://github.com/ontodev/valve.rs/raw/main/scripts/export.py"


### IEDB Data

# Tax ID -> epitope count
build/counts.tsv: src/get-counts.sql | build
	$(MIRROR_QUERY) < $< > $@

build/ncbi_include.tsv: src/ncbi-active-taxa.sql | build
	$(MIRROR_QUERY) < $< > $@

# build/counts_full.tsv: src/combine_taxids.py build/counts.tsv build/ncbi_include.tsv
# 	python3 $^ $@

# Custom IEDB taxa
# build/iedb_taxa.tsv: src/get-iedb-taxa.sql | build
# 	$(MIRROR_QUERY) < $< > $@


### Nanobot Database

$(DB): | build
	rm -f $@ build/*.built
	$(NANOBOT) init

.PHONY: save
save: $(EXPORT) $(DB)
	python3 $(EXPORT) data $(DB) src/schema/ table column datatype
	python3 $(EXPORT) data $(DB) src/ organism_core
	python3 $(EXPORT) data $(DB) build/ $$(grep build src/schema/table.tsv | cut -f1 | tr '\n' ' ')
	python3 src/sort_organism_core.py src/organism_core.tsv

.PHONY: reload
reload:
	sqlite3 $(DB) "DROP TABLE IF EXISTS 'column_conflict'"
	sqlite3 $(DB) "DROP TABLE IF EXISTS 'column'"
	sqlite3 $(DB) "DROP TABLE IF EXISTS 'datatype_conflict'"
	sqlite3 $(DB) "DROP TABLE IF EXISTS 'datatype'"
	$(NANOBOT) init


### NCBI Taxonomy

build/taxdmp.zip: | build
	# curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	# curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2023-06-01.zip
	curl -L -o $@ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2023-07-01.zip

build/ncbitaxon.built: build/taxdmp.zip | $(DB)
	sqlite3 $(DB) "DROP TABLE IF EXISTS ncbitaxon"
	python src/ncbitaxon2ldtab.py $< $(DB)
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_subject ON ncbitaxon(subject)"
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_predicate ON ncbitaxon(predicate)"
	sqlite3 $(DB) "CREATE INDEX idx_ncbitaxon_object ON ncbitaxon(object)"
	sqlite3 $(DB) "ANALYZE"
	touch $@


### Old Organism Tree

# Load an existing organism-tree.owl
build/organism-tree-old.built: build/ncbitaxon.built build/organism-tree-old.owl
	sqlite3 $(DB) "DROP TABLE IF EXISTS organism_tree_old"
	sed s/statement/organism_tree_old/g src/statement.sql | sqlite3 $(DB)
	$(LDTAB) import $(DB) $(word 2,$^) -t organism_tree_old
	# sqlite3 $(DB) "ANALYZE"
	touch $@


### Organism Tree

build/organism_core.html: src/build_organism_core.py src/organism_core.tsv
	python3 $^ $@

build/organism-tree.tsv: build/ncbitaxon.built src/assign_species.py src/organism_core.tsv build/iedb_taxa.tsv build/counts_full.tsv
	python $(word 2, $^) $(DB) $(filter %.tsv, $^) $@

# Build a new organism tree
build/organism-tree.built: build/ncbitaxon.built src/build_organism_tree.py build/organism-tree.tsv build/organism_core.html
	sqlite3 $(DB) "DROP TABLE IF EXISTS organism_tree"
	python $(word 2, $^) $(DB) $(filter %.tsv, $^)
	# sqlite3 $(DB) "ANALYZE"
	touch $@

build/organism-tree.ttl: build/organism-tree.built
	rm -f $@
	$(LDTAB) export $(DB) --table organism_tree $@



# .PHONY: load
# load: # src/schema/ build/organism-tree.owl build/subspecies-tree.owl build/protein-tree.owl
# 	rm -rf $(NANOBOTDB)
# 	$(NANOBOT) init
# 	sqlite3 $(NANOBOTDB) < src/prefixes.sql
# 	sed s/statement/organism_tree/g src/statement.sql | sqlite3 $(NANOBOTDB)
# 	$(LDTAB) import $(NANOBOTDB) build/organism-tree.owl -t organism_tree
# 	sed s/statement/subspecies_tree/g src/statement.sql | sqlite3 $(NANOBOTDB)
# 	$(LDTAB) import $(NANOBOTDB) build/subspecies-tree.owl -t subspecies_tree
# 	sed s/statement/protein_tree/g src/statement.sql | sqlite3 $(NANOBOTDB)
# 	$(LDTAB) import $(NANOBOTDB) build/protein-tree.owl -t protein_tree
# 	sqlite3 $(NANOBOTDB) ANALYZE

# build/organism-tree.tsv: build/ldtab.db
# 	rm -rf $@
# 	$(LDTAB) export $< $@

# build/organism-tree.sorted.tsv: build/organism-tree.tsv
# 	sort $< > $@


### Trees
# Subset of NCBITaxonomy terms rebuilt to be easier for immunologists to browse

#build/organism-tree.db: src/prefixes.sql src/build_tree.py build/ncbitaxon.db build/counts_full.tsv build/iedb_taxa.tsv build/ncbi_taxa.tsv build/taxon_parents.tsv build/top_level.tsv
#	rm -rf $@
#	sqlite3 $@ < $<
#	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# build/subspecies-tree.db: src/prefixes.sql src/add_all_taxa.py build/organism-tree.db build/ncbitaxon.db build/ncbi_include.tsv
# 	rm -rf $@
# 	sqlite3 $@ < $<
# 	python3 $(filter-out src/prefixes.sql,$^) $@ || (rm -rf $@ && exit 1)

# build/%-tree.ttl: src/db2ttl.py build/%-tree.db
# 	python3 $^ > $@

# remove other organism, unclassifed, and bad SARS coronavirus node
build/%-tree.owl: build/%-tree.ttl src/predicates.ttl
	$(ROBOT) merge --input $< --input $(word 2,$^) \
	annotate \
	--ontology-iri https://ontology.iedb.org/ontology/$(notdir $@) \
	--output $@


### Other Products

build/active-species.tsv: src/get_active_species.py build/organism-tree.db build/counts.tsv
	python3 $^ $@

