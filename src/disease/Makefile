ROBOT := java -jar ../../bin/robot.jar --prefix "ONTIE: https://ontology.iedb.org/ontology/ONTIE_"

TODAY := $(shell date +%Y-%m-%d)
BASE = https://ontology.iedb.org/ontology

all: build/disease-tree.owl

clean:
	rm -rf build

build:
	mkdir $@

# ---------- Build the ONTIE subset ---------- 

build/ontie.owl: | build
	curl -Lk -o $@ https://raw.githubusercontent.com/IEDB/ONTIE/master/ontie.owl

build/ontie-disease.ttl: build/ontie.owl ontie-remove.txt
	$(ROBOT) filter \
	 --input $< \
	 --term ONTIE:0003543 \
	 --term DOID:4 \
	 --select "self descendants annotations" \
	 --trim false \
	 remove \
	 --term-file $(word 2,$^) \
	 --output $@

# ---------- Build the DOID subset ---------- 

# Get DOID terms used in IEDB data
# Plus the manually-added DOID terms

# Active DOIDs are found in the curation database
# but the Google sheet (build/diseases.tsv) should contain all of them

build/diseases.tsv: | build
	curl -L -o $@ 'https://docs.google.com/spreadsheets/d/16M0RyUBEw_fW09x0U2X1vIXrek_dl2qMvqcOA7xI0WU/export?format=tsv&gid=1130995017'

build/doid-terms.txt: build/diseases.tsv
	cat $< | cut -f1 > $@

build/doid.owl:
	$(ROBOT) merge \
	--input-iri http://purl.obolibrary.org/obo/doid.owl \
	reason \
	relax \
	reduce \
	remove \
	--select "equivalents" \
	--output $@

# get everything under 'neoplasm'
build/cancer-module.owl: build/doid.owl
	$(ROBOT) filter \
	--input $< \
	--term DOID:14566 \
	--select "self descendants annotations" \
	--output $@

# get the structure of external terms
build/doid-external.owl: build/doid.owl doid-terms.txt
	$(ROBOT) filter \
	--input $< \
	--term-file $(word 2,$^) \
	--select "self parents equivalents" \
	--trim false \
	--output $@

# get the IDs of all external terms needed
build/ext-doid-terms.tsv: build/doid-external.owl ext-doid-terms.rq
	$(ROBOT) query \
	--input $< \
	--query $(word 2,$^) $@
	tail -n +2 $@ | awk -F'"' '{print $$2}' > $@.tmp 
	mv $@.tmp $@

build/ext-doid-terms.txt: build/ext-doid-terms.tsv manual-ext.txt
	cat $^ > $@

build/doid-external-plus.owl: build/doid.owl build/ext-doid-terms.txt
	$(ROBOT) extract \
	--input $< \
	--method MIREOT \
	--lower-terms $(word 2,$^) \
	--output $@

build/doid-subset.owl: build/doid.owl build/doid-terms.txt build/doid-external-plus.owl \
build/doid-external.owl doid-remove.txt 
	$(ROBOT) extract \
	--input $< \
	--method MIREOT \
	--lower-terms $(word 2,$^) \
	merge \
	--input $(word 3,$^) \
	--input $(word 4,$^) \
	remove \
	--term-file $(word 5,$^) \
	remove \
	--term HP:0012823 \
	--select "self descendants" \
	--output $@

build/remove-doid-classes.txt: build/doid-subset.owl remove-doid-classes.rq
	$(ROBOT) query \
	--input $< \
	--query $(word 2,$^) $(basename $@).tsv
	tail -n +2 $(basename $@).tsv | awk -F'"' '{print $$2}' > $@
	rm $(basename $@).tsv

build/doid-module.owl: build/doid-subset.owl build/remove-doid-classes.txt doid-remove.txt
	$(ROBOT) remove \
	--input $< \
	--term-file $(word 2,$^) \
	remove \
	--term-file $(word 3,$^) \
	reduce \
	--output $@

build/disease-tree.owl: build/doid-module.owl build/ontie-disease.ttl build/cancer-module.owl update-doid.ru
	$(ROBOT) merge \
	 --input $< \
	 --input $(word 2,$^) \
	 --input $(word 3,$^) \
	query \
	 --update $(word 4,$^) \
	reason \
	reduce \
	annotate \
	 --ontology-iri $(BASE)/$(notdir $@) \
	 --version-iri $(BASE)/$(TODAY)/$(notdir $@) \
	 --output $@
