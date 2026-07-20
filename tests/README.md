# Protein-tree fast local test loop

Kill the `edit -> push -> git pull on dev -> full make -> wait` loop. Two ways
to iterate on a single species in seconds.

## 1. Offline regression (any box, no data/network required)

Reproduces the taxon 1642 crash (bus brief #124) against a frozen fixture in
`tests/fixtures/species_1642/`. No MySQL, no organism rebuild, no UniProt.

```
cd tests && make test-1642
```

(or `python3 -m pytest -q test_species_1642.py`)

This is the seed of the protein-tree carve-out: a single-species snapshot +
`assign.py` orchestration, exercised in isolation. To freeze another species,
copy the `species_1642` fixture layout: a one-row `active-species.tsv`, the
matching `iedb/peptide.tsv` + `peptide_source.tsv` rows, and a
`species/<id>/` directory holding whatever state you want to reproduce.

## 2. Single-species real run (on arborist-dev, inputs already present)

When the build inputs already exist on the dev box, iterate on one species
directly instead of running the whole pipeline:

```
python3 src/protein/protein_tree/assign.py -t 1642 -b build
```

`-t/--taxon_id` restricts the run to one species; `-b/--build_path` points at
the existing build tree. Seconds, not a full `make protein`.
