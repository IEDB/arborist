# Arborist protein-tree tests

Two tiers. The default suite is **siloed and fully offline** — no external
binaries, no MySQL, no network — so it runs anywhere, asahidake included. The
heavy end-to-end suite runs on arborist-dev where the aligner binaries and
build data live.

## Offline suite (default — run this while iterating)

```
cd tests && make test        # == python3 -m pytest -q
```

Runs in seconds. Each test builds its own hermetic `build/` tree in a tmp dir
(see the `silo` fixture in `conftest.py`) and exercises one unit of logic:

| file | subsystem | how it stays offline |
|------|-----------|----------------------|
| `test_species_1642.py`    | assign orchestration (on-demand select / skip) | `ProteomeSelector` spied |
| `test_source_assignment.py` | `SourceProcessor` (source antigen -> protein) | frozen blast `alignments.csv` |
| `test_peptide_assignment.py`| `PeptideProcessor` (epitope -> parent protein) | pepmatch (pip lib, no bin) |
| `test_proteome_tsv.py`      | `generate_proteome_tsv` FASTA parsing | pure Biopython/regex |
| `test_proteome_report.py`, `test_stats.py`, `test_alert.py`, `test_email.py` | reporting / alerting | pure polars |

The external aligners (blastp/mmseqs2) and ARC never run here: their outputs
are frozen as tiny fixtures, or replaced by pepmatch, which does exact peptide
matching in-process.

The bin/data-dependent tests (`test_assignments`, `test_leidos`) **self-skip**
when their prerequisites are absent, so `make test` is green on asahidake.

## End-to-end suite (arborist-dev only)

Needs the aligner binaries in `../bin` and real build data.

```
cd tests && make test-e2e
```

Generates a real `assign.py -b build -n 4` run, then asserts the golden
`expected-assignments.tsv` per species. `test_leidos.py` additionally validates
`build/proteins/latest/*` after a full pipeline build.

## Single-species real run (arborist-dev, inputs already present)

Iterate on one species without the whole pipeline:

```
python3 src/protein/protein_tree/assign.py -t 1642 -b build
```

`-t/--taxon_id` restricts to one species; `-b/--build_path` points at the
existing build tree.
