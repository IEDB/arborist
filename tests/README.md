# Arborist protein-tree tests

Two tiers, split by the `e2e` pytest marker (registered in `pytest.ini`). The
default suite is **siloed and fully offline** — no external binaries, no MySQL,
no network — so it runs anywhere, asahidake included. The `e2e` tier runs on
arborist-dev where the aligner binaries and build data live.

## Offline suite (default — run this while iterating)

```
pytest            # from tests/  (or `make test` from the repo root)
```

`pytest.ini` sets `addopts = -m "not e2e"`, so bare `pytest` runs only the
offline suite. It finishes in seconds. Each test builds its own hermetic
`build/` tree in a tmp dir (see the `silo` fixture in `conftest.py`) and
exercises one unit of logic:

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

## End-to-end tier (arborist-dev only)

```
pytest -m e2e     # from tests/  (or `make test-e2e` from the repo root)
```

Needs the aligner binaries in `../bin` and real build data. `-m e2e` overrides
the default deselect. The `e2e` tests also **self-skip** when their
prerequisites are absent, so a stray selection off-box is safe.

- `test_assignments.py` — generates a real `assign.py -b build -n 4` run against
  the committed golden fixtures (a session fixture, no Makefile wrapper), then
  diffs each species against its `expected-assignments.tsv`.
- `test_protein_output.py` — validates the protein-tree output files delivered
  to the IEDB backend (`build/proteins/latest/*`); driven by the `make leidos`
  target after a full pipeline build.

## Single-species real run (arborist-dev, inputs already present)

Iterate on one species without the whole pipeline:

```
python3 src/protein/protein_tree/assign.py -t 1642 -b build
```

`-t/--taxon_id` restricts to one species; `-b/--build_path` points at the
existing build tree.
