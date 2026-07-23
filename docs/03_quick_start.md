# Quick Start

All commands from the **repository root**. Complete [Installation](./02_installation.md) first (`sqlite3`/`mysql`/`unzip`, `_venv`, `VENV_PYTHON`, `make deps`).

```bash
export VENV_PYTHON=$PWD/_venv/bin/python
export ALERT_EMAIL=0
mkdir -p cache
```

## 1. Secret-free smoke

No IEDB credentials. Validates downloads, Nanobot, SQLite, and the taxonomy conversion:

```bash
make ncbitaxon
```

Success: `build/arborist/ncbitaxon.built` exists and `current/taxdmp.zip` points at a cached taxdump.

## 2. IEDB + organism

Requires `IEDB_MYSQL_*` and network reachability to the database:

```bash
make iedb
make organism
```

Check for active-species / organism outputs under `build/` (see `make help` and organism docs).

## 3. Full build

Also requires `nonpeptide-tree-20240305.owl` in the repo root (see installation).

```bash
ALERT_EMAIL=0 make all
```

`make all` runs, in order: `deps`, `iedb`, `ncbitaxon`, `organism`, `proteome`, `protein`, `molecule`, `disease`, `leidos`.

Notes:

- **`proteome`** is slow (many UniProt requests). Be polite with parallel jobs; prefer a serial first full run.
- **`protein`** is large (multi‑GiB intermediates, high RAM on big species).
- **`leidos`** rearranges publication outputs and runs output E2E tests — use an isolated worktree if you share a machine with production data.
- First successful E2E means **through `leidos`**, not only “human peptides assigned.”

## 4. Web UI

```bash
make serve
```

Open <http://localhost:3000>.

## Other targets

| Command | When to use |
|---------|-------------|
| `make weekly` | Rebuild without re-selecting proteomes (reuses `build/species/`). **Not** a fresh-host E2E. |
| `make protein` | After organism/proteome data exist |
| `make clean` | Drop most build products; keeps species downloads when applicable |
| `make clobber` | Hard reset (`bin/`, `build/`, `cache/`, `current/`, …) |
| `make help` | List tasks |

## Offline / unit tests

```bash
"$VENV_PYTHON" -m pytest tests/ -q
```

End-to-end tests marked `e2e` need aligner binaries and a completed build; run only when that data is present:

```bash
"$VENV_PYTHON" -m pytest -m e2e
```
