# Arborist

Arborist builds trees for the IEDB.
The trees are used for the user interface on <https://iedb.org>
and the IEDB curation interface,
and also for validating IEDB data.
They combine data from the IEDB
with community ontologies such as the NCBI Taxonomy
and open scientific databases such as UniProt and Genbank.

WARN: This version of Arborist is still work-in-progress.
It makes extensive use of [Nanobot](https://github.com/ontodev/nanobot.rs),
which is also work-in-progress.

See [docs](./docs/) for detailed guides.
**New machine?** Start with [Installation](./docs/02_installation.md) then [Quick Start](./docs/03_quick_start.md).

## Requirements

- Linux **x86_64** only
- Host packages: `make`, Java (JRE), Python 3.8+, `curl`, `sqlite3`, MySQL/MariaDB client, `unzip`
- Network for public downloads (NCBI, UniProt, GitHub, …) and, for IEDB pulls, reachability to the IEDB MySQL host
- Disk: plan on **hundreds of GiB** free for a full fresh build (caches, per-species data, assignments)

## Fresh host (short path)

From the repo root:

```bash
# 1. Host packages (Debian/Ubuntu example)
sudo apt install -y make default-jre-headless python3 curl sqlite3 mariadb-client unzip

# 2. Python env (root requirements.txt, not only src/protein)
python3 -m venv _venv
# or: uv venv _venv
export VENV_PYTHON=$PWD/_venv/bin/python
$VENV_PYTHON -m pip install -U pip
$VENV_PYTHON -m pip install -r requirements.txt
# uv alternative: uv pip install --python "$VENV_PYTHON" -r requirements.txt

# 3. Project tools into bin/
make deps

# 4. Secret-free smoke (no IEDB credentials)
mkdir -p cache
ALERT_EMAIL=0 make ncbitaxon

# 5. Full build (needs IEDB_* env + nonpeptide OWL — see docs)
# export IEDB_MYSQL_HOST=... IEDB_MYSQL_PORT=... IEDB_MYSQL_USER=...
# export IEDB_MYSQL_PASSWORD=... IEDB_MYSQL_DATABASE=...
# place nonpeptide-tree-20240305.owl in the repo root
ALERT_EMAIL=0 make all
make serve   # http://localhost:3000
```

Always export `VENV_PYTHON` (or activate `_venv`) before `make`.
Set `ALERT_EMAIL=0` on machines that should not send SMTP alerts through the default LJI relay.

## Usage

The `Makefile` defines the build. Run `make help` for tasks.

You can run `make` on the host or in Docker (`./run_image.sh make`).
Native host + `_venv` is the better-tested path for a full E2E.

Suggested workflow once installed:

1. Configure IEDB MySQL via `IEDB_MYSQL_*` (see installation docs).
2. Place `nonpeptide-tree-20240305.owl` in the repo root (required for `make molecule` / `make all`).
3. `ALERT_EMAIL=0 make all`
4. `make serve` → <http://localhost:3000>

### Main build targets (dependency order)

| Target | Role |
|--------|------|
| `make deps` | Download/build tools into `bin/` |
| `make iedb` | Pull IEDB tables into `cache/` + local DB |
| `make ncbitaxon` | NCBI taxonomy → local SQLite |
| `make organism` | Organism / subspecies trees, active species |
| `make proteome` | Select + download proteomes (slow; many UniProt calls) |
| `make protein` | Peptide/source assignment + protein tree |
| `make molecule` | Merge protein + non-peptide trees |
| `make disease` | Disease tree |
| `make leidos` | Publish layout + output E2E tests |
| `make all` | `deps` through `leidos` |

Other:

- `make weekly` — full-ish rebuild **without** `proteome` (reuses species data). Not a fresh-machine E2E.
- `make clean` / `make clobber` — remove build artifacts (`clobber` is the hard reset)
- `make help`

## Files and directories

- `bin/` — tools installed by `make deps` (not committed)
- `build/` — generated trees, DBs, per-species work
- `cache/` — compressed upstream snapshots (IEDB, NCBI taxdump, …)
- `current/` — symlinks to the active cache entries
- `_venv/` — local Python env (gitignored)
- `nonpeptide-tree-20240305.owl` — **required external input** for molecule (gitignored; obtain from your team’s data store)
- `src/` — build scripts and schemas
- `tests/` — unit tests; `pytest -m e2e` needs a completed build + aligner binaries
- `docs/` — installation, quick start, design notes

## Python dependencies

Pinned in root `requirements.txt` (includes editable `./src/protein` and `pepmatch`).
Use that file for a complete env. Keep `pepmatch` on the version pinned there; older builds can fail on mixed empty/hit assignment frames.
