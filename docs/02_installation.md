# Installation

Linux **x86_64** only. Run all commands from the **repository root** unless noted.

## 1. System packages

These must be on the host `PATH` before Make can run. `make deps` does **not** install them.

| Package | Why |
|---------|-----|
| `make` | Build driver |
| Java JRE | `robot.jar`, `ldtab.jar` |
| Python 3.8+ | Build scripts; create `_venv` |
| `curl` | Downloads |
| `sqlite3` | **Required at Makefile parse time** — without it, even `make help` fails |
| MySQL/MariaDB client (`mysql`) | **Required at Makefile parse time**; used to pull IEDB |
| `unzip` | Extract QSV during `make deps` |

Debian/Ubuntu example:

```bash
sudo apt install -y make default-jre-headless python3 curl sqlite3 mariadb-client unzip
```

Confirm:

```bash
make help    # must print tasks, not "Please install SQLite 3"
```

Also useful: `git` (for some pip VCS deps).

## 2. Python environment

Install from the **root** `requirements.txt` (includes the editable protein package and pinned `pepmatch`).

```bash
python3 -m venv _venv
export VENV_PYTHON=$PWD/_venv/bin/python
"$VENV_PYTHON" -m pip install -U pip
"$VENV_PYTHON" -m pip install -r requirements.txt
```

With [uv](https://github.com/astral-sh/uv):

```bash
uv venv _venv
export VENV_PYTHON=$PWD/_venv/bin/python
uv pip install --python "$VENV_PYTHON" -r requirements.txt
```

`setup_env.sh` is an older helper that creates `_venv` via `pip`; prefer the commands above if you do not want its `sudo chown` step.

**Important:** the Makefile invokes `$(VENV_PYTHON)` for Python steps. Export it in every shell before `make`, or activate the venv (`source _venv/bin/activate`) so the same interpreter is used consistently.

```bash
export VENV_PYTHON=$PWD/_venv/bin/python
```

## 3. Project tools (`make deps`)

Downloads/builds ROBOT, LDTab, Nanobot, valve-export, QSV, BLAST+, HMMER (`hmmscan`), and MMseqs2 into `bin/`. The Makefile prepends `bin/` to `PATH` for its recipes.

```bash
export VENV_PYTHON=$PWD/_venv/bin/python
make deps
```

Sanity check (optional):

```bash
ls bin/robot bin/nanobot bin/blastp bin/hmmscan bin/mmseqs bin/qsv
```

## 4. IEDB MySQL (needed for `make iedb` and full `make all`)

Set all five (names are fixed):

```bash
export IEDB_MYSQL_HOST='<host>'
export IEDB_MYSQL_PORT='<port>'
export IEDB_MYSQL_USER='<user>'
export IEDB_MYSQL_PASSWORD='<password>'
export IEDB_MYSQL_DATABASE='<database>'
```

- Prefer a private env file or secret store; do not commit credentials.
- Some deployments use a mode-`0600` password file (e.g. `.iedb_query`) with wrapper scripts such as `arborist.sh` — that is optional; the Make/`mysql2tsv` path is the env vars above.
- The build host must be able to open a TCP connection to that MySQL instance (on-network, VPN, or firewall allowlist as your site requires). Verify with a simple client check that does not log the password.

`make ncbitaxon` and `make deps` do **not** need these variables.

## 5. External ontology input

`make molecule` (and therefore `make all`) expects this file at the **repo root**:

```text
nonpeptide-tree-20240305.owl
```

It is gitignored. Copy it from your team’s Arborist data host or artifact store. Without it, earlier stages can still run; molecule/all will fail when they reach that step.

## 6. Runtime knobs

| Variable | Purpose |
|----------|---------|
| `VENV_PYTHON` | Python used by Make (set this) |
| `ALERT_EMAIL=0` | Disable SMTP alert side effects (recommended on non-prod / personal hosts; default may try an institutional relay) |
| `IEDB_MYSQL_*` | IEDB pull |

Example:

```bash
export VENV_PYTHON=$PWD/_venv/bin/python
export ALERT_EMAIL=0
```

## 7. Disk and network

- Fresh full builds pull NCBI taxdump, IEDB extracts, many UniProt proteomes, ontologies, and large assignment intermediates. Keep **ample free disk** (on the order of hundreds of GiB for a comfortable full E2E).
- Public HTTPS/FTP egress is required for NCBI, UniProt, GitHub releases, OBO, etc.
- Optional speedup: copy a coherent `cache/`, `current/`, and `build/species/` snapshot from an existing build host (preserve symlink targets). That is not a substitute for documenting a from-scratch path.

## 8. Docker note

`Dockerfile` / `run_image.sh` exist but are not a complete turnkey E2E by themselves (secrets, external OWL, and full source layout still matter). Prefer native `_venv` + `make deps` unless you maintain the image path yourself.
