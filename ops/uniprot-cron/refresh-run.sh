#!/usr/bin/env bash
# Full proteome refresh on arborist-dev, end to end, contained to dev.
#
#   sudo ops/uniprot-cron/refresh-run.sh [DATESTAMP]
#
# Wipes build/species/, re-selects every active species against current UniProt,
# and runs the full build through `protein` so arborist-dev is a usable snapshot
# and the assignment numbers exist for a later gate. Promotes nothing.
#
# `arborist.sh` cannot be used for this: it ends in `make weekly`, and `weekly`
# omits `proteome`, so it would reassign against stale FASTA and never call
# select_proteome.py.
set -euo pipefail

BUILD_HOST_FQDN=arborist-dev.lji.org
REPO=/mnt/data/arborist
LOCK=/var/lock/arborist-build.lock
STATE_DIR=/var/lib/arborist-uniprot-watch
MYSQL_HOST_DEFAULT=web02.internal.iedb.org
DATESTAMP_LOOKBACK=4        # Sundays to probe before giving up
MIN_FREE_GB=200

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }
die() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] FATAL: $*" >&2; exit 1; }

# --- containment. This script deletes a directory tree; it runs on dev or nowhere.
# Exact full hostname, no glob: a short-name match would also accept an
# arborist-dev in some other domain.
host=$(hostname)
[[ "$host" == "$BUILD_HOST_FQDN" ]] || die "refusing: not $BUILD_HOST_FQDN (this is $host)"
[[ "$(id -u)" == "0" ]] || die "refusing: run as root (sudo)"
[[ -d "$REPO/.git" && -f "$REPO/arborist.sh" ]] || die "refusing: $REPO is not the arborist checkout"

cd "$REPO"
MYSQL_HOST="${MYSQL_HOST:-$MYSQL_HOST_DEFAULT}"

# --- env preamble, same as arborist.sh lines 4-19.
# shellcheck source=/dev/null
. _venv/bin/activate || die "missing _venv, run setup_env.sh"
export VENV_PYTHON="$REPO/_venv/bin/python"
export IEDB_MYSQL_HOST="$MYSQL_HOST"
export IEDB_MYSQL_PORT=3306
export IEDB_MYSQL_USER=iedb_query
IEDB_MYSQL_PASSWORD="$(cat .iedb_query)"
export IEDB_MYSQL_PASSWORD

# --- DATESTAMP: most recent Sunday, walking back until that snapshot exists.
# The build runs monthly at best, so a missing snapshot must degrade to last
# week's data, not abort the whole run.
snapshot_exists() {
  mysql --host="$IEDB_MYSQL_HOST" --port="$IEDB_MYSQL_PORT" \
        --user="$IEDB_MYSQL_USER" --password="$IEDB_MYSQL_PASSWORD" \
        --batch --skip-column-names \
        -e "SHOW DATABASES LIKE 'iedb_query_$1'" 2>/dev/null | grep -q .
}

DATESTAMP="${1:-}"
if [[ -n "$DATESTAMP" ]]; then
  snapshot_exists "$DATESTAMP" || die "iedb_query_$DATESTAMP not found on $MYSQL_HOST"
else
  for (( week = 0; week < DATESTAMP_LOOKBACK; week++ )); do
    candidate=$(date -d "-$(( $(date +%w) + week * 7 )) days" +%Y%m%d)
    if snapshot_exists "$candidate"; then
      DATESTAMP="$candidate"
      (( week > 0 )) && log "WARNING: fell back $week week(s); iedb_query_$candidate is the newest snapshot"
      break
    fi
    log "iedb_query_$candidate absent, trying the Sunday before"
  done
  [[ -n "$DATESTAMP" ]] || die "no iedb_query_* snapshot in the last $DATESTAMP_LOOKBACK Sundays on $MYSQL_HOST"
fi
export IEDB_MYSQL_DATABASE="iedb_query_$DATESTAMP"
log "IEDB snapshot: $IEDB_MYSQL_DATABASE on $MYSQL_HOST"

export IEDB_MYSQL_QUERY="mysql --host=${IEDB_MYSQL_HOST} --port=${IEDB_MYSQL_PORT} --user=${IEDB_MYSQL_USER} --password=${IEDB_MYSQL_PASSWORD} ${IEDB_MYSQL_DATABASE}"
$IEDB_MYSQL_QUERY -e "SELECT 1" >/dev/null || die "cannot reach $IEDB_MYSQL_DATABASE on $MYSQL_HOST"

# --- disk. assign.py reuses any per-species output with size > 0, so a run that
# dies on a full disk leaves truncated files the NEXT run silently trusts.
free_gb=$(df -BG --output=avail /mnt/data | tail -1 | tr -dc '0-9')
(( free_gb >= MIN_FREE_GB )) || die "only ${free_gb}G free on /mnt/data, need ${MIN_FREE_GB}G"
log "disk: ${free_gb}G free on /mnt/data"

# --- the run itself, under the shared build lock so a weekly cannot overlap.
release=$(python3 -c "import json,sys; print(json.load(open('$STATE_DIR/state.json'))['last_seen_release'])" 2>/dev/null || echo unknown)
log "starting refresh for UniProt release $release"

# state.json tracks what UniProt is on; these fields track what we actually
# built, so "is dev a release behind?" is answerable without guessing.
set_state() {
  [[ -f "$STATE_DIR/state.json" ]] || return 0
  STATE_FILE="$STATE_DIR/state.json" STATE_KEY="$1" STATE_VALUE="$2" python3 - <<'PY'
import json, os, pathlib
path = pathlib.Path(os.environ['STATE_FILE'])
state = json.loads(path.read_text())
state[os.environ['STATE_KEY']] = os.environ['STATE_VALUE']
tmp = path.with_suffix('.json.tmp')
tmp.write_text(json.dumps(state, indent=2, sort_keys=True) + '\n')
tmp.replace(path)
PY
}

exec 9>"$LOCK"
flock -n 9 || die "build lock held (weekly running?); refresh skipped"

set_state last_refresh_started "$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Everything under build/species is derived: taxa.txt, epitopes.tsv, sources.tsv
# from the IEDB stages, and proteome.fasta/species-data.tsv/*.pepidx from
# select_proteome.py. A release moves all of it together, so delete rather than
# invalidate piecemeal. Literal path, no variable, on purpose.
log "wiping /mnt/data/arborist/build/species/"
rm -rf /mnt/data/arborist/build/species/

make weekly_clean
make deps iedb ncbitaxon organism proteome protein
log "build complete"

# --- stamp. .pepidx carries no release field, so this file is the only way to
# answer "which UniProt release built this tree?"
stamp="$REPO/build/reports/release-stamp.tsv"
mkdir -p "$(dirname "$stamp")"
{
  printf 'uniprot_release\tiedb_snapshot\tmysql_host\tbuilt_at\tgit_sha\n'
  printf '%s\t%s\t%s\t%s\t%s\n' "$release" "$IEDB_MYSQL_DATABASE" "$MYSQL_HOST" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(git -C "$REPO" rev-parse HEAD)"
} > "$stamp"
log "wrote $stamp"
set_state last_refresh_completed "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
set_state refreshed_release "$release"
log "state: dev now built against $release"

# --- the gate. Nothing to block yet (promotion does not exist), so this records
# the verdict and marks a HOLD rather than pretending to guard anything.
set +e
"$VENV_PYTHON" src/protein/protein_tree/diff_metrics.py -b build/
gate=$?
set -e
case "$gate" in
  0)  log "gate: numbers identical -- eligible for promotion once that exists" ;;
  10) log "gate: numbers changed within thresholds -- needs review" ;;
  20) log "gate: ESCALATED -- volume drop, status regression or proteome ID flip" ;;
  *)  log "gate: diff-metrics failed with exit $gate" ;;
esac
if (( gate != 0 )); then
  hold="$STATE_DIR/HOLD-$release"
  printf 'gate_exit=%s\nat=%s\ndiff=%s\n' "$gate" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$REPO/build/reports/diff-metrics.tsv" > "$hold"
  log "wrote $hold"
fi
set_state last_gate_exit "$gate"

log "refresh done. Nothing promoted; prod is untouched."
