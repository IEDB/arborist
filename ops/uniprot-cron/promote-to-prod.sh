#!/usr/bin/env bash
# Push a refreshed arborist-dev tree to the two hosts that consume it.
#
#   sudo ops/uniprot-cron/promote-to-prod.sh --dry-run
#   sudo ops/uniprot-cron/promote-to-prod.sh --confirm
#
#   arborist.lji.org          <- all of build/species/, a full mirror
#   pepmatch-curation.lji.org <- every .pepidx, flattened into one directory
#
# Refuses unless the gate passed or Dan approved. This is the only script in the
# pipeline that writes outside arborist-dev.
set -euo pipefail

BUILD_HOST_FQDN=arborist-dev.lji.org
REPO=/mnt/data/arborist
SPECIES_DIR="$REPO/build/species"
STATE_DIR=/var/lib/arborist-uniprot-watch

PROD_HOST=arborist
PROD_SPECIES=/mnt/data/arborist/build/species
PEPMATCH_HOST=dmarrama@pepmatch-curation
PEPMATCH_DIR=/media/data/pepmatch/proteomes

# root@arborist-dev already has a key to prod. pepmatch-curation is reached with
# the existing dmarrama key rather than minting a second one; root can read it.
PEPMATCH_KEY=/home/dmarrama/.ssh/pepmatch_curation_rsa

SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=15)
PEPMATCH_SSH=("${SSH_OPTS[@]}" -i "$PEPMATCH_KEY" -o IdentitiesOnly=yes)

pepmatch_ssh_cmd() { echo "ssh ${PEPMATCH_SSH[*]}"; }

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }
die() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] FATAL: $*" >&2; exit 1; }

confirm=0
dry_run=0
for arg in "$@"; do
  case "$arg" in
    --confirm) confirm=1 ;;
    --dry-run) dry_run=1 ;;
    *) die "unknown argument: $arg (expected --confirm or --dry-run)" ;;
  esac
done
(( confirm || dry_run )) || die "refusing: pass --dry-run to preview or --confirm to promote"

host=$(hostname)
[[ "$host" == "$BUILD_HOST_FQDN" ]] || die "refusing: not $BUILD_HOST_FQDN (this is $host)"
[[ "$(id -u)" == "0" ]] || die "refusing: run as root (sudo)"
[[ -d "$SPECIES_DIR" ]] || die "refusing: no $SPECIES_DIR to promote"

cd "$REPO"
release=$(python3 -c "import json;print(json.load(open('$STATE_DIR/state.json'))['last_seen_release'])" 2>/dev/null || echo unknown)
gate_exit=$(python3 -c "import json;print(json.load(open('$STATE_DIR/state.json')).get('last_gate_exit','none'))" 2>/dev/null || echo none)

# --- the gate. Numbers unchanged promotes itself; anything else needs Dan.
approval="$STATE_DIR/APPROVED-$release"
if [[ "$gate_exit" == "0" ]]; then
  log "gate: diff-metrics exited 0 for $release"
elif [[ -f "$approval" ]]; then
  log "gate: $approval present, proceeding on Dan's approval"
else
  die "refusing: gate exit is '$gate_exit' and no $approval. Review build/reports/diff-metrics.tsv, then: mv $STATE_DIR/HOLD-$release $approval"
fi

# --- reachability, before anything expensive or destructive.
[[ -r "$PEPMATCH_KEY" ]] || die "missing $PEPMATCH_KEY (the pepmatch-curation key)"
ssh "${SSH_OPTS[@]}" "$PROD_HOST" true || die "cannot ssh $PROD_HOST as root"
ssh "${PEPMATCH_SSH[@]}" "$PEPMATCH_HOST" true \
  || die "cannot ssh $PEPMATCH_HOST with $PEPMATCH_KEY"
log "both destinations reachable"

# --- k=3 indexes. Arborist only builds k=5; pepmatch-curation serves 3-mers for
# a short list of heavily queried species.
if (( dry_run )); then
  log "dry-run: skipping 3-mer build"
else
  "$REPO/_venv/bin/python" src/protein/protein_tree/preprocess_3mers.py -b build/ \
    || die "3-mer preprocessing failed"
fi

rsync_common=(-aH --info=progress2 --partial-dir=.rsync-partial)
(( dry_run )) && rsync_common+=(-n)

# --- prod: hardlink snapshot first. Costs directory entries, not 83G, and rsync
# replaces files by rename so the snapshot keeps the old inodes intact.
snapshot="$PROD_SPECIES.prev-$release"
if (( dry_run )); then
  log "dry-run: would snapshot $PROD_HOST:$snapshot"
else
  log "snapshotting $PROD_HOST:$snapshot"
  ssh "${SSH_OPTS[@]}" "$PROD_HOST" \
    "rm -rf '$snapshot' && cp -al '$PROD_SPECIES' '$snapshot'" \
    || die "could not snapshot prod; nothing transferred"
fi

log "mirroring build/species/ -> $PROD_HOST:$PROD_SPECIES"
rsync "${rsync_common[@]}" --delete "$SPECIES_DIR/" "$PROD_HOST:$PROD_SPECIES/" \
  || die "prod rsync failed; snapshot $snapshot retained for rollback"

# --- pepmatch-curation: only the indexes, flattened to <taxon>_<k>mers.pepidx.
# No --partial: a reader mmapping a half-written index is the failure this avoids.
staging=$(mktemp -d /tmp/pepidx-promote.XXXXXX)
trap 'rm -rf "$staging"' EXIT
count=0
for index in "$SPECIES_DIR"/*/proteome_*mers.pepidx; do
  [[ -e "$index" ]] || continue
  taxon=$(basename "$(dirname "$index")")
  kmer=$(basename "$index" .pepidx); kmer=${kmer#proteome_}
  ln "$index" "$staging/${taxon}_${kmer}.pepidx" 2>/dev/null \
    || cp "$index" "$staging/${taxon}_${kmer}.pepidx"
  count=$(( count + 1 ))
done
log "staged $count index files"

log "syncing indexes -> $PEPMATCH_HOST:$PEPMATCH_DIR"
pepmatch_rsync=(-a --info=progress2 --temp-dir="$PEPMATCH_DIR/.staging")
(( dry_run )) && pepmatch_rsync+=(-n)
(( dry_run )) || ssh "${PEPMATCH_SSH[@]}" "$PEPMATCH_HOST" "mkdir -p '$PEPMATCH_DIR/.staging'"
rsync "${pepmatch_rsync[@]}" -e "$(pepmatch_ssh_cmd)" \
  "$staging/" "$PEPMATCH_HOST:$PEPMATCH_DIR/" \
  || die "pepmatch rsync failed; prod is already updated, snapshot $snapshot retained"

# --- verify. rsync -c re-reads and compares content on both sides. Only
# transfers and deletions count; directory mtimes do not.
differs() {
  rsync -acn --itemize-changes "$@" | grep -Eq '^(>|<|c|\*)'
}

if (( dry_run )); then
  log "dry-run: skipping checksum verification"
else
  log "verifying prod by checksum"
  differs --delete "$SPECIES_DIR/" "$PROD_HOST:$PROD_SPECIES/" \
    && die "prod differs after sync; snapshot $snapshot retained"
  log "verifying pepmatch by checksum"
  differs -e "$(pepmatch_ssh_cmd)" "$staging/" "$PEPMATCH_HOST:$PEPMATCH_DIR/" \
    && die "pepmatch differs after sync"

  log "dropping snapshot $snapshot"
  ssh "${SSH_OPTS[@]}" "$PROD_HOST" "rm -rf '$snapshot'"

  python3 - <<PY
import json, pathlib, datetime
path = pathlib.Path('$STATE_DIR/state.json')
if path.exists():
  state = json.loads(path.read_text())
  state['last_sync_completed'] = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
  state['last_sync_release'] = '$release'
  tmp = path.with_suffix('.json.tmp')
  tmp.write_text(json.dumps(state, indent=2, sort_keys=True) + '\n')
  tmp.replace(path)
PY
fi

(( dry_run )) && { log "dry-run complete. Nothing was written."; exit 0; }
log "promoted $release: prod mirrored, $count indexes on pepmatch-curation."
