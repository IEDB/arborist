#!/bin/sh
# Weekly Arborist run: ./arborist.sh DATESTAMP MYSQL_HOST
#
# Any failure below -- a failed pre-check, a make error, a make child killed by the
# kernel (make prints "Killed" and exits 2: the OOM kill of 2026-10-18), or a
# missing/stale output after a "successful" make -- sends an urgent email via protein_tree.notify_failure
# (default to dmarrama@lji.org only; ARBORIST_FAILURE_EMAIL overrides) with the
# tail of the run log. The exit status is the build's own. This script cannot
# report its own SIGKILL.
echo "--------------------------------------"

LOG="${ARBORIST_RUN_LOG:-$(pwd)/arborist-last-run.log}"
: > "$LOG"
STARTED=$(mktemp)

# alert STAGE EXIT_CODE [DETAIL]. Never stops the run; a send failure is logged.
alert() {
  PYTHONPATH="src/protein${PYTHONPATH:+:$PYTHONPATH}" "${VENV_PYTHON:-python3}" \
    -m protein_tree.notify_failure --stage "$1" --exit-code "$2" --log "$LOG" --detail "${3:-}" \
    || echo "ARBORIST: the failure email for '$1' could not be sent" | tee -a "$LOG"
}

# run_stage STAGE COMMAND...: run, copy all output to $LOG, alert on non-zero.
# The exit status goes through a file because POSIX sh has no pipefail.
run_stage() {
  stage=$1
  shift
  rcfile=$(mktemp)
  { "$@"; echo $? > "$rcfile"; } 2>&1 | tee -a "$LOG"
  rc=$(cat "$rcfile")
  rm -f "$rcfile"
  if [ "$rc" -ne 0 ]; then
    alert "$stage" "$rc"
  fi
  return "$rc"
}

echo "Activating environment..."
if [ ! -f _venv/bin/activate ]; then
  echo "Missing _venv, please setup_env.sh" | tee -a "$LOG"
  alert "activate _venv" 1 "Missing _venv, please setup_env.sh"
  exit 1
fi
. _venv/bin/activate
export VENV_PYTHON=$(pwd)/_venv/bin/python

DATESTAMP=$1
export MYSQL_HOST=$2

export IEDB_MYSQL_HOST="$MYSQL_HOST"
export IEDB_MYSQL_PORT="3306"
export IEDB_MYSQL_USER="iedb_query"
export IEDB_MYSQL_PASSWORD=$(cat .iedb_query)
export IEDB_MYSQL_DATABASE="iedb_query_$DATESTAMP"

export IEDB_MYSQL_QUERY="mysql --host=${IEDB_MYSQL_HOST} --port=${IEDB_MYSQL_PORT} --user=${IEDB_MYSQL_USER} --password=${IEDB_MYSQL_PASSWORD} ${IEDB_MYSQL_DATABASE}"

echo "Checking database connection..."
if ! run_stage "database connection check" $IEDB_MYSQL_QUERY -e "SELECT 1"; then
  echo "Could not connect to database ${IEDB_MYSQL_DATABASE} on host ${IEDB_MYSQL_HOST}:${IEDB_MYSQL_PORT}" | tee -a "$LOG"
  exit 1
fi

echo "--------------------------------------"
echo "Running Arborist..."

# A failed clean alerts but, as before, does not stop the build.
run_stage "make weekly_clean" make weekly_clean

run_stage "make weekly" make weekly
status=$?

if [ "$status" -eq 0 ]; then
  missing=""
  for output in disease-tree.owl molecule-tree.owl molecule-tree-with-gene.owl \
                parent-proteins.tsv source-parents.tsv epitope-mappings.tsv; do
    path="build/proteins/latest/$output"
    if [ ! -s "$path" ] || [ ! "$path" -nt "$STARTED" ]; then
      missing="$missing $path"
    fi
  done
  if [ -n "$missing" ]; then
    echo "make weekly exited 0 but outputs are missing, empty or older than this run:$missing" | tee -a "$LOG"
    alert "output check after make weekly" 1 "Missing, empty or stale:$missing"
    status=1
  fi
fi

rm -f "$STARTED"
exit "$status"
