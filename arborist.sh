#!/bin/sh
echo "--------------------------------------"

echo "Activating environment..."
. _venv/bin/activate || (echo "Missing _venv, please setup_env.sh"; exit 1)
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
$IEDB_MYSQL_QUERY -e "SELECT 1" || (echo "Could not connect to database ${IEDB_MYSQL_DATABASE} on host ${IEDB_MYSQL_HOST}:${IEDB_MYSQL_PORT}"; exit 1)

echo "--------------------------------------"
echo "Running Arborist..."

make weekly_clean
make weekly
