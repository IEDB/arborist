#!/bin/sh
echo "--------------------------------------"

echo "Activating environment..."
. _venv/bin/activate || (echo "Missing _venv, please setup_env.sh"; exit 1)

read -p "Enter MySQL Username: " IEDB_MYSQL_USER
read -s -p "Enter MySQL Password: " IEDB_MYSQL_PASSWORD
echo ""

DATESTAMP=$1
export IEDB_MYSQL_HOST="iedb-mysql.liai.org"
export IEDB_MYSQL_PORT="33306"
export IEDB_MYSQL_DATABASE="iedb_query"
export IEDB_MYSQL_USER="$IEDB_MYSQL_USER"
export IEDB_MYSQL_PASSWORD="$IEDB_MYSQL_PASSWORD"

echo "--------------------------------------"
echo "Running Arborist..."
sudo make weekly