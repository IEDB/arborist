#!/bin/sh

mkdir -p build/ cache/ current/
docker build --tag arborist:latest .
docker run --rm -it \
  --publish 3000:3000 \
  --env IEDB_MYSQL_HOST \
  --env IEDB_MYSQL_PORT \
  --env IEDB_MYSQL_USER \
  --env IEDB_MYSQL_PASSWORD \
  --env IEDB_MYSQL_DATABASE \
  --mount type=bind,source=$(pwd)/build,target=/arborist/build \
  --mount type=bind,source=$(pwd)/cache,target=/arborist/cache \
  --mount type=bind,source=$(pwd)/current,target=/arborist/current \
  --mount type=bind,source=$(pwd)/src,target=/arborist/src \
  arborist:latest \
  $@
