# arborist


## How it Works

- `bin/` contains any required binaries that aren't already installed
- `build/`
  - `iedb/`
  - `arborist/`
  - `<taxon_id>/` results for each active species
- `cache/` contains compressed data from upstream sources
  - `iedb/`
  - `ncbitaxon/`
- `current/` contains links to current cache to use
- `result/` date-stamped directories of results, and `latest` link
- `src/`
  - `iedb/`
  - `arborist/`
  - `species/`
  - `organism/`
  - `util/`
  - `templates/`

## Usage

See the `Makefile` for specific steps.

Optionally use Docker, or install dependencies to `bin/`:

```sh
$ make deps
```

Fetch source data from IEDB, which probably requires credentials:

```sh
$ src/cache-iedb iedb_query_20231029
```

Build IEDB tables

```sh
$ make iedb
```

Run the web interface:

```sh
$ make serve
```
