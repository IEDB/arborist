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

See [docs](./docs/) for more detailed information and instructions.

## Usage

The `Makefile` defines and documents all the specific steps for Arborist.
Run `make help` to see the list of main tasks.

You can either run `make` directly or inside a Docker container.
For Docker, run `./run_image.sh make` or `sudo -E ./run_image.sh make`.
If you aren't using Docker,
first install the required software by running `make deps`.
NOTE: Arborist currently supports only Linux on the x86_64 architecture.

The suggested workflow is:

1. Update the cache with the latest IEDB tables
  by running `src/iedb/update-cache`.
  This requires MySQL/MariaDB connection parameters to be set
  as IEDB_MYSQL_* environment variables:
  IEDB_MYSQL_HOST,
  IEDB_MYSQL_PORT,
  IEDB_MYSQL_USER,
  IEDB_MYSQL_PASSWORD,
  IEDB_MYSQL_DATABASE.
2. Run `make all` to build all trees.
3. Run `make serve` to start the web interface on <http://localhost:3000>.

These are the key Make tasks for building trees,
in their dependency order:

- `make iedb` load IEDB data:
  This runs the `src/iedb/update-cache` script
- `make ncbitaxon` build the NCBI Taxonomy
- `make organism` build the organism and subspecies trees:
  This also creates the list of "active species" used by IEDB,
  and the "active taxa" that fall under these species.
- `make proteome` select a proteome for each active species
- `make protein` build the protein tree
- `make all` build all trees

TODO: build more trees: peptide, molecule, assay, disease, geolocation, ...

Here are some other important Make tasks:

- `make deps` install required software
- `make serve` run the web interface on <http://localhost:3000>
- `make clean` remove all build files
- `make clobber` remove all generated files
- `make help` print this message


## Files and Directories

- `bin/` contains any required binaries that aren't already installed
- `build/` all sorts of generated files
  - `iedb/` selected tables from IEDB for use here
  - `arborist/` general build files
  - `<species_id>/` species-specific build files
- `cache/` compressed data from various sources
  - `iedb/` selected tables from IEDB
  - `ncbitaxon/` NCBI Taxonomy's `taxdmp.zip` files
- `current/` links to the cached data to use for builds
  - `iedb` links to a subdirectory of `cache/iedb/`
  - `taxdmp.zip` links to a file in `cache/ncbitaxon/`
- `result/` TODO date-stamped directories of results, and `latest` link
- `src/`
  - `iedb/` config and schemas for IEDB data
  - `arborist/` config and schemas for Arborist tables
  - `species/` config and schemas for species proteomes and protein trees
  - `organism/` scripts for building the organism tree
  - `proteome/` scripts for selecting proteomes
  - `util/` utility scripts for working with databases
  - `templates/` Nanobot HTML templates
