# Quick Start

This guide provides the essential commands to build all the Arborist data and run the local web server. All commands should be run from the root of the `arborist` project directory.

## 1. Install Dependencies

Before running any build processes, you must install the required command-line tools and software. See the [installation](./02_installation.md) guide for all dependencies. The `Makefile` includes a command to automate this for most dependencies.
```bash
make deps
```

This command will check for tools like ROBOT, Nanobot, QSV, BLAST, and others. If they are not found in your `PATH`, it will download them into the local `bin/` directory.

## 2. Run the Full Build

The simplest way to generate all the trees is to use the `all` command. This will execute all the necessary steps in the correct order, from fetching IEDB data to building the final files.

Note: This is a long-running process that requires a connection to the IEDB MySQL database initially (see [installation](./02_installation.md) guide for database configuration).
```bash
make all
```

This single command is a convenient wrapper that runs the following steps sequentially:
* `make iedb`: Fetches the latest data from the IEDB MySQL database and populates a local cache.
* `make ncbitaxon`: Downloads the NCBI Taxonomy dump and converts it into a local SQLite database.
* `make organism`: Builds the main organism and subspecies trees.
* `make proteome`: Selects the best UniProt proteome for each active species and pulls all relevant metadata. This is a particularly time-consuming step.
* `make protein`: Builds the protein tree by assigning IEDB epitopes and antigens to proteins from the selected proteomes.
* `make molecule`: Builds the molecule tree by combining the protein tree with non-peptidic molecule data.
* `make disease`: Builds the disease tree from DOID and ONTIE.
* `make leidos`: Prepares and copies the final output files for the IEDB backend and runs a few tests.

## 3. Run the Web Interface

After the build is complete, you can explore the generated trees and data tables using the local web server.
```bash
make serve
```

This will start the Arborist web interface, which you can access in your browser at `http://localhost:3000`.

## Other Useful Commands

* `make weekly`: Runs the entire build process except for the `proteome` selection step. This is run each week on our production server to capture new data without having to re-evaluate the best proteome for each species.
* `make clean`: Removes most temporary build files but preserves the `build/species/` directory, which contains downloaded proteomes and metadata. This allows for a quicker rebuild.
* `make clobber`: A more aggressive cleanup command that removes all generated files and directories, including `bin/`, `build/`, `cache/`, and `current/`. This will force a complete fresh build on the next run.