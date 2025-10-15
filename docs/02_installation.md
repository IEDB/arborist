# Installation

To run the Arborist project, you need to set up an environment with several system-level dependencies, Python packages, and database connection details. This guide will walk you through the process.

## 1. Dependencies

Arborist relies on a number of command-line tools for data processing, ontology manipulation, and bioinformatics tasks. Please ensure the following tools are installed and available in your system's `PATH`.

### System-Level Binaries and Software
* **make**: For executing the build commands defined in Makefiles.
* **Java**: A Java Runtime Environment (JRE) is required to run `robot.jar`.
* **Python**: The Python interpreter needs to be installed using at least Python 3.7 or greater.
* **curl**: For downloading files from the web, used in the `disease` module Makefile.
* **mysql-client**: The `mysql` command-line client is needed by the `src/util/mysql2tsv` script to query the IEDB database.
* **sqlite3**: The `sqlite3` command-line tool is used for interacting with the local SQLite databases.

### Bioinformatics and General Tools
These tools are primarily used by the `protein` module for sequence alignment and analysis.
* **NCBI BLAST+**: The `blastp` and `makeblastdb` executables are required. Download them from the [NCBI BLAST website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
* **MMseqs2**: Used for fast, sensitive protein searching. Installation instructions are on the [MMseqs2 GitHub page](https://github.com/soedinglab/MMseqs2).
* **HMMER**: The `hmmscan` tool is used for protein domain analysis. You can get it from the [HMMER website](http://eddylab.org/software/hmmer/).
* **ROBOT**: A tool for working with ontologies, used for filtering and merging. Get it from the [ROBOT GitHub page](https://github.com/ontodev/robot).
* **LDTab**: A tool for converting OWL ontologies into linked data tables. Get it from the [LDTab GitHub page](https://github.com/ontodev/ldtab.clj).
* **Nanobot**: The backend tool for serving and browsing trees/tables. Get it from the [Nanobot GitHub page](https://github.com/ontodev/nanobot.rs).
* **valve-export**: A lightweight ontology validation engine written in Rust. Find it in the [VALVE GitHub repository](https://github.com/ontodev/valve.rs).
* **QSV**: A high-performance command-line toolkit for CSV/TSV files. Get it from the [QSV GitHub page](https://github.com/jqnatividad/qsv).

Note: `make deps` will install these tools, but not the system-level software, which you will have to do manually. After installing these tools, make sure their binaries are located in a directory included in your system's `PATH`.

## 2. Python Environment Setup

It is highly recommended to use a Python virtual environment to manage the project's dependencies. You can run `setup_env.sh` which will install packages to `_venv` within Arborist root directory or you can run the following:

**Create and activate a virtual environment:**
```bash
python3 -m venv _venv
source _venv/bin/activate
pip install -r src/protein/requirements.txt
```

This will install all necessary libraries, including `polars`, `pandas`, `pepmatch`, and the ARC toolkit directly from its Git repository.

## 3. Database Configuration

The Arborist project needs to connect to an IEDB MySQL database to fetch the latest data for its builds. The `src/util/mysql2tsv` script requires the following environment variables to be set either manually or in your `.bashrc` file:

```bash
export IEDB_MYSQL_HOST='database_host'
export IEDB_MYSQL_PORT='port_#'
export IEDB_MYSQL_USER='your_username'
export IEDB_MYSQL_PASSWORD='your_password'
export IEDB_MYSQL_DATABASE='database_name'
```