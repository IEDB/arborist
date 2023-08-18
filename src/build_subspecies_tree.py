import logging
import sqlite3

from argparse import ArgumentParser
from build_organism_tree import (
    create_statement_table, index_statement_table
)


def main():
    parser = ArgumentParser('Build the organism tree')
    parser.add_argument('ldtab', help='Path to the LDTab SQLite database')
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Turn on logging'
    )
    args = parser.parse_args()

    # Set up logging
    log_format = '%(levelname)s: %(message)s'
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format=log_format)
    else:
        logging.basicConfig(level=logging.WARNING, format=log_format)

    con = sqlite3.connect(args.ldtab)
    create_statement_table(con, 'subspecies_tree')

    # Copy all triples from organism_tree
    con.execute('''
    INSERT INTO subspecies_tree
    SELECT
        1 AS asserted,
        0 AS retracted,
        "iedb-taxon:subspecies_tree" AS graph,
        subject,
        predicate,
        object,
        datatype,
        annotation
    FROM organism_tree''')

    # For all species in the organism_tree
    # copy all their triples from ncbitaxon into subspecies_tree.
    cur = con.execute(f'''
WITH RECURSIVE descendants(node) AS (
      SELECT subject AS node
      FROM organism_tree
      WHERE predicate = "iedb-taxon:level"
        AND object = "species"
    UNION
      SELECT subject AS node
      FROM ncbitaxon, descendants
      WHERE descendants.node = ncbitaxon.object
        AND ncbitaxon.predicate = "rdfs:subClassOf"
), new_descendants(node) AS (
      SELECT * FROM descendants
    EXCEPT
      SELECT subject FROM organism_tree
)
INSERT INTO subspecies_tree
SELECT
    1 AS asserted,
    0 AS retracted,
    "iedb-taxon:subspecies_tree" AS graph,
    subject,
    predicate,
    object,
    datatype,
    annotation
FROM ncbitaxon, new_descendants
WHERE ncbitaxon.subject = new_descendants.node
    ''')
    con.commit()

    index_statement_table(con, 'subspecies_tree')
    con.execute('ANALYZE')

    # Add ONTIE:0003619 'has species' annotation property
    con.execute('''
    INSERT INTO subspecies_tree VALUES
    (1, 0, 'iedb-taxon:subspecies_tree', 'ONTIE:0003619', 'rdf:type', 'owl:AnnotationProperty', '_IRI', NULL),
    (1, 0, 'iedb-taxon:subspecies_tree', 'ONTIE:0003619', 'rdfs:label', 'has species', 'xsd:string', NULL);
    ''')

    # Add ONTIE:0003619 'has species' annotations
    cur = con.execute(f'''
WITH RECURSIVE species_descendants(species, node) AS (
      SELECT subject AS species, subject AS node
      FROM subspecies_tree
      WHERE predicate = 'iedb-taxon:level'
        AND object = 'species'
    UNION
      SELECT species, subject AS node
      FROM subspecies_tree, species_descendants
      WHERE species_descendants.node = subspecies_tree.object
        AND subspecies_tree.predicate = 'rdfs:subClassOf'
)
INSERT INTO subspecies_tree
SELECT
    1 AS asserted,
    0 AS retracted,
    'iedb-taxon:subspecies_tree' AS graph,
    node AS subject,
    'ONTIE:0003619' AS predicate,
    species AS object,
    '_IRI' AS datatype,
    NULL AS annotation
FROM species_descendants
WHERE species != node
    ''')

    con.commit()


if __name__ == "__main__":
    main()
