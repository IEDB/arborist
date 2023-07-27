import csv
import logging
import sqlite3

from argparse import ArgumentParser, FileType
from assign_species import get_curie
from collections import defaultdict


datatypes = ['_IRI', 'xsd:string', 'xsd:boolean', 'xsd:integer']


def write_triples(con, table, triples):
    '''Given a SQLite database connection, a table name,
    and a list of triples to insert,
    insert those triples into the database table and commit.'''
    for triple in triples:
        if len(triple) < 3:
            raise Exception('Not enough slots in triple: ' + str(triple))
        elif len(triple) == 3:
            if ':' not in triple[2]:
                raise Exception(
                    'Object without datatype must be an IRI or CURIE: ', triple[2])
            triple.append('_IRI')
        elif len(triple) == 4 and triple[3] not in datatypes:
            raise Exception('Unexpected datatype in: ' + str(triple))

        if len(triple) == 4:
            triple.append(None)

        if len(triple) > 5:
            raise Exception('Too many slots in triple: ' + str(triple))

        if triple[1] == 'rdfs:label' and triple[3] != 'xsd:string':
            raise Exception(f'Bad datatype for label {triple[2]}^^{triple[3]}')

        # print('\t'.join([str(x) for x in triple]))

    cur = con.cursor()
    cur.executemany(
        f'INSERT INTO {table} '
        f'VALUES(1, 0, "iedb-taxon:{table}", ?, ?, ?, ?, ?)',
        triples
    )
    con.commit()


def create_statement_table(con, table):
    '''Given a database connection and a table name,
    create the LDTab table strucure for that table.'''
    con.execute(f'DROP TABLE IF EXISTS {table}')
    con.execute(f'''CREATE TABLE {table} (
  'assertion' INTEGER NOT NULL,
  'retraction' INTEGER NOT NULL DEFAULT 0,
  'graph' TEXT NOT NULL,
  'subject' TEXT NOT NULL,
  'predicate' TEXT NOT NULL,
  'object' TEXT NOT NULL,
  'datatype' TEXT NOT NULL,
  'annotation' TEXT
)''')
    con.commit()


def index_statement_table(con, table):
    '''Given a database connection and a table name
    create the usual LDTab indexes for that table.'''
    cur = con.cursor()
    cur.execute(f'CREATE INDEX idx_{table}_subject ON {table}(subject)')
    cur.execute(f'CREATE INDEX idx_{table}_predicate ON {table}(predicate)')
    cur.execute(f'CREATE INDEX idx_{table}_object ON {table}(object)')
    # cur.execute('ANALYZE')
    con.commit()


def check_lower(con, active_taxa):
    '''Given a database connection and a file handle to read active_taxa.tsv
    check that they can all be found in the 'organism_tree' table.'''
    all_found = True
    for row in csv.DictReader(active_taxa, delimiter='\t'):
        curie = row['curie']
        cur = con.execute(
            'SELECT count(*) FROM organism_tree '
            'WHERE subject = ? AND predicate = "rdfs:subClassOf"',
            (curie,)
        )
        res = cur.fetchone()
        if res[0] == 0:
            print(f'Missing node for active taxon {curie}')
            all_found = False
    return all_found


def copy_triples(con):
    '''Given a database connection,
    for each subject in the 'organism_tree' table,
    copy the relevant triples from the 'ncbitaxon' table.'''
    con.execute('''
    INSERT INTO organism_tree
    SELECT
        1 AS asserted,
        0 AS retracted,
        "iedb-taxon:organism_tree" AS graph,
        subject,
        predicate,
        object,
        datatype,
        annotation
    FROM ncbitaxon
    WHERE subject IN (SELECT subject FROM organism_tree)
      AND predicate NOT IN ("rdfs:subClassOf", "rdf:type", "rdfs:label")
    ''')


def main():
    parser = ArgumentParser('Build the organism tree')
    parser.add_argument('ldtab', help='Path to the LDTab SQLite database')
    parser.add_argument(
        'organism_tree',
        help='Path to organism-tree.tsv to read',
        type=FileType('r'),
    )
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

    with sqlite3.connect(args.ldtab) as con:
        create_statement_table(con, 'organism_tree')
        triples = []
        for row in csv.DictReader(args.organism_tree, delimiter='\t'):
            curie = row['curie']
            triples += [
                [curie, 'rdf:type', 'owl:Class'],
                [curie, 'rdfs:label', row['label'], 'xsd:string'],
                [curie, 'iedb-taxon:level', row['level'], 'xsd:string'],
                [curie, 'iedb-taxon:source-table', row['source_table'], 'xsd:string'],
            ]
            if row['epitope_count']:
                [curie, 'iedb-taxon:epitope-count', row['epitope_count'], 'xsd:integer'],
            if row['rank']:
                triples.append([curie, 'ncbitaxon:has-rank', row['rank'], 'xsd:string'])
            if row['parent']:
                triples.append([curie, 'rdfs:subClassOf', row['parent']])
        write_triples(con, 'organism_tree', triples)
        copy_triples(con)
        index_statement_table(con, 'organism_tree')
        # args.active_taxa.seek(0)
        # check_lower(con, args.active_taxa)


if __name__ == "__main__":
    main()
