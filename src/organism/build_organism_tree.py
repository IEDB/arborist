import csv
import logging
import sqlite3

from argparse import ArgumentParser, FileType
from assign_species import get_taxon_id


datatypes = ['_IRI', 'xsd:string', 'xsd:boolean', 'xsd:integer']


def insert_triples(con, table, triples):
    '''Given a SQLite database connection, a table name,
    and a list of triples to insert,
    insert those triples into the database table and commit.'''
    for triple in triples:
        if len(triple) < 3:
            raise Exception('Not enough slots in triple: ' + str(triple))
        elif len(triple) == 3:
            if ':' not in triple[2]:
                raise Exception(
                    'Object without datatype must be an IRI or CURIE: ',
                    triple[2]
                )
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


def insert_annotations(con, table):
    '''Insert triples for annotation properties.'''
    annotations = {
        'ONTIE:0003615': 'has taxon ID',
        'ONTIE:0003616': 'NCBI Taxonomy browser',
        'ONTIE:0003617': 'has taxonomic rank',
        'ONTIE:0003618': 'used in IEDB',
        'oio:hasAlternativeId': 'has_alternative_id',
        'oio:hasExactSynonym': 'has_exact_synonym',
        'oio:hasBroadSynonym': 'has_broad_synonym',
        'oio:hasRelatedSynonym': 'has_related_synonym',
        'oio:hasLabelSource': 'has_label_source',
        'oio:hasSynonymType': 'has_synonyms_type',
    }
    triples = []
    for curie, label in annotations.items():
        triples += [
            [curie, 'rdf:type', 'owl:AnnotationProperty'],
            [curie, 'rdfs:label', label, 'xsd:string'],
        ]
    insert_triples(con, table, triples)


def index_statement_table(con, table):
    '''Given a database connection and a table name
    create the usual LDTab indexes for that table.'''
    cur = con.cursor()
    cur.execute(f'CREATE INDEX idx_{table}_subject ON {table}(subject)')
    cur.execute(f'CREATE INDEX idx_{table}_predicate ON {table}(predicate)')
    cur.execute(f'CREATE INDEX idx_{table}_object ON {table}(object)')
    cur.execute(f'ANALYZE {table}')
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
    exclude = [
        'rdfs:subClassOf',
        'rdf:type',
        'rdfs:label',
        'oio:hasDbXref',
        'oio:hasOBONamespace',
        'ncbitaxon:has_rank',
    ]
    exclude = '","'.join(exclude)
    exclude = f'"{exclude}"'
    con.execute(f'''
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
      AND predicate NOT IN ({exclude})
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

    synonym_type = (
        '{"oio:hasSynonymType":['
        '{"datatype":"xsd:string","meta":"owl:Axiom","object":"IEDB"}'
        ']}'
    )
    with sqlite3.connect(args.ldtab) as con:
        create_statement_table(con, 'organism_tree')
        insert_annotations(con, 'organism_tree')
        triples = []
        for row in csv.DictReader(args.organism_tree, delimiter='\t'):
            curie = row['curie']
            triples += [
                [curie, 'rdf:type', 'owl:Class'],
                [curie, 'rdfs:label', row['label'], 'xsd:string'],
                [curie, 'iedb-taxon:label-source', row['label_source'],
                    'xsd:string'],
                [curie, 'iedb-taxon:level', row['level'], 'xsd:string'],
                [curie, 'iedb-taxon:source-table',
                    row['source_table'], 'xsd:string'],
            ]
            taxon_id = get_taxon_id(curie)
            if taxon_id:
                triples.append(
                    [curie, 'ONTIE:0003615', taxon_id, 'xsd:string'])
                if curie.startswith('NCBITaxon:'):
                    url = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/'
                    url += f'wwwtax.cgi?id={taxon_id}'
                    triples.append([curie, 'ONTIE:0003616', url])
            if row['iedb_synonyms']:
                synonyms = [s.strip() for s in row['iedb_synonyms'].split(';')]
                for synonym in synonyms:
                    triples.append([
                        curie,
                        'oio:hasExactSynonym',
                        synonym,
                        'xsd:string',
                        synonym_type,
                    ]),
            if row['epitope_count']:
                triples += [
                    [curie, 'iedb-taxon:epitope-count',
                     row['epitope_count'], 'xsd:integer'],
                    [curie, 'ONTIE:0003618', 'true', 'xsd:boolean'],
                ]
            if row['rank']:
                triples.append(
                    [curie, 'ONTIE:0003617', row['rank'], 'xsd:string'])
            if row['parent']:
                triples.append([curie, 'rdfs:subClassOf', row['parent']])
            if row['parent2']:
                triples.append([curie, 'rdfs:subClassOf', row['parent2']])
        insert_triples(con, 'organism_tree', triples)
        copy_triples(con)
        index_statement_table(con, 'organism_tree')
        # args.active_taxa.seek(0)
        # check_lower(con, args.active_taxa)


if __name__ == "__main__":
    main()
