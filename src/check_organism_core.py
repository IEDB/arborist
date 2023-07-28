import logging
import sqlite3

from argparse import ArgumentParser
from assign_species import (
    get_predicate, get_synonyms, get_rank, get_level, get_species, get_label,
    synonym_types
)


def check_tree(con, messages, row):
    curie = row['curie']
    if curie == 'NCBITaxon:1':
        return

    current = curie
    found_root = False
    while current:
        # TODO: use WITH RECURSIVE query
        parent = get_predicate(
            con,
            'rdfs:subClassOf',
            current,
            table='organism_tree'
        )
        if not parent:
            break
        if parent == 'NCBITaxon:1':
            found_root = True
            break
        if parent == current:
            print(f'Circular reference for {current}')
            break
        current = parent

    if found_root:
        return
    messages.append({
        'table': 'organism_core',
        'row': row['row_number'],
        'column': 'curie',
        'value': curie,
        'level': 'error',
        'rule': 'check_organism_core: no path to root',
        'message': 'Could not find a path to the root',
    })


def check_level(con, messages, row):
    if row['level']:
        return

    curie = row['curie']
    rank = get_rank(con, curie)
    level = get_level(rank)
    if rank and level:
        messages.append({
            'table': 'organism_core',
            'row': row['row_number'],
            'column': 'level',
            'value': '',
            'level': 'info',
            'rule': 'check_organism_core: level suggestion',
            'message': f'Found rank "{rank}" so level should be: {level}',
        })
        return

    species = get_species(con, {}, curie)
    if species:
        species_label = get_label(con, {}, species)
        messages.append({
            'table': 'organism_core',
            'row': row['row_number'],
            'column': 'level',
            'value': '',
            'level': 'info',
            'rule': 'check_organism_core: level suggestion',
            'message': (
                f'Found ancestor species "{species_label}" '
                'so level should be: lower'
            ),
        })


def check_label(con, messages, row):
    curie = row['curie']
    label = get_predicate(con, 'rdfs:label', curie)
    if label:
        return
    messages.append({
        'table': 'organism_core',
        'row': row['row_number'],
        'column': 'curie',
        'value': curie,
        'level': 'error',
        'rule': 'check_organism_core: missing taxon',
        'message': (
            'This taxon is not in the latest NCBI Taxonomy. '
            'Maybe it was merged.'
        ),
    })


def check_label_source_missing(con, messages, row):
    if row['label_source']:
        return

    curie = row['curie']
    label = get_predicate(con, 'rdfs:label', curie)
    if row['label'] == label:
        messages.append({
            'table': 'organism_core',
            'row': row['row_number'],
            'column': 'label_source',
            'value': '',
            'level': 'info',
            'rule': 'check_organism_core: label_source suggestion',
            'message': 'Label source should be: NCBI Taxonomy scientific name',
        })
        return

    synonyms = get_synonyms(con, curie)
    for synonym_type, syns in synonyms.items():
        for synonym in syns:
            new_label = f'{label} ({synonym})'
            syn_type = synonym_type
            if synonym_type in synonym_types:
                syn_type = synonym_types[synonym_type]
            if new_label == row['label']:
                messages.append({
                    'table': 'organism_core',
                    'row': row['row_number'],
                    'column': 'label_source',
                    'value': row['label_source'],
                    'level': 'info',
                    'rule': 'check_organism_core: label suggestion',
                    'message': (
                        'Label source should be: '
                        'NCBI Taxonony scientific name '
                        f'({syn_type})'
                    ),
                })


def check_label_source_ncbi(con, messages, row):
    if row['label_source'] != 'NCBI Taxonomy scientific name':
        return

    curie = row['curie']
    label = get_predicate(con, 'rdfs:label', curie)
    if row['label'] == label:
        return

    messages.append({
        'table': 'organism_core',
        'row': row['row_number'],
        'column': 'label',
        'value': row['label'],
        'level': 'warn',
        'rule': 'check_organism_core: label mismatch',
        'message': (
            'Mismatching NCBI Taxonomy scientific name: '
            f'expected "{row["label"]}" found "{label}"'
        )
    })

    synonyms = get_synonyms(con, curie)
    for synonym_type, syns in synonyms.items():
        for synonym in syns:
            new_label = f'{label} ({synonym})'
            syn_type = synonym_type
            if synonym_type in synonym_types:
                syn_type = synonym_types[synonym_type]
            if new_label == row['label']:
                messages.append({
                    'table': 'organism_core',
                    'row': row['row_number'],
                    'column': 'label_source',
                    'value': row['label_source'],
                    'level': 'info',
                    'rule': 'check_organism_core: label suggestion',
                    'message': (
                        'Label source should be: '
                        'NCBI Taxonony scientific name '
                        f'({syn_type})'
                    ),
                })


def main():
    parser = ArgumentParser()
    parser.add_argument('ldtab', help='Path to the LDTab SQLite database')
    parser.add_argument('-u', '--update',
                        action='store_true', help='Accept all suggestions')
    parser.add_argument('-v', '--verbose',
                        action='store_true', help='Turn on logging')
    args = parser.parse_args()

    # Set up logging
    log_format = '%(levelname)s: %(message)s'
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format=log_format)
    else:
        logging.basicConfig(level=logging.WARNING, format=log_format)

    con = sqlite3.connect(args.ldtab)

    con.execute('DELETE FROM message WHERE rule LIKE "check_organism_core:%"')
    con.commit()

    messages = []
    cur = con.execute('SELECT * FROM organism_core_view')
    for row in cur.fetchall():
        fields = [column[0] for column in cur.description]
        row = {key: value for key, value in zip(fields, row)}

        curie = row['curie']
        if not curie.startswith('NCBITaxon:'):
            continue

        check_tree(con, messages, row)
        check_level(con, messages, row)
        check_label(con, messages, row)
        check_label_source_missing(con, messages, row)
        check_label_source_ncbi(con, messages, row)

    updates = 0
    if args.update:
        remaining_messages = []
        for message in messages:
            if message['level'] != 'info':
                remaining_messages.append(message)
                continue
            if 'should be: ' not in message['message']:
                remaining_messages.append(message)
                continue
            updates += 1
            replacement = message['message'].split('should be: ', 1)[1]
            con.execute(
                'UPDATE organism_core '
                f'SET \'{message["column"]}\' = "{replacement}" '
                f'WHERE row_number = {message["row"]}'
            )
            con.commit()
        print('Updates:', updates)
        messages = remaining_messages

    print('Messages:', len(messages))
    if args.verbose:
        import csv
        fieldnames = [
            'table',
            'row',
            'column',
            'value',
            'level',
            'rule',
            'message',
        ]
        with open('build/messages.tsv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames,
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(messages)

    con.executemany(
        'INSERT INTO message '
        "('table', 'row', 'column', 'value', 'level', 'rule', 'message') "
        'VALUES (:table, :row, :column, :value, :level, :rule, :message)',
        messages
    )
    con.commit()

    # writer = csv.DictWriter(args.output, fieldnames,
    #                         delimiter='\t', lineterminator='\n',
    #                         extrasaction='ignore')
    # writer.writeheader()
    # writer.writerows(tree.values())
    # writer.writerows(assigned.values())


if __name__ == "__main__":
    main()
