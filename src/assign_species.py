import csv
import json
import logging
import sqlite3

from argparse import ArgumentParser, FileType


synonym_types = {
    'common_name': 'NCBI Taxonomy common name',
    'genbank_common_name': 'Genbank common name',
    'equivalent_name': 'NCBI Taxonomy equivalent name',
}

upper_ranks = [
    'clade',
    'class',
    'family',
    'genus',
    'species_group',
    'infraorder',
    'kingdom',
    'order',
    'parvorder',
    'phylum',
    'species_group',
    'subclass',
    'subfamily',
    'suborder',
    'subphylum',
    'superclass',
    'superkingdom',
]

lower_ranks = [
    'biotype',
    'cohort',
    'section',
    'series',
    'serogroup',
    'serotype',
    'strain',
    'subcohort',
    'subsection',
    'subspecies',
    'subtribe',
    'tribe',
    "varietas",
]


def get_curie(id):
    '''Given an id string, return a NCBITaxon, iedb-taxon, or OBO CURIE.'''
    try:
        if int(id) >= 10000000:
            return f'iedb-taxon:{id}'
        else:
            return f'NCBITaxon:{id}'
    except ValueError:
        if id.endswith("-other"):
            return f'iedb-taxon:{id}'
        return id


def get_taxon_id(curie):
    curie = curie.replace('NCBITaxon:', '')
    curie = curie.replace('iedb-taxon:', '')
    curie = curie.replace('-other', '')
    return curie


def get_other_curie(curie):
    taxon_id = get_taxon_id(curie)
    return f'iedb-taxon:{taxon_id}-other'


def get_predicate(con, predicate, curie, table='ncbitaxon'):
    cur = con.execute(
        f"SELECT object FROM '{table}' WHERE subject = ? AND predicate = ?",
        (curie, predicate)
    )
    res = cur.fetchone()
    if res:
        return res[0]


def get_parents(con, tree, curie):
    parents = []
    if curie in tree:
        row = tree[curie]
        if 'parent' in row and row['parent']:
            parents.append(row['parent'])
        if 'parent2' in row and row['parent2']:
            parents.append(row['parent2'])
    if not parents:
        cur = con.execute(
            "SELECT object FROM 'ncbitaxon' "
            'WHERE subject = ? AND predicate = "rdfs:subClassOf"',
            (curie,)
        )
        for res in cur.fetchall():
            parents.append(res[0])
    return parents


def get_parent(con, tree, curie):
    ps = get_parents(con, tree, curie)
    if ps and len(ps) > 0:
        return ps[0]
    return None


def get_parent_in_tree(con, tree, curie):
    current = curie
    while current:
        parent = get_parent(con, tree, current)
        if not parent:
            return None
        if parent in tree and 'use_other' in tree[parent] \
                and tree[parent]['use_other'] == 'TRUE':
            return get_other_curie(parent)
        if parent in tree:
            return parent
        current = parent


def get_species(con, tree, curie):
    current = curie
    while current:
        parent = get_parent(con, tree, current)
        if not parent:
            return None
        if parent in tree and 'level' in tree[parent] \
                and tree[parent]['level'] == 'species':
            return parent
        rank = get_rank(con, parent)
        if rank == 'species':
            return parent
        current = parent


def get_synonyms(con, curie):
    predicates = ['rdfs:label', 'oio:hasExactSynonym',
                  'oio:hasRelatedSynonym', 'oio:hasBroadSynonym']
    predicates = '","'.join(predicates)
    predicates = f'"{predicates}"'
    cur = con.execute(
        'SELECT object, annotation '
        'FROM ncbitaxon '
        f'WHERE subject = ? AND predicate IN ({predicates})'
        'ORDER BY object',
        (curie,)
    )

    # Collect the synonyms by type and in alphabetical order
    synonyms = {}
    for result in cur.fetchall():
        synonym = result[0]
        try:
            annotation = json.loads(result[1])
            synonym_type = annotation['oio:hasSynonymType'][0]['object']
        except Exception:
            continue
        if not synonym_type:
            continue
        synonym_type = synonym_type.replace('ncbitaxon:', '')
        if synonym_type not in synonyms:
            synonyms[synonym_type] = [synonym]
        else:
            synonyms[synonym_type].append(synonym)

    return synonyms


def get_label_and_source(con, tree, curie):
    label = None
    label_source = None
    if curie in tree:
        label = tree[curie]['label']
        label_source = tree[curie]['label_source']
    else:
        # TODO: common names, etc.
        label = get_predicate(con, 'rdfs:label', curie)
        label_source = 'NCBI Taxonomy scientific name'

        # Pick the best synonym type available
        synonyms = get_synonyms(con, curie)
        for type, type_label in synonym_types.items():
            if type in synonyms:
                synonyms[type].sort()
                synonym = synonyms[type][0]
                label = f'{label} ({synonym})'
                label_source = (
                    f'NCBI Taxonomy scientific name ({type_label})'
                )
                break

    if not label:
        raise Exception(f'No label found for {curie}')

    return label, label_source


def get_label(con, tree, curie):
    label, label_source = get_label_and_source(con, tree, curie)
    return label


def get_rank(con, curie):
    rank = get_predicate(con, 'ncbitaxon:has_rank', curie)
    if rank:
        rank = rank.replace('NCBITaxon:', '').replace(
            '<http://purl.obolibrary.org/obo/NCBITaxon#_', '').replace('>', '')
    return rank


def get_level(rank):
    if rank in upper_ranks:
        return 'upper'
    elif rank == 'species':
        return 'species'
    elif rank in lower_ranks:
        return 'lower'


def main():
    parser = ArgumentParser()
    parser.add_argument('ldtab', help='Path to the LDTab SQLite database')
    parser.add_argument(
        'core',
        help='Path to organism_core TSV',
        type=FileType('r'),
    )
    parser.add_argument(
        'iedb',
        help='Path to IEDB taxon file',
        type=FileType('r'),
    )
    parser.add_argument(
        'count',
        help='Path to count.tsv file',
        type=FileType('r'),
    )
    parser.add_argument(
        'output',
        help='Path to write',
        type=FileType('w'),
    )
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
    tree = {}

    # Add terms from organism_core.tsv
    others_to_check = set()
    for row in csv.DictReader(args.core, delimiter='\t'):
        curie = row['curie']
        tree[curie] = row
        tree[curie].update({
            'rank': get_rank(con, curie),
            'source_table': 'organism_core',
        })
        # TODO: more checks
        if row['use_other'] == 'TRUE':
            others_to_check.add(curie)

    for parent in others_to_check:
        curie = get_other_curie(parent)
        if curie not in tree:
            raise Exception(
                f'Missing "Other" node for {parent} {tree[parent]["label"]}'
            )

    # Add IEDB taxa
    iedb_taxa = []
    for row in csv.DictReader(args.iedb, delimiter='\t'):
        curie = get_curie(row['iedb_id'])
        # Maybe the IEDB taxon is already in organism_core
        if curie in tree:
            if row['label'] != tree[curie]['label']:
                print(
                    f'Mismatching labels for {curie}:',
                    f'"{tree[curie]["label"]}" vs "{row["label"]}"'
                )
            # if get_curie(row['parent_ids']) != tree[curie]['parent']:
            #     taxon_id = get_taxon_id(tree[curie]['parent'])
            #     print(
            #         f'Mismatching parents for {curie}:',
            #         f'{taxon_id} vs {row["parent_ids"]}'
            #     )
            continue
        iedb_taxa.append(curie)
        rank = row['rank']

        parents = []
        for id in row['parent_ids'].split(','):
            parents.append(get_curie(id))
        parent = parents[0]

        tree[curie] = {
            'curie': curie,
            'label': row['label'],
            'label_source': 'IEDB',
            'rank': rank,
            'level': get_level(rank),
            'parent': parent,
            'iedb_synonyms': row['synonyms'],
            'source_table': 'iedb_taxa',
        }
        if len(parents) == 2:
            tree[curie]['parent2'] = parents[1]
        elif len(parents) != 1:
            print(
                f'Wrong number of parents for {curie} {row["label"]}:',
                parents
            )

    # Add IEDB taxon parents and update parent labels
    for curie in iedb_taxa:
        row = tree[curie]
        parent = row['parent']
        if parent not in tree:
            label, label_source = get_label_and_source(con, tree, parent)
            rank = get_rank(con, parent)
            tree[parent] = {
                'curie': parent,
                'label': label,
                'label_source': label_source,
                'rank': rank,
                'level': get_level(rank),
                'source_table': 'iedb_taxa.parent',
            }
        row['parent_label'] = tree[parent]['label']

        if 'parent2' not in row or not row['parent2']:
            continue

        parent2 = row['parent2']
        if parent2 not in tree:
            label, label_source = get_label_and_source(con, tree, parent2)
            rank = get_rank(con, parent2)
            tree[parent2] = {
                'curie': parent2,
                'label': label,
                'label_source': label_source,
                'rank': rank,
                'level': get_level(rank),
                'source_table': 'iedb_taxa.parent2',
            }
        row['parent2_label'] = tree[parent2]['label']

    # Add all active taxa, and assign epitope counts
    for row in csv.DictReader(args.count, delimiter='\t'):
        curie = get_curie(row['source_organism_org_id'])
        if curie in tree:
            tree[curie]['epitope_count'] = row['count']
        elif curie.startswith('iedb-taxon:'):
            raise Exception(
                'IEDB taxon should already have been added '
                f'{curie} {row.get("label") or "UNKNOWN"}'
            )
        else:
            label, label_source = get_label_and_source(con, tree, curie)
            rank = get_rank(con, curie)
            tree[curie] = {
                'curie': curie,
                'label': label,
                'label_source': label_source,
                'rank': rank,
                'level': get_level(rank),
                'epitope_count': row['count'],
                'source_table': 'count',
            }

    # Assign species to lower taxa
    species_to_add = set()
    for curie, row in tree.items():
        if 'level' in row and row['level'] != 'lower':
            continue

        if 'rank' in row:
            rank = row['rank']
        else:
            rank = get_rank(con, curie)

        if rank in upper_ranks:
            continue
        elif rank == 'species':
            species = curie
        else:
            species = get_species(con, tree, curie)

        if species:
            row['species'] = species
            row['species_label'] = get_label(con, tree, species)
            if species not in tree:
                species_to_add.add(species)

    # Add any missing species
    for species in species_to_add:
        label, label_source = get_label_and_source(con, tree, species)
        rank = get_rank(con, species)
        tree[species] = {
            'curie': species,
            'label': label,
            'label_source': label_source,
            'rank': rank,
            'level': get_level(rank),
            'species': species,
            'species_label': label,
            'source_table': 'species not otherwise in use',
        }

    # Assign parents that have not yet been assigned
    for curie, row in tree.items():
        if 'parent' in row and row['parent']:
            continue
        parent = get_parent_in_tree(con, tree, curie)
        if not parent:
            if curie != 'NCBITaxon:1':
                print(f'No parent found for {curie} {row["label"]} {row}')
            continue
        tree[curie].update({
            'parent': parent,
            'parent_label': get_label(con, tree, parent),
        })

    fieldnames = [
        'curie', 'label', 'label_source',
        'iedb_synonyms',
        'rank', 'level',
        'epitope_count',
        'parent', 'parent_label',
        'parent2', 'parent2_label',
        'species', 'species_label',
        'source_table', 'use_other',
    ]
    writer = csv.DictWriter(args.output, fieldnames,
                            delimiter='\t', lineterminator='\n',
                            extrasaction='ignore')
    writer.writeheader()
    writer.writerows(tree.values())
    # writer.writerows(assigned.values())


if __name__ == "__main__":
    main()
