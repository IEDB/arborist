import csv

from argparse import ArgumentParser
from collections import defaultdict


def read_core(core_path):
    '''Read upper.tsv and return:
    a 'rows' dictionary from ID to row
    and a 'tree' dictionary from ID to a list of children.'''
    rows = {}
    tree = defaultdict(list)
    with open(core_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            curie = row['curie']
            rows[curie] = row
            parent = row['parent']
            if parent:
                tree[parent].append(curie)
    return rows, tree


def sort_rows(rows, tree, curie):
    '''Given rows, tree, and a CURIE,
    walk the tree depth-first
    and return the sorted list of rows.'''
    output = [rows[curie]]
    if len(tree[curie]) > 0:
        for child in tree[curie]:
            output += sort_rows(rows, tree, child)
    return output


def main():
    parser = ArgumentParser('Sort organism_core.tsv in tree-order')
    parser.add_argument('core', help='Path to organism_core.tsv', type=str)
    args = parser.parse_args()

    rows, tree = read_core(args.core)
    sorted_rows = sort_rows(rows, tree, 'NCBITaxon:1')
    fieldnames = list(sorted_rows[0].keys())

    missing = set(rows.keys()) - set([r['curie'] for r in sorted_rows])
    missing_parents = {}
    if len(missing) > 0:
        for curie in missing:
            row = rows[curie]
            missing_parents[row['parent']] = row['parent_label']
            print(f'Row not under root {curie} {row["label"]}')
            sorted_rows.append(rows[curie])

    print('MISSING PARENTS')
    for curie, label in missing_parents.items():
        print(curie, label)

    with open(args.core, 'w') as f:
        writer = csv.DictWriter(f, fieldnames,
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(sorted_rows)


if __name__ == "__main__":
    main()
