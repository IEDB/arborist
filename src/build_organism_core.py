from argparse import ArgumentParser, FileType
from assign_species import get_curie
from sort_organism_core import read_core


def render_tree(rows, tree, curie, depth):
    '''Given rows, tree, a CURIE, and a depth,
    walk the tree depth-first and
    return a list of lines for a nested HTML unordered list.'''
    indent = '  ' * depth
    label = rows[curie]['label']
    output = [f'{indent}<li><a href="../organism_tree/{curie}">{label}</a></li>']
    if len(tree[curie]) > 0:
        output += [f'{indent}<ul>']
        for child in tree[curie]:
            output += render_tree(rows, tree, child, depth + 1)
        output += [f'{indent}</ul>']
    return output


def main():
    parser = ArgumentParser('Render upper.tsv to a nested HTML list')
    parser.add_argument('core', help='Path to organism_core.tsv', type=str)
    parser.add_argument(
        'output',
        help='Path to HTML file to write',
        type=FileType('w'),
    )
    args = parser.parse_args()

    rows, tree = read_core(args.core)
    output = [
        '<html>',
        '<head>',
        '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">',
        '<style>a:link, a:visited { text-decoration: none !important; }</style>',
        '</head>',
        '<body>',
        '<div class="container">',
        '<h1>Upper Organism Tree</h1>',
        '<p>Generated from the <a href="../upper">upper</a> table to build the new <a href="../organism_tree/NCBITaxon:1">organism tree</a>.</p>',
        '<ul>'
    ]
    output += render_tree(rows, tree, 'NCBITaxon:1', 1)
    output += ['</ul>', '</div>', '</body>', '</html>']
    args.output.write('\n'.join(output))


if __name__ == "__main__":
    main()
