#!/usr/bin/env python3

# Serve HCC-KB development and curation site,
# calling Nanobot to do the hard work.

import csv
import os
import subprocess

from bottle import get, post, request, response, run, HTTPError
from jinja2 import Environment, FileSystemLoader

env = Environment(
    loader=FileSystemLoader('src/templates/')
)


@get('/')
def index():
    # return template('<b>Hello {{name}}</b>!', name=name)
    output = [
        '<h3>Welcome to Arborist!</h3>',
        '<p>Arborist builds trees for the IEDB '
        'using IEDB data, community standard ontologies, '
        'and upstream databases</p>',
        '<ul>',
        '  <li>',
        '    <a href="/iedb/table">IEDB</a> upstream tables:',
        '    <ul>',
    ]
    iedb = [
        'ncbi_include',
        'iedb_taxa',
        'source',
        'object',
        'epitope',
        'peptide',
        'peptide_source'
    ]
    for entry in iedb:
        output.append(f'    <li><a href="/iedb/{entry}">{entry}</a></li>')
    output += [
        '    </ul>',
        '  </li>',
        '  <li>',
        '    <a href="/arborist/table">Arborist</a> tables and trees: ',
        '    <ul>',
    ]
    arborist = {
        'NCBI Taxonomy': 'ncbitaxon/NCBITaxon:1',
        'Organism Core Table': 'organism_core',
        'Organism Tree Table': 'organism_tree_tsv',
        'Organism Tree': 'organism_tree/NCBITaxon:1',
        'Subspecies Tree': 'subspecies_tree/NCBITaxon:1',
        'Active Species': 'active_species',
        'Proteome': 'proteome',
        'Protein Tree': 'protein_tree'
    }
    for name, href in arborist.items():
        output.append(f'      <li><a href="/arborist/{href}">{name}</a></li>')
    output += [
        '    </ul>',
        '  </li>',
        '</ul>',
    ]
    return render('\n'.join(output))


@get('/iedb/<path:path>')
def get_iedb(path):
    return nanobot('GET', 'iedb', path)


@post('/iedb/<path:path>')
def post_iedb(path):
    return nanobot('POST', 'iedb', path)


@get('/arborist/<path:path>')
def get_arborist(path):
    return nanobot('GET', 'arborist', path)


@post('/arborist/<path:path>')
def post_arborist(path):
    return nanobot('POST', 'arborist', path)


def nanobot(method, dataset, path):
    """Call Nanobot as a CGI script
    for the given dataset, and path."""
    if method not in ['GET', 'POST']:
        raise HTTPError(f'Bad method {method}')

    result = subprocess.run(
        [os.path.join(os.getcwd(), 'bin/nanobot')],
        cwd=f'build/{dataset}/',
        env={
            'GATEWAY_INTERFACE': 'CGI/1.1',
            'REQUEST_METHOD': method,
            'PATH_INFO': path,
            'QUERY_STRING': request.query_string,
        },
        input=request.body.getvalue().decode('utf-8'),
        text=True,
        capture_output=True
    )
    reading_headers = True
    body = []
    for line in result.stdout.splitlines():
        if reading_headers and line.strip() == '':
            reading_headers = False
            continue
        if reading_headers:
            name, value = line.split(': ', 1)
            if name == 'status':
                response.status = value
            else:
                response.set_header(name, value)
        else:
            body.append(line)
    return '\n'.join(body)


def render(content):
    """Given an HTML content string,
    render the `page.html` template
    and return the resulting HTML string."""
    template = env.from_string('''{% extends "page.html" %}
{% block content %}
    {{ content }}
{% endblock %}''')
    return template.render(
        page={
            'project_name': 'Arborist',
            'tables': {
            },
        },
        content=content
    )


run(host='0.0.0.0', port=3000)
