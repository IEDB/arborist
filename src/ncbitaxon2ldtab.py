#!/usr/bin/env python3

import argparse
import io
import sqlite3
import zipfile

from collections import defaultdict
from datetime import date

oio = {
    "SynonymTypeProperty": "synonym_type_property",
    "hasAlternativeId": "has_alternative_id",
    "hasBroadSynonym": "has_broad_synonym",
    "hasDbXref": "database_cross_reference",
    "hasExactSynonym": "has_exact_synonym",
    "hasOBOFormatVersion": "has_obo_format_version",
    "hasOBONamespace": "has_obo_namespace",
    "hasRelatedSynonym": "has_related_synonym",
    "hasScope": "has_scope",
    "hasSynonymType": "has_synonym_type",
}

exact_synonym = "oio:hasExactSynonym"
related_synonym = "oio:hasRelatedSynonym"
broad_synonym = "oio:hasBroadSynonym"

predicates = {
    "acronym": broad_synonym,
    "anamorph": related_synonym,
    "blast name": related_synonym,
    "common name": exact_synonym,
    "equivalent name": exact_synonym,
    "genbank acronym": broad_synonym,
    "genbank anamorph": related_synonym,
    "genbank common name": exact_synonym,
    "genbank synonym": related_synonym,
    "in-part": related_synonym,
    "misnomer": related_synonym,
    "misspelling": related_synonym,
    "synonym": related_synonym,
    "scientific name": exact_synonym,
    "teleomorph": related_synonym,
}

ranks = [
    "class",
    "cohort",
    "family",
    "forma",
    "genus",
    "infraclass",
    "infraorder",
    "kingdom",
    "order",
    "parvorder",
    "phylum",
    "section",
    "series",
    "species group",
    "species subgroup",
    "species",
    "subclass",
    "subcohort",
    "subfamily",
    "subgenus",
    "subkingdom",
    "suborder",
    "subphylum",
    "subsection",
    "subspecies",
    "subtribe",
    "superclass",
    "superfamily",
    "superkingdom",
    "superorder",
    "superphylum",
    "tribe",
    "varietas",
]

nodes_fields = [
    "tax_id",  # node id in GenBank taxonomy database
    "parent_tax_id",  # parent node id in GenBank taxonomy database
    "rank",  # rank of this node (superkingdom, kingdom, ...)
    "embl_code",  # locus-name prefix; not unique
    "division_id",  # see division.dmp file
    "inherited_div_flag",  # (1 or 0) 1 if node inherits division from parent
    "genetic_code_id",  # see gencode.dmp file
    "inherited_GC_flag",  # (1 or 0) 1 if node inherits genetic code from parent
    "mitochondrial_genetic_code_id",  # see gencode.dmp file
    "inherited_MGC_flag",  # (1 or 0) 1 if node inherits mitochondrial gencode from parent
    "GenBank_hidden_flag",  # (1 or 0) 1 if name is suppressed in GenBank entry lineage
    "hidden_subtree_root_flag",  # (1 or 0) 1 if this subtree has no sequence data yet
    "comments",  # free-text comments and citations
]


def escape_literal(text):
    return text.replace('"', '\\"')


def label_to_id(text):
    return text.replace(" ", "_").replace("-", "_")


def convert_synonyms(tax_id, synonyms):
    """Given a tax_id and list of synonyms,
    return a Turtle string asserting triples and OWL annotations on them."""
    pairs = []
    for synonym, unique, name_class in synonyms:
        if name_class in predicates:
            synonym = escape_literal(synonym)
            predicate = predicates[name_class]
            synonym_type = name_class.replace('genbank', 'GenBank')
            pairs.append([
                predicate,
                synonym,
                "xsd:string",
                f"""{{"oio:hasSynonymType":[{{"datatype":"xsd:string","meta":"owl:Axiom","object":"{synonym_type}"}}]}}"""
            ])
    return pairs


def convert_node(node, label, merged, synonyms, citations):
    """Given a node dictionary, a label string, and lists for merged, synonyms, and citations,
    return a Turtle string representing this tax_id."""
    tax_id = node["tax_id"]
    s = f"NCBITaxon:{tax_id}"
    pairs = [["rdf:type", "owl:Class"]]

    label = escape_literal(label)
    pairs.append(["rdfs:label", label, "xsd:string"])

    parent_tax_id = node["parent_tax_id"]
    if parent_tax_id and parent_tax_id != "" and parent_tax_id != tax_id:
        pairs.append(["rdfs:subClassOf", f"NCBITaxon:{parent_tax_id}"])

    rank = node["rank"]
    if rank and rank != "" and rank != "no rank":
        # if rank not in ranks:
        #     print(f"WARN Unrecognized rank '{rank}'")
        rank = label_to_id(rank)
        # WARN: This is a special case for backward compatibility
        if rank in ["species_group", "species_subgroup"]:
            pairs.append([
                "ncbitaxon:has_rank",
                f"<http://purl.obolibrary.org/obo/NCBITaxon#_{rank}>"
            ])
        else:
            pairs.append(["ncbitaxon:has_rank", f"NCBITaxon:{rank}"])

    gc_id = node["genetic_code_id"]
    if gc_id:
        pairs.append(["oio:hasDbXref", f"GC_ID:{gc_id}", "xsd:string"])

    for merge in merged:
        pairs.append(["oio:hasAlternativeId", f"NCBITaxon:{merge}", "xsd:string"])

    for pubmed_id in citations:
        pairs.append(["oio:hasDbXref", f"PMID:{pubmed_id}", "xsd:string"])

    pairs.append(["oio:hasOBONamespace", "ncbi_taxonomy", "xsd:string"])

    pairs += convert_synonyms(tax_id, synonyms)

    return [[s] + pair for pair in pairs]


def split_line(line):
    """Split a line from a .dmp file"""
    return [x.strip() for x in line.split("	|")]


datatypes = ["_IRI", "xsd:string"]


def insert_triples(con, triples):
    for triple in triples:
        if len(triple) < 3:
            raise Exception("Not enough slots in triple: " + str(triple))
        elif len(triple) == 3:
            triple.append("_IRI")
        elif len(triple) == 4 and triple[3] not in datatypes:
            raise Exception("Unexpected datatype in: " + str(triple))

        if len(triple) == 4:
            triple.append(None)

        if len(triple) > 5:
            raise Exception("Too many slots in triple: " + str(triple))

    cur = con.cursor()
    cur.executemany("""INSERT INTO ncbitaxon VALUES(1, 0, "obo:ncbitaxon", ?, ?, ?, ?, ?)""", triples)
    con.commit()


def convert(taxdmp_path, output_path, taxa=None):
    """Given the paths to the taxdmp.zip file and an output Turtle file,
    and an optional set of tax_id strings to extract,
    read from the taxdmp.zip file, collect annotations,
    convert nodes to Turtle strings,
    and write to the output file."""
    scientific_names = defaultdict(list)
    labels = {}
    synonyms = defaultdict(list)
    merged = defaultdict(list)
    citations = defaultdict(list)

    con = sqlite3.connect(output_path)
    cur = con.cursor()

    # WARN: These settings are fast, but not safe.
    cur.execute("PRAGMA journal_mode = OFF")
    cur.execute("PRAGMA synchronous = 0")
    cur.execute("PRAGMA cache_size = 1000000")
    cur.execute("PRAGMA locking_mode = EXCLUSIVE")
    cur.execute("PRAGMA temp_store = MEMORY")

    prefixes = [
        ["rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#"],
        ["rdfs", "http://www.w3.org/2000/01/rdf-schema#"],
        ["xsd", "http://www.w3.org/2001/XMLSchema#"],
        ["owl", "http://www.w3.org/2002/07/owl#"],
        ["obo", "http://purl.obolibrary.org/obo/"],
        ["oio", "http://www.geneontology.org/formats/oio#"],
        ["dct", "http://purl.org/dc/dct/"],
        ["ncbitaxon", "http://purl.obolibrary.org/obo/ncbitaxon#"],
        ["NCBITaxon", "http://purl.obolibrary.org/obo/NCBITaxon_"],
    ]
    # cur.execute("CREATE TABLE prefix(prefix TEXT PRIMARY KEY, base TEXT NOT NULL)")
    # cur.executemany("INSERT INTO prefix VALUES(?, ?)", prefixes)
    # con.commit()

    cur.execute("""CREATE TABLE ncbitaxon (
    assertion INT NOT NULL,
    retraction INT NOT NULL DEFAULT 0,
    graph TEXT NOT NULL,
    subject TEXT NOT NULL,
    predicate TEXT NOT NULL,
    object TEXT NOT NULL,
    datatype TEXT NOT NULL,
    annotation TEXT
)""")
    con.commit()

    isodate = date.today().isoformat()
    s = "<http://purl.obolibrary.org/obo/ncbitaxon.owl>"
    triples = [
        [s, "rdf:type", "owl:Ontology"],
        [s, "owl:versionIRI", f"<http://purl.obolibrary.org/obo/ncbitaxon/{isodate}/ncbitaxon.owl>"],
        [s, "dct:title", "NCBI organismal classification", "xsd:string"],
        [s, "dct:description", "An ontology representation of the NCBI organismal taxonomy", "xsd:string"],
        [s, "dct:license", "<https://creativecommons.org/publicdomain/zero/1.0/>"],
        [s, "rdfs:comment", "Built by https://github.com/obophenotype/ncbitaxon", "xsd:string"],
    ]
    insert_triples(con, triples)

    s = "obo:IAO_0000115"
    triples = [
        [s, "rdf:type", "owl:AnnotationProperty", "_IRI"],
        [s, "rdfs:label", "definition", "xsd:string"],
    ]
    insert_triples(con, triples)

    s = "ncbitaxon:has_rank"
    triples = [
        [s, "rdf:type", "owl:AnnotationProperty"],
        [s, "obo:IAO_0000115", "A metadata relation between a class and its taxonomic rank (eg species, family)", "xsd:string"],
        [s, "rdfs:label", "has_rank", "xsd:string"],
        [s, "rdfs:comment", "This is an abstract class for use with the NCBI taxonomy to name the depth of the node within the tree. The link between the node term and the rank is only visible if you are using an obo 1.3 aware browser/editor; otherwise this can be ignored", "xsd:string"],
        [s, "oio:hasOBONamespace", "ncbi_taxonomy", "xsd:string"],
    ]
    insert_triples(con, triples)

    for predicate, label in oio.items():
        s = f"oio:{predicate}"
        triples = [
            [s, "rdf:type", "owl:AnnotationProperty"],
            [s, "rdfs:label", label, "xsd:string"],
        ]
        insert_triples(con, triples)

    for label, parent in predicates.items():
        predicate = label_to_id(label)
        parent = parent.replace("oio", "oio")
        s = f"oio:{predicate}"
        triples = [
            [s, "rdf:type", "owl:AnnotationProperty"],
            [s, "rdfs:label", label, "xsd:string"],
            [s, "oio:hasScope", parent, "xsd:string"],
            [s, "rdfs:subPropertyOf", "oio:SynonymTypeProperty"],
        ]
        insert_triples(con, triples)

    with zipfile.ZipFile(taxdmp_path) as taxdmp:
        with taxdmp.open("names.dmp") as dmp:
            for line in io.TextIOWrapper(dmp):
                tax_id, name, unique, name_class, _ = split_line(line)
                if name_class == "scientific name":
                    labels[tax_id] = name
                    scientific_names[name].append([tax_id, unique])
                else:
                    synonyms[tax_id].append([name, unique, name_class])

        # use unique name only if there's a conflict
        for name, values in scientific_names.items():
            tax_ids = [x[0] for x in values]
            if len(tax_ids) > 1:
                uniques = [x[1] for x in values]
                if len(tax_ids) != len(set(uniques)):
                    print("WARN: Duplicate unique names", tax_ids, uniques)
                for tax_id, unique in values:
                    labels[tax_id] = unique
                    # Reason for the line below 
                    # issue #56: https://github.com/obophenotype/ncbitaxon/issues/56
                    if name != 'environmental samples':
                        synonyms[tax_id].append(
                            [name, unique, "scientific name"]
                        )

        with taxdmp.open("merged.dmp") as dmp:
            for line in io.TextIOWrapper(dmp):
                old_tax_id, new_tax_id, _ = split_line(line)
                merged[new_tax_id].append(old_tax_id)

        with taxdmp.open("citations.dmp") as dmp:
            for line in io.TextIOWrapper(dmp):
                (
                    cit_id,
                    cit_key,
                    pubmed_id,
                    medline_id,
                    url,
                    text,
                    tax_id_list,
                    _,
                ) = split_line(line)
                # WARN: the pubmed_id is always "0", we treat medline_id as pubmed_id
                if medline_id == "0":
                    continue
                for tax_id in tax_id_list.split():
                    if taxa and tax_id not in taxa:
                        continue
                    citations[tax_id].append(medline_id)

        with taxdmp.open("nodes.dmp") as dmp:
            for line in io.TextIOWrapper(dmp):
                node = {}
                fields = split_line(line)
                for i in range(0, min(len(fields), len(nodes_fields))):
                    node[nodes_fields[i]] = fields[i]
                tax_id = node["tax_id"]
                if taxa and tax_id not in taxa:
                    continue
                triples = convert_node(
                    node,
                    labels[tax_id],
                    merged[tax_id],
                    synonyms[tax_id],
                    citations[tax_id],
                )
                insert_triples(con, triples)

        # TODO: delnodes

    s = "<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>"
    triples = [
        [s, "rdf:type", "owl:Class"],
        [s, "rdfs:label", "taxonomic rank", "xsd:string"],
        [s, "rdfs:comment", "This is an abstract class for use with the NCBI taxonomy to name the depth of the node within the tree. The link between the node term and the rank is only visible if you are using an obo 1.3 aware browser/editor; otherwise this can be ignored.", "xsd:string"],
        [s, "oio:hasOBONamespace", "ncbi_taxonomy", "xsd:string"],
    ]
    insert_triples(con, triples)

    for label in ranks:
        rank = label_to_id(label)
        if rank in ["species_group", "species_subgroup"]:
            s = f"<http://purl.obolibrary.org/obo/NCBITaxon#_{rank}>"
        else:
            s = f"NCBITaxon:{rank}"
        triples = [
            [s, "rdf:type", "owl:Class"],
            [s, "rdfs:label", label, "xsd:string"],
            [s, "rdfs:subClassOf", "<http://purl.obolibrary.org/obo/NCBITaxon#_taxonomic_rank>"],
            [s, "oio:hasOBONamespace", "ncbi_taxonomy", "xsd:string"],
        ]
        insert_triples(con, triples)


def main():
    parser = argparse.ArgumentParser(
        description="Convert NCBI Taxonomy taxdmp.zip to Turtle format"
    )
    parser.add_argument("taxdmp", type=str, help="The taxdmp.zip file to read")
    parser.add_argument("taxa", type=str, nargs="?", help="A list of taxa to build")
    # TODO: upper, lower
    parser.add_argument("ldtab", type=str, help="The SQLite database to write")
    args = parser.parse_args()

    taxa = None
    if args.taxa:
        taxa = set()
        with open(args.taxa) as taxalist:
            for line in taxalist:
                taxa.add(line.split()[0])

    convert(args.taxdmp, args.ldtab, taxa)


if __name__ == "__main__":
    main()
