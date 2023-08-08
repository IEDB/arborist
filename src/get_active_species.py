import csv
import sqlite3

from argparse import ArgumentParser, FileType


def get_ancestors(cur, term_id):
    """Return a set of ancestors for a given term ID, all the way to the top-level (below owl:Thing)

    :param cur: database Cursor object to query
    :param term_id: term ID to get the ancestors of"""
    cur.execute(
        f"""WITH RECURSIVE ancestors(node) AS (
                VALUES ('{term_id}')
                UNION
                 SELECT object AS node
                FROM organism_tree
                WHERE predicate IN ('rdfs:subClassOf', 'rdfs:subPropertyOf')
                  AND object = '{term_id}'
                UNION
                SELECT object AS node
                FROM organism_tree, ancestors
                WHERE ancestors.node = organism_tree.subject
                  AND organism_tree.predicate IN ('rdfs:subClassOf', 'rdfs:subPropertyOf')
                  AND organism_tree.object NOT LIKE '_:%'
            )
            SELECT * FROM ancestors""",
    )
    return set([x[0] for x in cur.fetchall()])


def get_descendants(cur, term_id):
    """Return a set of descendants for a given term ID."""
    cur.execute(
        f"""WITH RECURSIVE descendants(node) AS (
                VALUES ('{term_id}')
                UNION
                 SELECT subject AS node
                FROM organism_tree
                WHERE predicate IN ('rdfs:subClassOf', 'rdfs:subPropertyOf')
                  AND subject = '{term_id}'
                UNION
                SELECT subject AS node
                FROM organism_tree, descendants
                WHERE descendants.node = organism_tree.object
                  AND organism_tree.predicate IN ('rdfs:subClassOf', 'rdfs:subPropertyOf')
            )
            SELECT * FROM descendants""",
    )
    return set([x[0] for x in cur.fetchall()])


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("counts", type=FileType("r"))
    parser.add_argument("output", type=FileType("w"))
    args = parser.parse_args()

    active_tax_ids = []
    reader = csv.DictReader(args.counts, delimiter="\t")
    for row in reader:
        active_tax_ids.append(row["source_organism_org_id"])

    rows = []
    with sqlite3.connect(args.db) as conn:
        cur = conn.cursor()
        cur.execute(
            """SELECT subject FROM organism_tree
            WHERE predicate = 'iedb-taxon:level' and object = 'species' ORDER BY subject"""
        )
        for res in cur.fetchall():
            curie = res[0]
            tax_id = int(curie.split(":")[1])
            if tax_id == 694009:
                continue

            # Get label
            cur.execute(
                "SELECT object FROM organism_tree WHERE predicate = 'rdfs:label' AND subject = ?",
                (curie,),
            )
            res = cur.fetchone()
            if not res:
                print("ERROR: Missing label for " + curie)
                continue
            label = res[0]
            species_key = (
                str(tax_id)
                + "-"
                + label.replace(" ", "-")
                .replace("(", "")
                .replace(")", "")
                .replace(".", "")
                .replace("'", "")
                .replace('"', "")
                .replace("/", "_")
            )

            # Get the active taxa
            descendants = get_descendants(cur, curie)
            active_taxa = [
                int(x.split(":")[1]) for x in descendants if x.split(":")[1] in active_tax_ids
            ]
            if tax_id == 10002316:
                # Add parent SARS Coronavirus to SARS-CoV1
                active_taxa.append(694009)
            if not active_taxa:
                # No taxa with epitopes, this is not active
                continue
            active_taxa.sort()

            # Get the group
            ancestors = get_ancestors(cur, curie)
            group = None
            if "OBI:0100026" not in ancestors:
                group = "other"
            elif "NCBITaxon:2157" in ancestors:
                group = "archeobacterium"
            elif "NCBITaxon:2" in ancestors:
                group = "bacterium"
            elif "NCBITaxon:10239" in ancestors:
                group = "virus"
            elif "NCBITaxon:2759" in ancestors:
                if "NCBITaxon:58024" in ancestors:
                    group = "plant"
                elif "NCBITaxon:7742" in ancestors:
                    group = "vertebrate"
                else:
                    group = "other-eukaryote"
            if not group:
                print("ERROR: Unable to find a group for " + curie)

            rows.append(
                {
                    "Species Key": species_key,
                    "Species ID": tax_id,
                    "Species Label": label,
                    "Active Taxa": ", ".join([str(x) for x in active_taxa]),
                    "Group": group,
                }
            )

    rows = sorted(rows, key=lambda k: k["Species ID"])

    writer = csv.DictWriter(
        args.output,
        fieldnames=["Species Key", "Species ID", "Species Label", "Active Taxa", "Group"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
