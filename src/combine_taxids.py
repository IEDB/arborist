from argparse import ArgumentParser, FileType


def main():
    parser = ArgumentParser()
    parser.add_argument("counts", type=FileType("r"))
    parser.add_argument("ncbi_include", type=FileType("r"))
    parser.add_argument("output", type=FileType("w"))
    args = parser.parse_args()

    counts = [x.strip() for x in args.counts.readlines()]
    included_ids = [x.split("\t")[0] for x in counts if x != "source_organism_org_id	count"]
    ncbi_include = []
    for line in args.ncbi_include.readlines():
        line = line.strip()
        if line == "taxon_id":
            continue
        if line not in included_ids:
            ncbi_include.append(line + "\t0")
    counts.extend(ncbi_include)
    args.output.writelines("\n".join(counts))


if __name__ == "__main__":
    main()
