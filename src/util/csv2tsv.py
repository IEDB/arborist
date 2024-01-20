#!/usr/bin/env python3

import argparse
import csv


def convert(input_path, output_path):
    input = open(input_path, 'r')
    output = open(output_path, 'w')
    rows = csv.reader(input)
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerows(rows)
    input.close()
    output.close()


def main():
    parser = argparse.ArgumentParser(description='Convert CSV to TSV')
    parser.add_argument('input', type=str, help='The input CSV fila path')
    parser.add_argument('output', type=str, help='The output TSV file path')
    args = parser.parse_args()

    convert(args.input, args.output)


if __name__ == '__main__':
    main()
