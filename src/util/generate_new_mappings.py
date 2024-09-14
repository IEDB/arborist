import polars as pl
import sys

def generate_new_mappings(latest_file, previous_file, output_file):
    latest_df = pl.read_csv(latest_file, separator='\t')
    previous_df = pl.read_csv(previous_file, separator='\t')

    new_rows = latest_df.join(previous_df, on='epitope_id', how='anti')

    new_rows.write_csv(output_file, separator='\t')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <latest_file> <previous_file> <output_file>")
        sys.exit(1)

    latest_file = sys.argv[1]
    previous_file = sys.argv[2]
    output_file = sys.argv[3]

    generate_new_mappings(latest_file, previous_file, output_file)
