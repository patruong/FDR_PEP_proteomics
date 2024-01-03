import pandas as pd
import argparse

def fix_target_decoy(input_file, output_file):
    df = pd.read_csv(input_file, sep = "\t")
    map_targetdecoy = lambda x:1 if x == "target" else -1
    df["target/decoy"] = df["target/decoy"].map(map_targetdecoy)
    df.to_csv(output_file, sep = "\t", index = False)

if __name__ == "__main__":
    print("Replacing target with 1 and decoy with -1.")
    parser = argparse.ArgumentParser(description='Fix target/decoy values in a file.')
    parser.add_argument('input_file', type=str, help='The input file')
    parser.add_argument('output_file', type=str, help='The output file')
    args = parser.parse_args()

    fix_target_decoy(args.input_file, args.output_file)
    print("Output:"+ args.output_file)