import argparse
import pandas as pd
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', dest='output')
    parser.add_argument('input_files', nargs='+')
    return parser.parse_args()


def merge_lioness_scatter(input_files):
    '''
    Merge the slices into a single matrix
    '''
    # this uses a multiIndex which combines the first two columns
    # containing the transcription factor and gene
    df = pd.read_csv(input_files[0], index_col=[0,1])

    for fname in input_files[1:]:
        other = pd.read_csv(fname, index_col=[0,1])
        df = pd.concat([df, other], axis=1)
    return df


if __name__ == "__main__":
    args = parse_args()
    df = merge_lioness_scatter(args.input_files)
    df.to_csv(args.output, float_format='%.5f')