import argparse
import glob
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', dest='directory')
    parser.add_argument('-o', '--output', dest='output')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    matching_files = glob.glob(f'{args.directory}/lioness.*.csv')

    # the files are TF x gene. Performing a column sum gives the 
    # targeting score for each gene
    master_df = pd.DataFrame()
    for f in matching_files:
        sample_id = os.path.basename(f).split('.')[1]
        df = pd.read_csv(f, index_col=0)
        s = df.sum(axis=0)
        s.name = sample_id
        master_df = pd.concat([master_df, s], axis=1)

    # master_df is a genes x samples matrix
    master_df.to_csv(args.output, sep='\t', float_format='%.5f')