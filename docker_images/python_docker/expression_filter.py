import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_matrix')
    parser.add_argument('-o', '--motif', dest='output')
    parser.add_argument('-n', '--ppi', dest='num_genes', type=int)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    exprs_df = pd.read_table(args.input_matrix, sep='\t', index_col=0)

    # take the row-wise mean, sort, and take the top N values
    exprs_df = exprs_df.loc[
        exprs_df.mean(axis=1).sort_values(ascending=False).index[:args.num_genes]
    ]
    exprs_df.to_csv(args.output, sep='\t')