import argparse
import glob
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', dest='directory')
    parser.add_argument('-g', '--genes', dest='genes_output')
    parser.add_argument('-t', '--tf', dest='tf_output')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    matching_files = glob.glob(f'{args.directory}/lioness.*.csv')

    # the files are TF x gene. Performing a column sum gives the 
    # targeting score for each gene. Row-sum gives targeting score
    # for tf
    gene_df = pd.DataFrame()
    tf_df = pd.DataFrame()
    for f in matching_files:
        sample_id = os.path.basename(f).split('.')[1]
        df = pd.read_csv(f, index_col=0)
        colsum = df.sum(axis=0)
        rowsum = df.sum(axis=1)
        colsum.name = sample_id
        rowsum.name = sample_id
        gene_df = pd.concat([gene_df, colsum], axis=1)
        tf_df = pd.concat([tf_df, rowsum], axis=1)

    # gene_df is a genes x samples matrix
    gene_df.to_csv(args.genes_output, sep='\t', float_format='%.5f')

    # tf_df is a genes x samples matrix
    tf_df.to_csv(args.tf_output, sep='\t', float_format='%.5f')