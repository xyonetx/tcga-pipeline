'''
This script extracts out the high and low quantiles
of a desired gene and writes them to separate files.
'''

import argparse
import sys

import pandas as pd


def parse_cl_args():
    '''
    Parse the commandline args and return a namespace
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input', dest='input_file', required=True)
    parser.add_argument('-g', '--gene', dest='gene', required=True)
    parser.add_argument('-pl', '--low', dest='low_quantile',
                        type=float, required=True)
    parser.add_argument('-ph', '--high', dest='high_quantile',
                        type=float, required=True)
    parser.add_argument(
        '--low-output', dest='low_q_output_file', required=True)
    parser.add_argument(
        '--high-output', dest='high_q_output_file', required=True)
    return parser.parse_args()


def get_gene_expressions(exp_df, gene):
    '''
    Return the expressions for a given gene.
    The returned value is a pandas Series instance
    '''
    try:
        return exp_df.loc[gene]
    except KeyError as ex:
        sys.stderr.write(f'Could not find gene: {ex}')
        sys.exit(1)


def split_data(gene_series, low_quantile, high_quantile):
    '''
    Given a pandas Series of expressions and the low/high quantiles,
    return the samples corresponding to the low and high
    expression groups. Returns a 2-item tuple where each item
    is a list of sample names (column IDs)
    '''
    quantiles = gene_series.quantile([low_quantile, high_quantile])
    low_expression_samples = gene_series.index[
        gene_series < quantiles[low_quantile]
    ]
    high_expression_samples = gene_series.index[
        gene_series > quantiles[high_quantile]
    ]
    return low_expression_samples, high_expression_samples


if __name__ == '__main__':
    args = parse_cl_args()
    df = pd.read_table(args.input_file, index_col=0, sep='\t')
    gene_expressions = get_gene_expressions(df, args.gene)
    low_exp_samples, high_exp_samples = split_data(
        gene_expressions,
        args.low_quantile,
        args.high_quantile
    )
    low_exp_df = df[low_exp_samples]
    high_exp_df = df[high_exp_samples]

    low_exp_df.to_csv(args.low_q_output_file, sep='\t')
    high_exp_df.to_csv(args.high_q_output_file, sep='\t')
