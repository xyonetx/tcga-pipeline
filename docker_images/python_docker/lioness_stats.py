import argparse

import pandas as pd
from scipy.stats import mannwhitneyu


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input')
    parser.add_argument('-o', '--output', dest='output')
    parser.add_argument('-a', '--ann', dest='annotations')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    ann_df = pd.read_table(args.annotations, index_col=0)

    low_samples = ann_df.loc[ann_df['expression_state'] == 'low'].index.values
    high_samples = ann_df.loc[ann_df['expression_state'] == 'high'].index.values

    lioness_df = pd.read_csv(args.input, index_col=[0,1])

    high_df = lioness_df[high_samples]
    low_df = lioness_df[low_samples]

    u, pvals = mannwhitneyu(low_df, high_df, axis=1)
    result = pd.DataFrame(pvals, index=lioness_df.index)
    result.to_csv(args.output)