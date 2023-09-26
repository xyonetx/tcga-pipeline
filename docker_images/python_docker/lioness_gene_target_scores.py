import argparse

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input')
    parser.add_argument('-t', '--testout', dest='testing_output')
    parser.add_argument('-s', '--scoreout', dest='scores_output')
    parser.add_argument('-a', '--ann', dest='annotations')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    ann_df = pd.read_table(args.annotations, index_col=0)

    low_samples = ann_df.loc[ann_df['expression_state'] == 'low'].index.values
    high_samples = ann_df.loc[ann_df['expression_state'] == 'high'].index.values

    # the lioness DF has a dual-index of (tf, gene)
    lioness_df = pd.read_csv(args.input, index_col=[0,1])
    lioness_df.index.names = ['tf','gene'] # just in case not named

    ts_df = pd.DataFrame(columns = lioness_df.columns)
    for gene, subdf in lioness_df.groupby(['gene']):
        target_scores = subdf.sum(axis=0)
        ts_df.loc[gene] = target_scores

    high_df = ts_df[high_samples]
    low_df = ts_df[low_samples]

    u, pvals = mannwhitneyu(low_df, high_df, axis=1)
    result = pd.DataFrame(pvals, index=ts_df.index, columns=['pval'])

    high_df_means = high_df.mean(axis=1)
    low_df_means = low_df.mean(axis=1)
    diff = high_df_means-low_df_means
    signs = np.sign(diff)
    result['rnk'] = -np.log(result['pval'])*signs

    ts_df.to_csv(args.scores_output, sep='\t')
    result['rnk'].sort_values(ascending=False).to_csv(args.testing_output, sep='\t', header=False)