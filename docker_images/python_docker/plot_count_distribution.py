import argparse
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='norm_counts')
    parser.add_argument('-g', '--gene', dest='gene')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    tcga_type = os.path.basename(args.norm_counts).split('.')[0]
    nc = pd.read_table(args.norm_counts, index_col=0)
    counts = nc.loc[args.gene]
    fig, ax = plt.subplots(figsize=(10,10))
    sns.histplot(counts, bins=20, ax=ax)
    ax.set_title(f'{tcga_type}')
    plt.show()
    fig.savefig(f'nc.{args.gene}.{tcga_type}.png', bbox_inches='tight')
    plt.close()