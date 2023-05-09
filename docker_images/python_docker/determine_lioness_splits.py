import sys
import argparse

import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--ann', dest='annotations')
    parser.add_argument('-o', '--output', dest='output')
    parser.add_argument('-n', '--nmax', dest='nmax', type=int)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    ann_df = pd.read_table(args.annotations, index_col=0)
    N = ann_df.shape[0]

    ann_df['split_id'] = np.arange(N) // args.nmax
    ann_df.to_csv(args.output, sep='\t')

    split_num = int(np.ceil(N/args.nmax))
    for i in range(split_num):
        with open(f'split_{i}.txt', 'w') as fout:
            fout.write(str(i))