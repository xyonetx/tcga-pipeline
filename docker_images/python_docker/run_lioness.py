import sys
import argparse
import pickle

import pandas as pd
import numpy as np
from netZooPy.lioness.lioness import Lioness


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_matrix')
    parser.add_argument('-m', '--motif', dest='motifs')
    parser.add_argument('-p', '--ppi', dest='ppi')
    parser.add_argument('-n', '--num_cores', dest='num_cores')
    parser.add_argument('-a', '--ann', dest='annotations')
    parser.add_argument('-k', '--pickle', dest='panda_pickle')
    return parser.parse_args()


def load_panda_obj(pickle_name):
    with open(pickle_name, 'rb') as f:
        panda_obj = pickle.load(f)
        return panda_obj


if __name__ == '__main__':
    args = parse_args()

    ann_df = pd.read_table(args.annotations, index_col=0)
    samples = ann_df.index.values.tolist()

    # need to convert the sample names to their integer locations (columns) 
    # in the expression matrix we are going to feed to PANDA
    exp_mtx = pd.read_table(args.input_matrix, index_col=0)
    idx = np.where([x in samples for x in exp_mtx.columns])[0]

    panda_obj = load_panda_obj(args.panda_pickle)

    Lioness(panda_obj,
            ncores=args.num_cores,
            save_single=True,
            ignore_final=True,
            subset_numbers=idx,
            precision='single',
            save_fmt='csv')
