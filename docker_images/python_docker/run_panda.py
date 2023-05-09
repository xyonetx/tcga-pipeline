import sys
import argparse
import pickle

import pandas as pd
import numpy as np
from netZooPy.panda.panda import Panda


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_matrix')
    parser.add_argument('-m', '--motif', dest='motifs')
    parser.add_argument('-p', '--ppi', dest='ppi')
    parser.add_argument('-a', '--ann', dest='annotations')
    parser.add_argument('-o', '--output', dest='output_pickle')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    ann_df = pd.read_table(args.annotations, index_col=0)
    samples = ann_df.index.values.tolist()

    # need to convert the sample names to their integer locations (columns) 
    # in the expression matrix we are going to feed to PANDA
    exp_mtx = pd.read_table(args.input_matrix, index_col=0)
    idx = np.where([x in samples for x in exp_mtx.columns])[0]
    
    panda_obj = Panda(args.input_matrix,
                args.motifs,
                args.ppi,
                keep_expression_matrix=True,
                modeProcess='legacy',
                remove_missing=True,
                with_header=True,
                save_memory=False)

    # Save PANDA object as pickle
    with open(args.output_pickle, "wb") as f:
        pickle.dump(panda_obj, f)