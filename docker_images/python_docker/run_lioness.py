import sys
import argparse
import pickle

import pandas as pd
import numpy as np
from netZooPy.lioness.lioness import Lioness


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--index', dest='index')
    parser.add_argument('-k', '--input_pkl', dest='input_pkl')
    parser.add_argument('-m', '--motif', dest='motifs')
    parser.add_argument('-p', '--ppi', dest='ppi')
    parser.add_argument('-n', '--num_cores', dest='num_cores')
    parser.add_argument('-a', '--ann', dest='annotations')
    parser.add_argument('-o', '--output', dest='output')
    return parser.parse_args()


def load_panda_obj(pickle_name):
    with open(pickle_name, 'rb') as f:
        panda_obj = pickle.load(f)
        return panda_obj


if __name__ == '__main__':
    args = parse_args()

    panda_obj = load_panda_obj(args.input_pkl)

    # read the split file (which simply contains an integer
    # indicating the split index) and parse as an integer:
    idx = int(open(args.index).read().strip())

    # parse the annotation file and subset to keep only those
    # samples in the current shard:
    ann_df = pd.read_table(args.annotations, index_col=0)
    ann_df = ann_df.loc[ann_df.split_id == idx]

    selected_samples = ann_df.index

    exp_samples = panda_obj.expression_samples
    locations = [exp_samples.get_loc(k) for k in selected_samples]

    output_file = f'{args.output}.{idx}.csv'

    Lioness(panda_obj,
            ncores=args.num_cores,
            subset_numbers=locations,
            export_filename=output_file)