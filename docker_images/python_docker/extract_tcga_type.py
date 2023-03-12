'''
This script extracts out the count matrix for a
given TCGA type (e.g. TCGA-BRCA)
'''

import argparse
import sys

import pandas as pd


def parse_cl_args():
    '''
    Parse the commandline args and return a namespace
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--hdf5', dest='input_hdf5', required=True)
    parser.add_argument('-t', '--type', dest='tcga_type', required=True)
    parser.add_argument('-o', '--output', dest='output_file', required=True)
    return parser.parse_args()


def subset_hdf5(input_hdf5, tcga_type):
    '''
    Extract out a single TCGA type from the HDF5 file.
    Returns a pandas dataframe
    '''
    tcga_type = tcga_type.replace('-', '_').lower()
    group_id = f'{tcga_type}/ds'
    with pd.HDFStore(input_hdf5, 'r') as fin:
        try:
            return fin.get(group_id)
        except KeyError as ex:
            sys.stderr.write(f'Cannot find key {ex}')
            sys.exit(1)


if __name__ == '__main__':
    args = parse_cl_args()
    df = subset_hdf5(args.input_hdf5, args.tcga_type)
    df.to_csv(args.output_file, sep='\t')
