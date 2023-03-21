import sys
import argparse

from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_matrix')
    parser.add_argument('-m', '--motif', dest='motifs')
    parser.add_argument('-p', '--ppi', dest='ppi')
    parser.add_argument('-n', '--num_cores', dest='num_cores')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    panda_obj = Panda(args.input_matrix,
                args.motifs,
                args.ppi,
                keep_expression_matrix=True,
                with_header=True,
                save_memory=False)

    Lioness(panda_obj,
            ncores=args.num_cores,
            save_single=True,
            ignore_final=True,
            save_fmt='csv')
