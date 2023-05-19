import argparse
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.compare import compare_survival


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('annotations', nargs='+')
    return parser.parse_args()



def make_kaplan_meier_plot(ann, tcga_type):
    '''
    Create a survival plot for the two cohorts and
    calculate the p-value from log-rank test.

    `ann` is a dataframe containing 'days_to_death', 'vital_status'
        and 'expression_state' columns
    '''
    print(f'here {tcga_type}')
    ann['status'] = ann['vital_status'].apply(lambda x: False if x=='Alive' else True)
    dtd_max = ann['days_to_death'].max()
    ann['time_exit'] = ann['days_to_death'].apply(lambda x: dtd_max if pd.isna(x) else x)

    l = []
    for item in ann[['status','time_exit']].iterrows():
        l.append(
            (item[1]['status'], item[1]['time_exit'])
        )
    l = np.array(l, dtype=[('status', 'bool'), ('time', '<f4')])
    chi_sqd, pval = compare_survival(l, ann['expression_state'].values)
   
    fig, ax = plt.subplots(figsize=(10,10))
    for g, subdf in ann.groupby('expression_state'):
        time, surv_prob = kaplan_meier_estimator(subdf['status'], subdf['time_exit'])
        ax.step(time, surv_prob, where='post', label=f'MUC1 {g}')
    ax.legend()
    ax.set_xlabel('Days to Death', fontsize=20)
    ax.set_ylabel('$\hat{S}$', fontsize=20)
    ax.set_title(f'{tcga_type} (pval={pval:.3f})')
    plt.show()
    fig.savefig(f'kaplan_meier.{tcga_type}.png', bbox_inches='tight')
    plt.close()

    
if __name__ == '__main__':
    args = parse_args()
    total_ann = pd.DataFrame()
    for ann_file in args.annotations:
        ann_df = pd.read_table(ann_file, index_col=0)
        tcga_type = os.path.basename(ann_file.split('.')[0])
        make_kaplan_meier_plot(ann_df, tcga_type)
        total_ann = pd.concat([total_ann, ann_df], axis=0)
    
    make_kaplan_meier_plot(total_ann, 'Overall')