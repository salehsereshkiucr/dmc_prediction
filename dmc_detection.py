import sys
import os

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

import pandas as pd
import argparse
from scipy.stats import fisher_exact

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--control', help='address of a bed file with the methylation levels of cytosines in control set', required=True)
parser.add_argument('-t', '--test', help='address of a bed file with the methylation levels of cytosines in test set', required=True)
parser.add_argument('-ct', '--coverage_threshold', help='The minimum coverage for both test and control to call a cytosine dmc (not_dmc)', required=False, default=10)
parser.add_argument('-md', '--minimum_difference', help='The minimum difference between methylation levels of control and test to call a cytosine as dmc', required=False, default=0.6)

#python dmc_detection.py -t ./data/test_mC_report.bed -c ./data/control_mC_report.bed


args = parser.parse_args()

control_add = args.control
test_add = args.test
coverage_threshold = int(args.coverage_threshold)
min_diff = float(args.minimum_difference)

# control_add = './data/control_mC_report.bed'
# test_add = './data/test_mC_report.bed'
# coverage_threshold = 10
# min_diff = 0.6

def calculate_p_value(row):
    contingency_table = [
        [round(row['meth_x'] * row['coverage_x']), round((1 - row['meth_x']) * row['coverage_x'])],
        [round(row['meth_y'] * row['coverage_y']), round((1 - row['meth_y']) * row['coverage_y'])]
    ]
    _, p_value = fisher_exact(contingency_table)
    return p_value


df1, df2 = pd.read_csv(control_add, sep='\t', names=['chr', 'start', 'end', 'meth', 'coverage']), pd.read_csv(test_add, sep='\t', names=['chr', 'start', 'end', 'meth', 'coverage'])
df1, df2 = df1[df1.coverage > coverage_threshold], df2[df2.coverage > coverage_threshold]
df1, df2 = df1.drop('end', axis=1), df2.drop('end', axis=1)
merged_df = pd.merge(df1, df2, on=['chr', 'start']) #20 sec run
print('Found the common Cytosines, {} from {} of control and {} of test'.format(len(merged_df), len(df1), len(df2)))
dmrs = merged_df[abs(merged_df.meth_x - merged_df.meth_y) > min_diff]
dmrs.loc[:, 'p_value'] = dmrs.apply(calculate_p_value, axis=1)
subset_df = pd.merge(merged_df, dmrs, on=['chr', 'start'], how='left', indicator=True)
not_dmrs = subset_df[subset_df['_merge'] == 'left_only'].drop(columns=['_merge']).iloc[:, :6]
not_dmrs.columns = dmrs.columns[:6]
dmrs, not_dmrs = dmrs[['chr', 'start']], not_dmrs[['chr', 'start']]
not_dmrs.to_csv('./data/dmcs.txt', sep='\t', index=False, header=False)
dmrs.to_csv('./data/not_dmcs.txt', sep='\t', index=False, header=False)

print('\n\n\n\n\n')
print('Results were saved in ./data/dmcs.txt and ./data/not_dmcs.txt')
print('DONE!')

