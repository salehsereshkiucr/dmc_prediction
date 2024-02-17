import pandas as pd
import numpy as np
from os import path
from utils import data_reader

def make_fcpg_loc_df(fcpgs, epic2):
    merged_df = fcpgs.merge(epic2[['Name', 'chr', 'position']], left_on='ID', right_on='Name', how='left')
    return merged_df[['Name', 'chr', 'position']].dropna()

def make_methseq_dic(accession, methylations, sequences, coverage_thrshld):
    fn = './dump_files/{}_meth_seq_{}.pkl'.format(accession, str(coverage_thrshld))
    if path.exists(fn):
        return data_reader.load_dic(fn)
    methylations['meth'] = methylations['meth'].fillna(0)
    methylations.loc[(methylations.meth == 0), 'meth'] = 0
    methylations.loc[(methylations.coverage < coverage_thrshld), 'meth'] = 0
    meth_seq = {}
    count = 0
    for chr in sequences.keys():
        meths = np.zeros(len(sequences[chr]))
        meth_subset = methylations[methylations['chr'] == chr]
        meths[[meth_subset['position']]] = meth_subset['meth']
        meth_seq[chr.lower()] = np.expand_dims(meths, axis=1)
        count += 1
        print(count, len(sequences.keys()))
    data_reader.save_dic(fn, meth_seq)
    return meth_seq

def clean_grch38_seq_dic_chros(seq_dic):
    chro_df = pd.read_csv('./GRCH38_chromosoms.txt', sep='\t', header=None)
    chro_dic = {row[3]: 'chr'+row[0].split()[-1] for index, row in chro_df.iterrows()}
    return {chro_dic[chro]: seq_dic[chro] for chro in chro_dic.keys()}

def fix_annot(annot_df):
    chro_df = pd.read_csv('./GRCH38_chromosoms.txt', sep='\t', header=None)
    replace_dic = {row[3]: 'chr'+row[0].split()[-1] for index, row in chro_df.iterrows()}
    annot_df['chr'] = annot_df['chr'].replace(replace_dic, inplace=False)
    return annot_df

def convert_assembely_to_onehot(organism_name, sequences, from_file=False):
    fn = './dump_files/' + organism_name + '_sequences_onehot.pkl'
    if from_file and path.exists(fn):
        return data_reader.load_dic(fn)
    one_hots = {}
    for key in sequences.keys():
        one_hots[key] = convert_seq_to_onehot(str(sequences[key]))
        print('done', key)
    data_reader.save_dic(fn, one_hots)
    return one_hots

# def convert_seq_to_onehot(sequence):
#     sequence = str(sequence)
#     extra_letters = {'M', 'B', 'S', 'W', 'R', 'Y', 'N', 'K'}
#     for ll in extra_letters:
#         sequence = sequence.replace(ll, 'N')
#     characters = ['A', 'C', 'G', 'T', 'N']
#     mlb = MultiLabelBinarizer(classes=characters)
#     one_hot_encoding = mlb.fit_transform([[char] for char in sequence])
#     return one_hot_encoding

def convert_seq_to_onehot(sequence):
    sequence = str(sequence)
    extra_letters = {'M', 'B', 'S', 'W', 'R', 'Y', 'N', 'K'}
    for ll in extra_letters:
        sequence = sequence.replace(ll, 'N')
    characters = ['A', 'C', 'G', 'T', 'N']
    seq_ser = pd.Series(list(sequence))
    res_series = []
    for ch in characters:
        new_ser = pd.Series(0, index=seq_ser.index)
        new_ser[seq_ser == ch] = 1
        res_series.append(new_ser)
    return pd.concat(res_series, axis=1).values.astype(bool)


def make_annotseq_dic(annot_df, annot_types, sequences, strand_spec=False):
    annot_seqs = {}
    strand_num = 2 if strand_spec else 1
    for chr in sequences.keys():
        annot_seqs[chr] = np.zeros((len(sequences[chr]), strand_num*len(annot_types)))
        for at_idx, at in enumerate(annot_types):
            annot_subset_df = annot_df[annot_df.type == at]
            annot_seq_p = np.zeros((len(sequences[chr]), 1), dtype='short')
            annot_seq_n = np.zeros((len(sequences[chr]), 1), dtype='short')
            annot_df_chr_subset = annot_subset_df[annot_subset_df['chr'] == chr]
            for index, row in annot_df_chr_subset.iterrows():
                if row['strand'] == '+' or not strand_spec:
                    annot_seq_p[int(row['start'] - 1): int(row['end'] - 1)] = 1
                else:
                    annot_seq_n[int(row['start'] - 1): int(row['end'] - 1)] = 1
            annot_seqs[chr][:, at_idx*strand_num] = annot_seq_p[:, 0]
            if strand_spec:
                annot_seqs[chr][:, at_idx*strand_num + 1] = annot_seq_n[:, 0]
        print(chr)
    return annot_seqs


# sequences = data_reader.readfasta(cnfg['assembly'])
#
# one_hots = {}
# for key in sequences.keys():
#     one_hots[key.lower()] = convert_seq_to_onehot(str(sequences[key]))
#     print('done', key)
