import pandas as pd
from Bio import SeqIO
import pickle

def read_fcpgs(address):
    fcpgs = pd.read_csv(address)
    fcpgs = fcpgs.rename(columns={'Unnamed: 0': 'ID'})
    return fcpgs

def read_epic2_conf(address):
    epic2 = pd.read_csv(address)
    epic2.columns = epic2.iloc[6]
    epic2 = epic2.iloc[7:]
    epic2 = epic2[epic2['CHR'].isin(['chr'+str(i+1) for i in range(23)])]
    epic2 = epic2.rename(columns={'CHR': 'chr', 'MAPINFO': 'position'})
    return epic2

def read_annot(address, chromosomes = None):
    annot_df = pd.read_table(address, sep='\t', comment='#')
    annot_df.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    if chromosomes != None:
        annot_chrs = annot_df.chr.unique()
        for chr in annot_chrs:
            if chr not in chromosomes:
                annot_df = annot_df[annot_df['chr'] != chr]
    annot_df['chr'] = annot_df['chr'].apply(str.lower)
    return annot_df

def read_wgbs_bed(address):
    df = pd.read_csv(address, sep='\t', header=None, comment='#')
    if df.shape[1] != 5:
        df = df.iloc[:, :5]
    df.columns = ['chr', 'start', 'end', 'meth', 'coverage']
    return df


def readfasta(address):
    recs = SeqIO.parse(address, "fasta")
    sequences = {}
    for chro in recs:
        sequences[chro.id.lower()] = chro.seq
    for i in sequences.keys():
        sequences[i] = sequences[i].upper()
    return sequences

def write_fasta_file(dictionary, file_path):
    with open(file_path, 'w') as file:
        for header, body in dictionary.items():
            file.write('>' + header + '\n')
            file.write(str(body) + '\n')

def write_fasta(sequences, file_name):
    with open(file_name, 'w') as file:
        for i, seq in enumerate(sequences, 1):
            file.write(f">{i}\n{seq}\n")

def save_dic(file_name, dic):
    f = open(file_name, "wb")
    pickle.dump(dic, f)
    f.close()
    print(file_name + ' saved as pickle ')

def load_dic(file_name):
    f = open(file_name, "rb")
    res = pickle.load(f)
    f.close()
    print(file_name + ' loaded from a pickle ')
    return res

