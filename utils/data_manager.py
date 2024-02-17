from torch.utils.data import Dataset
import pandas as pd
import torch

class DNADataset(Dataset):
    def __init__(self, data, labels, tokenizer, max_length, kmer):
        self.data = data
        self.labels = labels
        self.tokenizer = tokenizer
        self.max_length = max_length
        self.include_meth = 'meth' in data.columns
        self.kmer = kmer
    def __len__(self):
        return len(self.data)
    def __getitem__(self, idx):
        sequence = self.data.iloc[idx]['seq']
        methylations = self.data.iloc[idx]['meth'] if self.include_meth else torch.Tensor([])
        label = self.labels[idx]
        # Tokenize sequence
        encoded_sequence = self.tokenizer(
            seq2kmer(sequence, int(self.kmer)),
            max_length=self.max_length,
            padding='max_length',
            truncation=True,
            return_tensors='pt'
        )
        input_ids = encoded_sequence["input_ids"].squeeze()
        attention_mask = encoded_sequence["attention_mask"].squeeze()
        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "label": label,
            "methylations": methylations,
            "sequence": sequence
        }

def seq2kmer(seq, k):
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers


def make_input_dataset(assembly, coordinate_df, label, meth_seq=None, window_size=128):
    res_df = pd.DataFrame()
    res_df['seq'] = coordinate_df.apply(lambda row: str(assembly[row['chr'].lower()][int(row['position']) - int(window_size//2): int(row['position']) + int(window_size//2)]), axis=1)
    if meth_seq != None: res_df['meth'] = coordinate_df.apply(lambda row: meth_seq[row['chr'].lower()][int(row['position']) - int(window_size//2): int(row['position']) + int(window_size//2)][:,0], axis=1)
    res_df['label'] = label
    return res_df

def load_df(dmc, not_dmc):
    return pd.read_csv(dmc, sep='\t', header=None, names=['chr', 'position']), pd.read_csv(not_dmc, sep='\t', header=None, names=['chr', 'position'])

# def load_df(KO):
#     if KO == 'TET' or KO == 'DNMT3' or KO == 'QKO' or KO == 'PKO':
#         if KO == 'TET':
#             KO_dmr = configs.TKO_DMRs_address
#             KO_notdmr = configs.TKO_notDMRs_address
#         elif KO == 'DNMT3':
#             KO_dmr = configs.DKO_DMRs_address
#             KO_notdmr = configs.DKO_notDMRs_address
#         elif KO == 'QKO':
#             KO_dmr = configs.QKO_DMRs_address
#             KO_notdmr = configs.QKO_notDMRs_address
#         elif KO == 'PKO':
#             KO_dmr = configs.PKO_DMRs_address
#             KO_notdmr = configs.PKO_notDMRs_address
#         dmrs = pd.read_csv(KO_dmr, header=None, sep='\t', names=['chr', 'position', 'meth1', 'coverage1', 'meth2', 'coverage2', 'p_value'])
#         not_dmrs = pd.read_csv(KO_notdmr, header=None, sep='\t', names=['chr', 'position', 'meth1', 'coverage1', 'meth2', 'coverage2'])
#     elif 'mouse' in KO:
#         adds = {'mouse_3aKO': (configs.mouse_3aKO_DMRs_address, configs.mouse_3aKO_notDMRs_address),
#                 'mouse_3bKO': (configs.mouse_3bKO_DMRs_address, configs.mouse_3bKO_notDMRs_address),
#                 'mouse_SI_TET23_KO': (configs.mouse_SI_TET23_DMRs_address, configs.mouse_SI_TET23_notDMRs_address)}
#         dmrs = pd.read_csv(adds[KO][0], header=None, sep='\t', names=['chr', 'position', 'meth1', 'coverage1', 'meth2', 'coverage2', 'p_value'])
#         not_dmrs = pd.read_csv(adds[KO][1], header=None, sep='\t', names=['chr', 'position', 'meth1', 'coverage1', 'meth2', 'coverage2'])
#     elif KO == 'fCpG':
#         fcpgs = dr.read_fcpgs(configs.fCpGs_address1)
#         epic2 = dr.read_epic2_conf(configs.epic_file)
#         dmrs = preprocess.make_fcpg_loc_df(fcpgs, epic2)
#         not_dmrs = epic2[~epic2['Name'].isin(fcpgs['ID'])][['Name', 'chr', 'position']].sample(n=len(dmrs))
#     return dmrs, not_dmrs



# class DNADataset(Dataset):
#     def __init__(self, sequences, labels, tokenizer, max_length, kmer):
#         self.sequences = sequences
#         self.labels = labels
#         self.tokenizer = tokenizer
#         self.max_length = max_length
#         self.kmer = kmer
#     def __len__(self):
#         return len(self.sequences)
#     def __getitem__(self, idx):
#         sequence = self.sequences[idx]
#         label = self.labels[idx]
#         # Tokenize sequence
#         input_ids = self.tokenizer(sequence)["input_ids"]
#         attention_mask = None
#         return {
#             "input_ids": input_ids,
#             "attention_mask": attention_mask,
#             "label": label,
#             "sequence": sequence
#         }
