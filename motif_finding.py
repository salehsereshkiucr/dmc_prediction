import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)
import utils.data_reader as data_reader
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import argparse
from transformers import BertTokenizer, BertModel
import utils.model_manager as mdlmngr
from utils.data_manager import DNADataset
import utils.motif_utils as mu
import subprocess


parser = argparse.ArgumentParser()

parser.add_argument('-dmc_seq', '--dmc_seq', help='address of context sequences of differentially methylated cytosines in fasta format', required=True)
parser.add_argument('-ndmc_seq', '--ndmc_seq', help='address of context sequences of not differentially methylated cytosines in fasta format', required=True)
parser.add_argument('-mn', '--model_name', help='the address to load model', required=False, default='./models/trained_clf')
parser.add_argument('-ws', '--window_size', help='window_size', required=False, default=512)
parser.add_argument('-kmer', '--kmer', help='kmer', required=False, default=6)

args = parser.parse_args()

dmc_seq = args.dmc_seq
ndmc_seq = args.ndmc_seq
model_address = args.model_name
kmer = int(args.kmer)
window_size = int(args.window_size)

#python motif_finding.py -dmc_seq ./data/motif_sample_pos.fa -ndmc_seq ./data/motif_sample_neg.fa

dmc_seq = './data/motif_sample_pos.fa'
ndmc_seq = './data/motif_sample_neg.fa'
model_address = './models/trained_clf'
kmer = 6
window_size = 512

pos_seqs = list(data_reader.readfasta(dmc_seq).values())
neg_seqs = list(data_reader.readfasta(ndmc_seq).values())

clf_model = mdlmngr.load_clf_model(model_address)
#clf_model = mdlmngr.load_clf_model("./dump_files/pretrained_models_/" + "kmer_"+str(kmer) + "_ws_" + str(window_size) + "_KO_" + str(train_KO))
tokenizer = BertTokenizer.from_pretrained("zhihan1996/DNA_bert_"+str(kmer), trust_remote_code=True)
#clf_model = mdlmngr.load_clf_model("./dump_files/pretrained_models_/" + "kmer_"+str(kmer) + "_ws_" + str(window_size) + "_KO_" + str(KO))
print('model loaded from '+ model_address)


def model_eval(X, Y, model, tokenizer, window_size, kmer):
    batch_size = 2
    labels = torch.tensor(Y, dtype=torch.float32)
    dataset = DNADataset(X, labels, tokenizer, window_size, kmer)
    data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model.to(device)
    model.eval()
    sequences = np.empty(len(dataset), dtype=object)
    labels, scores = np.zeros(len(dataset)), np.zeros((len(dataset), window_size))
    idx = 0
    with torch.no_grad():
        for batch in data_loader:
            input_ids = batch["input_ids"].to(device)
            attention_mask = batch["attention_mask"].to(device)
            o = model(input_ids=input_ids, attention_mask=attention_mask, output_attentions=True)
            sequences[idx:idx+len(input_ids)] = batch["sequence"]
            labels[idx:idx+len(input_ids)] = batch["label"].argmax(dim=1).numpy()
            scores[idx:idx+len(input_ids)] = get_avg_att_scores(o[-1][-1].cpu().numpy(), kmer)
            idx += len(input_ids)
    return sequences, labels.astype(int), scores

def get_avg_att_scores(attention_scores, kmer):
    scores = np.zeros([attention_scores.shape[0], attention_scores.shape[-1]])
    for index, attention_score in enumerate(attention_scores):
        attn_score = []
        for i in range(1, attention_score.shape[-1]-kmer+2):
            attn_score.append(float(attention_score[:,0,i].sum()))
        for i in range(len(attn_score)-1):
            if attn_score[i+1] == 0:
                attn_score[i] = 0
                break
        counts = np.zeros([len(attn_score)+kmer-1])
        real_scores = np.zeros([len(attn_score)+kmer-1])
        for i, score in enumerate(attn_score):
            for j in range(kmer):
                counts[i+j] += 1.0
                real_scores[i + j] += score
        real_scores = real_scores / counts
        real_scores = real_scores / np.linalg.norm(real_scores)
        scores[index] = real_scores
    return scores


X = pd.DataFrame({'seq': [str(s) for s in pos_seqs + neg_seqs]})
Y = np.concatenate([np.tile([0, 1], (len(pos_seqs), 1)), np.tile([1, 0], (len(neg_seqs), 1))], axis=0)

sequences, labels, scores = model_eval(X, Y, clf_model.bert, tokenizer, window_size, kmer)

dev = pd.DataFrame()
dev['seq'] = sequences
dev['label'] = labels
dev_pos = dev[dev['label'] == 1]
dev_neg = dev[dev['label'] == 0]
pos_atten_scores = scores[dev_pos.index.values]
neg_atten_scores = scores[dev_neg.index.values]

if not os.path.exists('./dump_files/'):
        os.makedirs('./dump_files/')

motif_seqs = mu.find_high_att_regions(pos_atten_scores, list(dev_pos['seq']), min_len=6, kwargs={})
data_reader.write_fasta(motif_seqs.keys(), './dump_files/high_atten_seqs_pos.fasta')

motif_seqs = mu.find_high_att_regions(neg_atten_scores, list(dev_neg['seq']), min_len=6, kwargs={})
data_reader.write_fasta(motif_seqs.keys(), './dump_files/high_atten_seqs_neg.fasta')

print('high attention regions saved in dump_files folder'.format(len(motif_seqs)))

print('STREME is running .....')
subprocess.run(["streme", "--p", "./dump_files/high_atten_seqs_pos.fasta", "--n", "./dump_files/high_atten_seqs_neg.fasta", "--verbosity", "1", "--minw", "6", "--maxw", "12", "--dna", "--nmotifs", "100", "--oc", "./dump_files/"])
print('STREME DONE!')

all_motifs = mu.parse_motif_file("./dump_files/streme.txt")

selected_motifs = sorted(all_motifs.items(), key=lambda item: float(item[1][0]))[:3]

mu.write_selected_motifs("./results/streme_selected.txt", selected_motifs)


print('SELECTED MOTIFS WROTE to file ./results/streme_selected.txt')

#python motif_finding.py -dmc_seq ./data/motif_sample_pos.fa -ndmc_seq ./data/motif_sample_neg.fa

