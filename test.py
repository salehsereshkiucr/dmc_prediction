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
import utils.data_manager as dmngr
parser = argparse.ArgumentParser()
parser.add_argument('-mc', '--mc_address', help='address of a tab separated file with the positions of cytosines', required=True)
parser.add_argument('-ga', '--genome_assembly', help='address of a tab separated file with the positions of not differentially methylated cytosines', required=True)
parser.add_argument('-mn', '--model_name', help='the address to load model', required=False, default='./models/trained_clf')
parser.add_argument('-ws', '--window_size', help='window_size', required=False, default=512)
parser.add_argument('-kmer', '--kmer', help='kmer', required=False, default=6)

args = parser.parse_args()

mc_address = args.mc_address
genome_assembly = args.genome_assembly
model_address = args.model_name
kmer = int(args.kmer)
window_size = int(args.window_size)


# kmer = 6
# window_size = 512
# mc_address = './data/sample_dmc_test.csv'
# genome_assembly = './data/sample_asembly.fa'
# model_address = './models/trained_clf'


clf_model = mdlmngr.load_clf_model(model_address)
tokenizer = BertTokenizer.from_pretrained("zhihan1996/DNA_bert_"+str(kmer), trust_remote_code=True)
print('model is loaded....', model_address)
assembly = data_reader.readfasta(genome_assembly)
mcs = pd.read_csv(mc_address, sep='\t', header=None, names=['chr', 'position'])
print('data is loaded... ')

dataset = dmngr.make_input_dataset(assembly, mcs, -1, window_size=window_size)
dataset = dataset.sample(frac=1).reset_index(drop=True)

X = dataset.drop(columns=['label'])
Y = np.asarray(pd.cut(dataset['label'], bins=2, labels=[0, 1], right=False))
b = np.zeros((Y.size, Y.max()+1))
b[np.arange(Y.size), Y] = 1
Y = b

batch_size = 32
labels = torch.tensor(Y, dtype=torch.float32)
dataset = DNADataset(X, labels, tokenizer, window_size, kmer)
data_loader_test = DataLoader(dataset, batch_size=batch_size, shuffle=True)

total_samples_test = len(data_loader_test.dataset)

device = 'cuda' if torch.cuda.is_available() else 'cpu'

clf_model.to(device)
clf_model.eval()
test_correct = 0
#keep track of the predicted and labels to draw the ROC curve
predicted_test = np.zeros((total_samples_test, 2)) - 1
idx = 0

with torch.no_grad():
    for batch in data_loader_test:
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        target_labels = batch["label"].to(device)
        input_ids.requires_grad = False
        attention_mask.requires_grad = False
        target_labels.requires_grad = False
        logits, _ = clf_model(input_ids, attention_mask)
        predicted_test[idx: idx + len(logits)] = logits.cpu().numpy()
        idx += len(logits)

assert -1 not in predicted_test
if not os.path.exists('./results/'):
        os.makedirs('./results/')

with open('./results/test_res.txt', 'w') as file:
    for number in np.argmax(predicted_test, axis=1):
        file.write(f"{number}\n")

print('results saved in {}'.format('./results/test_res.txt'))

#python test.py -mc ./data/sample_dmc_test.csv -ga ./data/sample_asembly.fa
