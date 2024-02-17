import sys
import os

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

import utils.data_reader as data_reader
import torch.optim as optim
from tqdm import tqdm
import pandas as pd
import torch
from torch.utils.data import DataLoader
import torch.nn as nn
from sklearn.model_selection import train_test_split
import numpy as np
import argparse
from transformers import BertTokenizer, BertModel
from utils.model_manager import DNASequenceClassifier
from utils.data_manager import DNADataset
from utils import data_manager as dmngr
from transformers import AutoModel, AutoTokenizer, AutoModelForMaskedLM
from utils import preprocess


parser = argparse.ArgumentParser()
parser.add_argument('-dmc', '--dmc_address', help='address of a tab separated file with the positions of differentially methylated cytosines', required=True)
parser.add_argument('-ndmc', '--not_dmc_address', help='address of a tab separated file with the positions of not differentially methylated cytosines', required=True)
parser.add_argument('-ga', '--genome_assembly', help='address of a tab separated file with the positions of not differentially methylated cytosines', required=True)
parser.add_argument('-kmer', '--kmer', help='kmer', required=False, default=6)
parser.add_argument('-ws', '--window_size', help='window_size', required=False, default=512)
parser.add_argument('-bs', '--batch_size', help='batch size', required=False, default=32)
parser.add_argument('-mn', '--model_name', help='the address to save model', required=False, default='./models/trained_clf')

args = parser.parse_args()

kmer = int(args.kmer)
window_size = int(args.window_size)
batch_size = int(args.batch_size)
dmc_address = args.dmc_address
not_dmc_address = args.not_dmc_address
genome_assembly = args.genome_assembly
model_address = args.model_name
#include_methylation = args.include_methylation

# kmer = 6
# window_size = 512
# batch_size = 32
# dmc_address = './data/sample_dmc_train.csv'
# not_dmc_address = './data/sample_notdmc_train.csv'
# genome_assembly = './data/sample_asembly.fa'
# model_address = './models/trained_clf'


num_classes = 2  # Binary classification

model_name = "zhihan1996/DNA_bert_"+str(kmer)
model_name_alias = 'brt'
model = BertModel.from_pretrained(model_name, num_labels=2, finetuning_task="dnaprom", cache_dir=None)
tokenizer = BertTokenizer.from_pretrained(model_name, trust_remote_code=True)
clf_model = DNASequenceClassifier(model, None, num_classes) #If meth_window_size specified here it automatically considers including methylation profiles

assembly = data_reader.readfasta(genome_assembly)
dmrs, not_dmrs = dmngr.load_df(dmc_address, not_dmc_address)

dataset = pd.concat([dmngr.make_input_dataset(assembly, dmrs, 1, window_size=window_size, meth_seq=None),
                     dmngr.make_input_dataset(assembly, not_dmrs, 0, window_size=window_size, meth_seq=None)])
dataset = dataset.sample(frac=1).reset_index(drop=True)
dataset = dataset[dataset['seq'].str.len() == window_size]

X = dataset.drop(columns=['label'])
Y = np.asarray(pd.cut(dataset['label'], bins=2, labels=[0, 1], right=False))
b = np.zeros((Y.size, Y.max()+1))
b[np.arange(Y.size), Y] = 1
Y = b
# X.shape (N, 1)|| X.columns ['seq'] || type(X) <class 'pandas.core.frame.DataFrame'> || Y.shape (N, 2) || type(Y) <class 'numpy.ndarray'>


labels = torch.tensor(Y, dtype=torch.float32)
dataset = DNADataset(X, labels, tokenizer, window_size, kmer)
data_loader_train = DataLoader(dataset, batch_size=batch_size, shuffle=True)

device = 'cuda' if torch.cuda.is_available() else 'cpu'

optimizer = optim.Adam(clf_model.parameters(), lr=1e-5)
criterion = nn.BCELoss()
clf_model.to(device)
# Set the number of epochs
num_epochs = 5
total_samples = len(data_loader_train.dataset)  # Total number of samples
batch_size = data_loader_train.batch_size  # Batch size
# Training loop
for epoch in range(num_epochs):
    clf_model.train()
    total_loss = 0
    total_correct = 0
    total_processed = 0  # To keep track of total processed samples
    progress_bar = tqdm(data_loader_train, desc=f"Epoch {epoch+1}/{num_epochs}", leave=False)
    for batch in progress_bar:
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        target_labels = batch["label"].to(device).to(torch.float32)
        methylations = batch["methylations"].to(device)
        optimizer.zero_grad()
        logits, _ = clf_model(input_ids, attention_mask, methylations=methylations)
        loss = criterion(logits, target_labels)
        loss.backward()
        optimizer.step()
        predictions = logits.argmax(dim=1)
        total_correct += (predictions == target_labels.argmax(dim=1)).sum().item()
        total_processed += batch_size
        if total_processed % (10*batch_size) == 0:
            b_acc = total_correct / total_processed
            progress_percentage = (total_processed / total_samples) * 100
            progress_bar.set_postfix(loss=loss.item(), batch_acc=b_acc, progress=f"{progress_percentage:.2f}%")
        total_loss += loss.item()
    epoch_accuracy = total_correct / total_samples
    average_loss = total_loss / len(data_loader_train)
    print(f"Epoch [{epoch+1}/{num_epochs}] - Average Loss: {average_loss:.4f}, Epoch Accuracy: {epoch_accuracy:.4f}")
    # Evaluate the model on test data

clf_model.save(model_address)
print('model saved in in {}'.format(model_address))


#python train.py -dmc ./data/sample_dmc_train.csv -ndmc ./data/sample_notdmc_train.csv -ga ./data/sample_asembly.fa
