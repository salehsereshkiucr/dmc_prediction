import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)
import torch
import torch.nn as nn
from transformers import BertModel
import glob

class DNASequenceClassifier(nn.Module):
    def __init__(self, bert_model, fc_ad, num_classes, **kwargs):
        super(DNASequenceClassifier, self).__init__()
        self.bert = bert_model
        fc_first_layer_size = bert_model.config.hidden_size + kwargs['meth_window_size'] if 'meth_window_size' in kwargs.keys() else bert_model.config.hidden_size
        self.fc = nn.Sequential(
            nn.Linear(fc_first_layer_size, 128),
            nn.Dropout(0.5),
            nn.ReLU(),
            nn.Linear(128, 24),
            nn.Dropout(0.5),
            nn.ReLU(),
            nn.Linear(24, num_classes),
            nn.Softmax()
        ).to(torch.float32)
        if fc_ad != None:
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
            self.fc.load_state_dict(torch.load(fc_ad, map_location=torch.device(device)))
        self.include_meth = 'meth_window_size' in kwargs.keys()
    def forward(self, input_ids, attention_mask, methylations=None):
        #output_attentions=True
        outputs = self.bert(input_ids=input_ids, attention_mask=attention_mask)
        pooled_output = outputs[1]  # Get pooled output pooled_output.shape = torch.Size([batch_size, 768])
        # if isinstance(pooled_output, tuple):
        #     pooled_output = outputs[0][:, 0, :]
        #pooled_output = outputs.last_hidden_state[:, 0, :]
        if self.include_meth: pooled_output = torch.cat((pooled_output, methylations), dim=1).to(torch.float32)
        logits = self.fc(pooled_output)
        return logits, outputs
    def save(self, file_name):
        self.bert.save_pretrained(file_name+'_bert')
        torch.save(self.fc.state_dict(), file_name+'_torchnn.pth')


def find_torch_and_bert_paths(folder_path):
    torch_pattern = os.path.join(folder_path, '*_torchnn.pth')
    bert_pattern = os.path.join(folder_path, '*_bert')
    torch_files = glob.glob(torch_pattern)
    bert_folders = [d for d in glob.glob(bert_pattern) if os.path.isdir(d)]
    if len(torch_files) != 1:
        raise ValueError(f"Expected exactly one torch model file, but found {len(torch_files)}.")
    if len(bert_folders) != 1:
        raise ValueError(f"Expected exactly one bert folder, but found {len(bert_folders)}.")
    return torch_files[0], bert_folders[0]


def load_clf_model(file_name):
    trc_add, brt_add = find_torch_and_bert_paths(file_name)
    return DNASequenceClassifier(BertModel.from_pretrained(brt_add), trc_add, 2)
