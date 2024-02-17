import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)
import torch
import torch.nn as nn
from transformers import BertModel

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
            self.fc.load_state_dict(torch.load(fc_ad))
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

def load_clf_model(file_name):
    return DNASequenceClassifier(BertModel.from_pretrained(file_name+'_bert'), file_name+'_torchnn.pth', 2)
