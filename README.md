
# L-Map: Detailed Module Documentation

## Overview

L-Map is a tool designed for the prediction and analysis of differentially methylated cytosines (DMCs) using DNABERT, a Large Language Model. It has several modules: Prediction, motif finding, Training, and DMC detection through dedicated modules.

## Installation

This project requires Python 3.6 or later. Before installing the required dependencies, ensure you have Python and pip installed on your system.

### Setting Up a Virtual Environment (Recommended)

It's recommended to create a virtual environment for this project to avoid conflicts with other packages. Use the following commands to create and activate a virtual environment:

# Create a virtual environment (replace 'myenv' with your preferred environment name)
```
python -m venv myenv
```

# Activate the virtual environment
# On Windows
```
myenv\Scripts\activate
```
# On macOS and Linux
```
source myenv/bin/activate
```

# Installing Dependencies
With the virtual environment activated, install the project dependencies by running:
```
pip install -r requirements.txt
```

## Modules

### Prediction

The prediction module enables the use of a trained model to classify a list of cytosines as DMCs (Differentially Methylated Cytosines) or not. For each of the seven knock-out datasets described in the paper "Sequence level TETs and DNMT3s domain prediction via a Language Model," we provide one pre-trained model. These models can be accessed and downloaded from [this Google Drive Link](https://drive.google.com/drive/folders/1jCDAdqcvpeI9nNrWiJAYQh1jgauOaO2E?usp=sharing)

To use these pre-trained models, download them and specify the downloaded folder's path in the module command's `model_name` attribute. This attribute should point to a directory containing the []_torchnn.pth file and a []_bert folder, which includes both a pytorch_model.bin and a config.json file.

Additionally, you must provide the genome assembly and a file listing the cytosine positions in the attributes. Sample datasets can be found in ./data/.

**Important:** The files located in ./models/pretrained_models/TET/ are placeholders due to GitHub's file size limitations and cannot be used directly to load a model. Before using this module, ensure you have downloaded a model from [here]((https://drive.google.com/drive/folders/1jCDAdqcvpeI9nNrWiJAYQh1jgauOaO2E?usp=sharing)) or trained your own model.

**Command:**
```
python test.py -mc <mc_address> -ga <genome_assembly> -mn <model_address>
```

**Options:**
- `-mc`, `--mc_address`: File with cytosines to classify. (Required)
- `-ga`, `--genome_assembly`: Genome assembly file in FASTA. (Required)
- `-mn`, `--model_name`: Location of the trained model. (Required)
- `-ws`, `--window_size`: Window size for sequence processing. (Default: 512)
- `-kmer`: K-mer size for sequence representation. (Default: 6)


**Example Command:**

```
python test.py -mc ./data/sample_dmc_test.csv -ga ./data/sample_asembly.fa -mn ./models/pretrained_models/TET/ -ws 512 -kmer 6
```

### Motif Finding 

This module is aimed at identifying motifs corresponding to DMCs using a trained model. Sample datasets are available at ./data/

**Command:**
```
python motif_finding.py -dmc_seq <dmc_seq_address> -ndmc_seq <ndmc_seq_address> [options]
```

**Options:**
- `-dmc_seq`, `--dmc_seq`: Context sequences of DMCs in FASTA. (Required)
- `-ndmc_seq`, `--ndmc_seq`: Context sequences of nDMCs in FASTA. (Required)
- `-mn`, `--model_name`: Location of the trained model. (Default: './models/trained_clf')
- `-ws`, `--window_size`: Window size for sequence processing. (Default: 512)
- `-kmer`: K-mer size for sequence representation. (Default: 6)

**Example Command:**

```
python motif_finding.py -dmc_seq ./data/motif_sample_pos.fa -ndmc_seq ./data/motif_sample_neg.fa -mn ./models/trained_clf -ws 512 -kmer 6
```




### Training

To use the training module, you'll need to prepare two tab-separated files: one for the positions of differentially methylated cytosines (DMCs) and another for the positions of not differentially methylated cytosines (nDMCs), as well as the genome assembly in FASTA format. Sample datasets are available at ./data/

Run the training module with the following command:

**Command:**
```
python train.py -dmc <dmc_address> -ndmc <not_dmc_address> -ga <genome_assembly> [options]
```

**Options:**
- `-dmc`, `--dmc_address`: Address of the file with DMCs. (Required)
- `-ndmc`, `--not_dmc_address`: Address of the file with nDMCs. (Required)
- `-ga`, `--genome_assembly`: Genome assembly file in FASTA format. (Required)
- `-kmer`: K-mer size for sequence representation. (Default: 6)
- `-ws`, `--window_size`: Window size for sequence processing. (Default: 512)
- `-bs`, `--batch_size`: Batch size for training. (Default: 32)
- `-mn`, `--model_name`: Save location for the trained model. (Default: './models/trained_clf')


**Example Command:**
```
python train.py -dmc ./data/sample_dmc_train.csv -ndmc ./data/sample_notdmc_train.csv -ga ./data/sample_asembly.fa -kmer 6 -ws 512 -bs 32 -mn ./models/trained_clf
```

This command will train a model using the specified DMC and nDMC files, with the provided genome assembly, and save it to the specified location.



### DMC Detection

This module is designed to analyze methylation levels and identify DMCs and non-DMCs based on the provided criteria. Sample datasets are available at ./data/

**Command:**
```
python dmc_detection.py -c <control_address> -t <test_address> [options]
```

**Options:**
- `-c`, `--control`: BED file with methylation levels in control. (Required)
- `-t`, `--test`: BED file with methylation levels in test. (Required)
- `-ct`, `--coverage_threshold`: Minimum coverage for classification. (Default: 10)
- `-md`, `--minimum_difference`: Minimum methylation level difference for DMC classification. (Default: 0.6)

**Example Command:**

```
python dmc_detection.py -t ./data/test_mC_report.bed -c ./data/control_mC_report.bed -ct 10 -md 0.6
```
