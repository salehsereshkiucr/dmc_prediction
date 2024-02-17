
# L-Map: Detailed Module Documentation

## Overview

L-Map is a tool designed for the prediction and analysis of differentially methylated cytosines (DMCs) using DNABERT, a Large Language Model. It has several modules: training, testing, motif finding, and DMC detection through dedicated modules.

## Prerequisites

Ensure you have the following libraries installed:
- pandas (2.0.3)
- numpy (1.24.2)
- sklearn (1.3.0)
- transformers (4.18.0)
- tensorflow (2.11.0)
- pytorch (2.0.1+cu117)

## Modules

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


### Testing

The testing module allows for the utilization of a trained model to predict a list of cytosines as DMCs or not. Sample datasets are available at ./data/

**Command:**
```
python test.py -mc <mc_address> -ga <genome_assembly> [options]
```

**Options:**
- `-mc`, `--mc_address`: File with cytosines to classify. (Required)
- `-ga`, `--genome_assembly`: Genome assembly file in FASTA. (Required)
- `-mn`, `--model_name`: Location of the trained model. (Default: './models/trained_clf')
- `-ws`, `--window_size`: Window size for sequence processing. (Default: 512)
- `-kmer`: K-mer size for sequence representation. (Default: 6)


**Example Command:**

```
python test.py -mc ./data/sample_dmc_test.csv -ga ./data/sample_asembly.fa -mn ./models/trained_clf -ws 512 -kmer 6
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
