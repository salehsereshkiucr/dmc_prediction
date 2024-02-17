
import numpy as np
import re


def contiguous_regions(condition, len_thres=5):
    """
    Modified from and credit to: https://stackoverflow.com/a/4495197/3751373
    Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    Arguments:
    condition -- custom conditions to filter/select high attention
            (list of boolean arrays)
    Keyword arguments:
    len_thres -- int, specified minimum length threshold for contiguous region
        (default 5)
    Returns:
    idx -- Index of contiguous regions in sequence
    """
    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()
    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit
    # Reshape the result into two columns
    idx.shape = (-1,2)
    # eliminate those not satisfying length of threshold
    idx = idx[np.argwhere((idx[:,1]-idx[:,0])>=len_thres).flatten()]
    return idx


#This method is from DNABERT motif finding module
def find_high_attention(score, min_len=5, **kwargs):
    """
    With an array of attention scores as input, finds contiguous high attention
    sub-regions indices having length greater than min_len.
    Arguments:
    score -- numpy array of attention scores for a sequence
    Keyword arguments:
    min_len -- int, specified minimum length threshold for contiguous region
        (default 5)
    **kwargs -- other input arguments:
        cond -- custom conditions to filter/select high attention
            (list of boolean arrays)
    Returns:
    motif_regions -- indices of high attention regions in sequence
    """
    cond1 = (score > np.mean(score))
    cond2 = (score > 10*np.min(score))
    cond = [cond1, cond2]
    cond = list(map(all, zip(*cond)))
    if 'cond' in kwargs: # if input custom conditions, use them
        cond = kwargs['cond']
        if any(isinstance(x, list) for x in cond): # if input contains multiple conditions
            cond = list(map(all, zip(*cond)))
    cond = np.asarray(cond)
    # find important contiguous region with high attention
    motif_regions = contiguous_regions(cond, min_len)
    return motif_regions



#This method is from DNABERT motif finding module
def find_high_att_regions(atten_scores, seqs, min_len, kwargs):
    motif_seqs = {}
    ## find the motif regions
    # It findind the high attention region in each of the seqeunces then fill the motif_seqs dictionary.
    # The dictionary contains the sequences of motifs as keys and then two lists as values. one to keep sequence index and the other to keep the coordinates of motif in the sequence.
    print("* Finding high attention motif regions")
    for i, score in enumerate(atten_scores):
        seq_len = len(seqs[i])
        score = score[0:seq_len]
        # handle kwargs
        if 'atten_cond' in kwargs:
            motif_regions = find_high_attention(score, min_len=min_len, cond=kwargs['atten_cond'])
        else:
            motif_regions = find_high_attention(score, min_len=min_len)
        for motif_idx in motif_regions:
            seq = seqs[i][motif_idx[0]:motif_idx[1]]
            if seq not in motif_seqs:
                motif_seqs[seq] = {'seq_idx': [i], 'atten_region_pos':[(motif_idx[0],motif_idx[1])]}
            else:
                motif_seqs[seq]['seq_idx'].append(i)
                motif_seqs[seq]['atten_region_pos'].append((motif_idx[0],motif_idx[1]))
    return motif_seqs


def parse_motif_file(filename):
    motifs = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
    motif_name = None
    e_value = None
    matrix = []
    for line in lines:
        # Check for motif line and extract motif name
        if line.startswith('MOTIF') or line.startswith('***'):
            if motif_name and e_value and matrix:
                motifs[motif_name] = (e_value, matrix)
            if line.startswith('MOTIF'): motif_name = line.split()[1]
            matrix = []
            e_value = None
        elif 'E= ' in line and motif_name and e_value is None:
            e_value_match = re.search(r'E=\s*([0-9.e+-]+)', line)
            if e_value_match:
                e_value = e_value_match.group(1)
        elif motif_name and e_value and line.strip() and not line.startswith('MOTIF'):
            matrix_row = [float(x) for x in line.strip().split()]
            matrix.append(matrix_row)
    if motif_name and e_value and matrix:
        motifs[motif_name] = (e_value, matrix)
    return motifs

def write_selected_motifs(out_file, s_mtfs):
    output = ''
    for mtf, (pval, mtrx) in s_mtfs:
        output += mtf + ' e_val = ' + str(pval) + '\n'
        for row in mtrx:
            output += ' '.join(str(item) for item in row) + '\n'
        output += '\n'
    with open(out_file, 'w') as outputfile:
        outputfile.write(output)



