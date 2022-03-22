import os
import subprocess
from operator import itemgetter
from collections import Counter



# lrft_consensus_seq 和 left_consensus_tsd 两个方法几乎相同，想想能不能合并到一个，这样就可以省些代码，不然显得累赘

def mafft(align_seq_file):
    ''' use MAFFT to create MSA '''

    out_file_msa = '.'.join(align_seq_file.split('.')[:-1]) + '.msa.fa'
    
    # 根据序列的长度进行参数的选择
    #     if
    args = [ 'mafft', '--randomseed', '1', align_seq_file ]
 
    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = FNULL)

    with open(out_file_msa, 'w') as out_fa:
        for line in p.stdout:
            line = line.decode()
            # if line[0] != '>':
            out_fa.write(line)
    return out_file_msa
   
def consensus_seq(sequences, fre_cuoff, type):
    ids = '*'
    if list(sequences.keys()) != [None]:
        ids = ';'.join(list(sequences.keys()))
    seq_len = len(list(sequences.values())[0])
    consensus_seq1 = ''
    consensus_seq2 = ''
    for i in range(seq_len):
        try:
            seq_i = [ seq[i] for seq in list(sequences.values()) ]
        except IndexError:
            print(ids)
        # 第一类结果，不算gap，每个取出现概率较大的那个，并且概率要在80%以上
        # consensus_seq
        seq_i_A = seq_i.count('A')
        seq_i_T = seq_i.count('T')
        seq_i_C = seq_i.count('C')
        seq_i_G = seq_i.count('G')
        seq_i_GAP = seq_i.count('-')
        seq_i_nucle = sorted([('A', seq_i_A), ('T', seq_i_T), ('C', seq_i_C), ('G', seq_i_G)], key = itemgetter(1))[-1]
        
        seq = seq_i_nucle[0]
        if type == "seq":
            if seq_i_nucle[1] == 0 or seq_i_GAP >= len(seq_i) * fre_cuoff:
                seq = ''
        elif type == 'tsd':
            if seq_i_nucle[1] < len(seq_i) * fre_cuoff or seq_i_GAP >= len(seq_i) * fre_cuoff:
                seq = '-'
        
        consensus_seq1 = consensus_seq1 + seq

        # 第二类结果，把“—”gap也算进去，取出现频次最高的那个
        seq_i_count = sorted( Counter(seq_i).items(), key=itemgetter(1) )
        if seq_i_count[-1][1] >= 2 :
            consensus_seq2 = consensus_seq2 + seq_i_count[-1][0] 
        else :
            consensus_seq2 = consensus_seq2 + '-'

        consensus_seq2 = consensus_seq2 + seq_i_count[-1][0] 

    # print(ids)
    return [ids, consensus_seq1, consensus_seq2]

def get_consensus_seq(align_seq_file, fre_cuoff, type):
    align_sequences = {}
    id = None
    seq = ''
    for line in open(mafft( align_seq_file ), 'r'):
        if line.startswith('>'):
            if id != None:
                align_sequences[id] = seq.upper()
            id = line.strip().split('>')[1]
            seq = ''
        else:
            seq = seq + line.strip()
    align_sequences[id] = seq.upper()
    
    if type == "seq":
        # align_seq = '\t|\t'.join(consensus_seq(align_sequences, fre_cuoff, type))
        align_seq = consensus_seq(align_sequences, fre_cuoff, type)[1]
    elif type == "tsd":
        # print(type)
        align_seq = consensus_seq(align_sequences, fre_cuoff,type)[1]
    return align_seq
