import os
import subprocess
from operator import itemgetter
from collections import Counter

def mafft(align_seq_file):
    ''' use MAFFT to create MSA '''
    # 通过记录的reads片段进行比对得到比对结果

    out_file_msa = '.'.join(align_seq_file.split('.')[:-1]) + '.msa.fa'
    
    # 根据序列的长度进行参数的选择
    #     if 
    args = ['mafft', '--randomseed', '1', align_seq_file]
 
    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = FNULL)

    with open(out_file_msa, 'w') as out_fa:
        for line in p.stdout:
            line = line.decode()
            # if line[0] != '>':
            out_fa.write(line)
    return out_file_msa
   




def get_consensus_seq(sequences, fre_cuoff):
    # 遍历比对结果，根据cutoff取出每个位置的consensus结果
    ids = '*'
    if list(sequences.keys()) != [None]:
        ids = ';'.join(list(sequences.keys()))
    seq_len = len(list(sequences.values())[0])
    consensus_seq = ''
    consensus_seq2 = ''
    for i in range(seq_len):
        try:
            seq_i = [ seq[i] for seq in list(sequences.values()) ]
        except IndexError:
            print(ids)

        set_seq_i = list(set(seq_i))
        
        # 第一类结果，不算gap，每个取出现概率较大的那个，并且概率要在80%以上
        # 必须完全consensus
        # 给consensus设置一个cutoff

        # consensus_seq
        seq_i_A = seq_i.count('A')
        seq_i_T = seq_i.count('T')
        seq_i_C = seq_i.count('C')
        seq_i_G = seq_i.count('G')
        seq_i_GAP = seq_i.count('-')
        seq_i_nucle = sorted([('A', seq_i_A), ('T', seq_i_T), ('C', seq_i_C), ('G', seq_i_G)], key = itemgetter(1))[-1]
        
        if seq_i_nucle[1] < len(seq_i) * fre_cuoff or seq_i_GAP >= len(seq_i) * fre_cuoff :
            seq = '-'
        else:
            seq = seq_i_nucle[0]
        consensus_seq = consensus_seq + seq

        # 第二类结果，把“—”gap也算进去，取出现频次最高的那个
        seq_i_count = sorted( Counter(seq_i).items(), key=itemgetter(1) )
        if seq_i_count[-1][1] >= 2:
            consensus_seq2 = consensus_seq2 + seq_i_count[-1][0] 
        else:
            consensus_seq2 = consensus_seq2 + '-'

    return consensus_seq


def consensus_tsd(align_seq_file, fre_cuoff):
    # 把mafft比对结果处理得到包含比对结果的列表
    # 然后从列表得到consensus结果
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
    align_seq = get_consensus_seq(align_sequences, fre_cuoff)
    return align_seq


