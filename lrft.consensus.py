import os
import sys
import subprocess
import pickle
from operator import itemgetter
from collections import Counter

PROJECT_PATH = "/data/tusers/boxu/lrft/result/nanopore/human"
BED_PATH = PROJECT_PATH + "/insertion/tmp"
PREFIX="SS"
SEQ_INFO_FILE = PROJECT_PATH + PREFIX.split('.')[0] + ".clip.reads.pkl"

# 多序列比对
def mafft(align_seq_file):
    ''' use MAFFT to create MSA '''

    out_fn = '.'.join(align_seq_file.split('.')[:-1]) + '.msa.fa'
	
	# 根据序列的长度进行参数的选择
	# 	if 
    args = ['mafft', '--randomseed', '1', align_seq_file]
 
    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = FNULL)

    with open(out_fn, 'w') as out_fa:
        for line in p.stdout:
            line = line.decode()
            # if line[0] != '>':
            out_fa.write(line)
    return out_fn
   
def get_consensus_seq(sequences):
	ids = '*'
	# print(sequences)
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

		seq_i_A = seq_i.count('A')
		seq_i_T = seq_i.count('T')
		seq_i_C = seq_i.count('C')
		seq_i_G = seq_i.count('G')
		seq_i_GAP = seq_i.count('-')
		seq_i_nucle = sorted([('A', seq_i_A), ('T', seq_i_T), ('C', seq_i_C), ('G', seq_i_G)], key = itemgetter(1))[-1]
		# 第一类结果，不算gap，每个取出现概率较大的那个，并且概率要在80%以上
		if seq_i_nucle[1] == 0 or seq_i_GAP >= len(seq_i) * 0.8:
			seq = ''
		else:
			seq = seq_i_nucle[0]
		consensus_seq = consensus_seq + seq

		# 一类结果，把“—”gap也算进去，取出现频次最高的那个
		seq_i_count = sorted( Counter(seq_i).items(), key=itemgetter(1) )
		consensus_seq2 = consensus_seq2 + seq_i_count[-1][0] 

	# print(ids)
	return [ids, consensus_seq, consensus_seq2]

with open(SEQ_INFO_FILE, "rb") as f3:
	READS_SEQUENCE = pickle.load(f3)



for BED_FILE in os.walk(BED_PATH):
    print(BED_FILE)
    for line in open(BED_FILE, 'r'):
        line = line.strip().split('\t')

