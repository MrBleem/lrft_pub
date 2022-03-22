import datetime
starttime = datetime.datetime.now()

import re
import sys
import pysam
import pickle
# import numpy as np
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt

# python ${BIN_PATH}/lrft_clip_reads.py ${BAM} ${PREFIX}.mapped.sorted.TE.clip.fq ${PREFIX}.mapped.sorted.TE.clip.bed ${PREFIX}.clip.reads.pkl

bamFile=sys.argv[1]
outFile_fq=sys.argv[2]
outFile_bed=sys.argv[3]
outFile_pickle=sys.argv[4]

# bamFile = "/data/tusers/boxu/lrft/result/nanopore/human/map/SRR11669560.mapped.sorted.TE.bam"
# outFile_fq = "/data/tusers/boxu/lrft/result/nanopore/human/map/t1.fq"
# outFile_bed = "/data/tusers/boxu/lrft/result/nanopore/human/map/t1.bed"
# outFile_pickle = "/data/tusers/boxu/lrft/result/nanopore/human/map/t1.pickle"


bamfile = pysam.AlignmentFile(bamFile, "rb")
g1=open(outFile_fq,'w')
g2=open(outFile_bed,'w')
g2.write('name\tref\tmap_pos\tl_clip_seq_l\tinsertion_len\tl_clip_seq_r\tinsertion_seq\tmap_qua\n')


class reads_info:
	def __init__(self,read):
		self.name = read.qname
		self.cigar = read.cigarstring
		self.mapq = read.mapq
		self.seq = read.query_sequence
		self.ref = read.reference_name
		self.map_pos = read.pos+1
		self.flag = str(read.flag)
		flag = hex(read.flag)
		flag_strand = int(flag[2:])%100
		if flag_strand != 10 :
			self.strand = '+'
		else:
			self.strand = '-'

		self.map_qua = ''
		if read.query_qualities != None:
			for m in read.query_qualities:
				m_tmp = int(m) + 33
				self.map_qua = self.map_qua + chr(m_tmp)


# SRR11669560.sra.12196843
# SRR11669560.sra.12196667


def comp_reverse(seq):
	base = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	new_seq = ''
	for i in range(len(seq)):
		s = seq[i]
		new_seq = base[s] + new_seq
	return new_seq


def record(read_info, g1, g2, READS_PICKLE, t='skip'):
	s_cigar = re.findall('\d+[HS]', read_info.cigar)
	cigar_end = read_info.cigar[-1]

	l_clip_seq_r = 0
	l_clip_seq_l = 0
	if len(s_cigar) == 0:
		return READS_PICKLE
	if len(s_cigar) == 1:
		l_clip_seq = int(re.split('[HS]', s_cigar[0])[0])
		if cigar_end == 'S' or cigar_end == 'H' :
			clip_seq = read_info.seq[ -l_clip_seq: ]
			clip_mq = read_info.map_qua[ -l_clip_seq: ]
			l_clip_seq_r = l_clip_seq
			insertion_seq = read_info.seq[ :-l_clip_seq_r ]
			insertion_seq_mq = read_info.map_qua[ :-l_clip_seq_r ]
		else:
			clip_seq = read_info.seq[ :l_clip_seq ]
			clip_mq = read_info.map_qua[ :l_clip_seq ]
			l_clip_seq_l = l_clip_seq
			insertion_seq = read_info.seq[ l_clip_seq_l: ]
			insertion_seq_mq = read_info.map_qua[ l_clip_seq_l: ]
	if len(s_cigar) == 2:
		l_clip_seq_l = int(re.split('[HS]', s_cigar[0])[0])
		l_clip_seq_r = int(re.split('[HS]', s_cigar[1])[0])
		clip_seq = read_info.seq[ :l_clip_seq_l ] + read_info.seq[ -l_clip_seq_r: ]
		clip_mq = read_info.map_qua[ :l_clip_seq_l ] + read_info.map_qua[ -l_clip_seq_r: ]
		insertion_seq = read_info.seq[ l_clip_seq_l: -l_clip_seq_r ]
		insertion_seq_mq = read_info.map_qua[ l_clip_seq_l: -l_clip_seq_r ]
	

	insertion_len = len(read_info.seq) - l_clip_seq_l - l_clip_seq_r
	new_reads_name = ','.join( [ '@'+read_info.name, read_info.ref, read_info.strand, str(read_info.map_pos), str(l_clip_seq_l), str(insertion_len), str(l_clip_seq_r) ]) # reads name的分隔符
	if insertion_len >= 200:
		g1.write(new_reads_name + ' ' + 'length=' + str(len(clip_seq)) + '\n')
		g1.write(clip_seq+'\n')
		g1.write('+\n')
		g1.write(clip_mq+'\n')
		
		new_reads=read_info.name+'\t'+read_info.ref+'\t'+str(read_info.map_pos)+'\t'+str(l_clip_seq_l)+'\t'+str(insertion_len)+'\t'+str(l_clip_seq_r)
		g2.write(new_reads+'\t'+insertion_seq+'\t'+insertion_seq_mq+'\n')

		reads_id = ','.join([read_info.name, read_info.ref, read_info.strand, str(read_info.map_pos), str(l_clip_seq_l), str(insertion_len), str(l_clip_seq_r)] ) # read_info.name+'.'+read_info.ref+'.'+read_info.strand+'.'+str(read_info.map_pos)+'.'+str(l_clip_seq_l)+'.'+str(insertion_len)+'.'+str(l_clip_seq_r)
		reads_seq = insertion_seq
		READS_PICKLE[reads_id] = reads_seq
	return READS_PICKLE



READS_SEQ_INFO = {}
READS_PICKLE = {}
bamfile = pysam.AlignmentFile(bamFile, "rb")
# 开始拆分reads
for read in bamfile:
	read_info = reads_info(read)
	# if read_info.seq != None:
	# 	print('seq:'+str(len(read_info.seq)))
	# else:
	# 	print('seq:0')

	# if read_info.seq != None and read_info.seq != '*' and (read_info.flag == 0 or read_info.flag == 16):
	if read_info.flag == '0' or read_info.flag == '16':
		record(read_info, g1, g2, READS_PICKLE)
		if read_info.name not in READS_SEQ_INFO.keys():
			READS_SEQ_INFO[read_info.name] = [[read_info.strand, read_info.seq, read_info.map_qua],[]]
		else:
			# recore_read in this name
			for supply_reads in READS_SEQ_INFO[read_info.name][1]:
				if supply_reads.strand == read_info.strand:
					supply_reads.seq = read_info.seq
					supply_reads.map_qua = read_info.map_qua
				else:
					supply_reads.seq = comp_reverse(read_info.seq)
					supply_reads.map_qua = "".join(list(reversed(read_info.map_qua)))
				record(supply_reads, g1, g2, READS_PICKLE)
				# recore
			READS_SEQ_INFO[read_info.name][0] = [ read_info.strand, read_info.seq, read_info.map_qua ]
			READS_SEQ_INFO[read_info.name][1] = []
	else:
		if read_info.name not in READS_SEQ_INFO.keys():
			READS_SEQ_INFO[read_info.name] = [ ['', '', ''], [read_info] ]
		else:
			if READS_SEQ_INFO[read_info.name][0][1] != '':
				if read_info.strand == READS_SEQ_INFO[read_info.name][0][0]:
					read_info.seq = READS_SEQ_INFO[read_info.name][0][1]
					read_info.map_qua = READS_SEQ_INFO[read_info.name][0][2]
				else:
					read_info.seq = comp_reverse(READS_SEQ_INFO[read_info.name][0][1])
					read_info.map_qua = "".join(list(reversed(READS_SEQ_INFO[read_info.name][0][2])))
				record(read_info, g1, g2, READS_PICKLE, 'ck3')
			else:
				READS_SEQ_INFO[read_info.name][1].append(read_info)

print(len(READS_SEQ_INFO.keys()))			
bamfile.close()
	

with open(outFile_pickle, "wb") as f1:
	pickle.dump(READS_PICKLE, f1)

g1.close()
g2.close()

endtime = datetime.datetime.now()
print (endtime - starttime)


