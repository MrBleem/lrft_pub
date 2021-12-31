
####

# 比v1添加了insertion类型的判断，暂时将insertion的长度作为判断标准 

# 比v2添加了trunked的判断，就是截断reads的insertion的判断

# 比v3:尝试全部模块化

####



import datetime
starttime = datetime.datetime.now()
#long running

import os
import re
import sys
import pickle
import getopt
import subprocess
import pysam
from operator import itemgetter, attrgetter
from collections import Counter


from lrft_consensus_seq import consensus_seq
# from lrft_get_insertion_position import get_insertion_position
from lrft_cluster import Cluster_raeds

# help文档
# if __name__ == '__main__':
#     helpdoc = '''
# Descrtption:
#     long reads for transposon

# Usage
#     python lrft_insertion -i path_to_bam_file -o out_path -p prefix

# Parameter
#     -h/--help
#         Print helpdoc
#     -i/--infile
#         Input file
#     -o/--outpath
#         Out path
#     -p/--prefix
#         Result prefix
#     '''
#     try:
#         opts,args = getopt.getopt(sys.argv[1:], "hi:o:p:", ["infile=", "outpath=", "prefix="])
#         if len(opts) != 2:
#             print( "Options Error!\n\n" + helpdoc)
#             sys.exit(2)
#     except getopt.GetoptError:
#         print( "Options Error!\n\n" + helpdoc )
#         sys.exit(2)
#     for opt,arg in opts:
#         if opt in ("-h","--help"):
#             print(helpdoc)
#             sys.exit()
#         elif opt in ("-i","--infile"):
#             BAMFILE = arg
#             PREFIX = BAMFILE.split('/')[-1].split('.bam')[0]
#         elif opt in ("-o","--outpath"):
#             OUTPATH = arg
#         elif opt in ("-p", "--prefix"):
#             PREFIX = arg


# 这个版本只用一种insertion的方法， max_deletion的方方法暂时不用
# 因为想在大部分物种中使用，主要依靠一种思想比较靠谱

# imput
PROJECT_PATH = sys.argv[1]    # "/data/tusers/boxu/lrft/result/nanopore/human"
BAMFILE = sys.argv[2]        # "SRR11669560.mapped.sorted.tcg.bam"
PREFIX = sys.argv[3]        # "SRR11669560"


DATA_PATH = PROJECT_PATH + "/map/"
RESULT_PATH = PROJECT_PATH + "/insertion/"
bamfile = pysam.AlignmentFile( DATA_PATH + BAMFILE, "rb" )
# PREFIX = bamFile.split('.')[0]

# filter_reads
# if sys.argv[4] :
    # CHROMS = sys.argv[4]
    # chrom = CHROMS.strip().split(',')
    # READS = bamfile.fetch( chrom[0], 10000, 1000000 )
# else:
READS = bamfile

# READS = bamfile.fetch( "1" )
# region="1:100,993,787-100,995,464"
# chr=region.split(':')[0]
# region_start=region.split(':')[1].split('-')[0].replace(',', '')
# region_end=region.split(':')[1].split('-')[1].replace(',', '')
# READS = bamfile.fetch( chr, int(region_start) , int(region_end))

outFile = RESULT_PATH + PREFIX + ".lrft.table.txt"
outfile = open(outFile, 'w')
outfile.write('chr\tinsertion_start\tinsertion_end\tTE\tstrand\tsurpporting_reads\tun_surpporting_reads\tfrequency\tinsertion_type\tisnertion_te_length\tinsertion_te_sequence\n')


outFile_pickle = RESULT_PATH + PREFIX + ".lrft.cluster.pkl"
SEQ_INFO_FILE = DATA_PATH + PREFIX.split('.')[0] + ".clip.reads.pkl"


# if os.path.exists(outFile_pickle):
#     print('?')
# rm_args = ['rm', outFile_pickle]
# subprocess.Popen(rm_args)


# insertion_group
insertion_types = [ 'perfect', 'left_clipped', 'right_clipped', 'secondary_alignment', 'supplymentary_alignment' ]


# 对序列进行反向互补
def comp_reverse(seq):
    base = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    new_seq = ''
    for i in range(len(seq)):
        s = seq[i]
        new_seq = base[s] + new_seq
    return new_seq


# 合并两个insertion区间，取并集
def merge_insertion(insertion1, insertion2):
    # --------      insertion1
    #    --------   insertion2
    # -----------   merged_insertion
    new_insertion_s = insertion1[0]
    new_insertion_e = insertion1[1]
    # if insertion2[0] < new_insertion_s:
    #     new_insertion_s = insertion2[0]
    if insertion2[1] > new_insertion_e:
        new_insertion_e = insertion2[1]
        
    return [new_insertion_s, new_insertion_e, insertion1[2]+";"+insertion2[2], sum([insertion1[3],insertion2[3]])/2,  insertion1[4]+";"+insertion2[4] ]


def record_insertion(ref, te, insertion):
    # 将insertion记录到文件中
    # 这一步包括：
    #       · 判断insertion的方向 
    #       · insertion的consensus sequence
    if insertion[3] == 0:
        insertion[3] = insertion[2]
    #                 outfile_confused.write(str(ref)+'\t'+str(insertion[0])+'\t'+str(insertion[1])+'\t'+te+'\t*\t'+ str(insertion[2])+'\t'+str(insertion[3]-insertion[2])+'\t'+str(float(insertion[2])/insertion[3])+'\t'+insertion[4]+'\t'+str(insertion[5])+'\t'+str(insertion[6])+'\n')
    READS_SEQ_INFO = insertion[6]
    READS_SEQ = []
    insertion_strand = []
    insertion_strand_te = []
    insertion_strand_genome = []
    align_seq_file = open( RESULT_PATH + PREFIX +".seq.temp.fq", 'w' )

    for read_seq in READS_SEQ_INFO.split(';'):                
        read_seq_info = read_seq.split('|')
                
        read_id = read_seq_info[0]
        te_strand = read_seq_info[1]
        genome_strand = read_seq_info[2]
        clip_left = read_seq_info[3]
        clip_right = read_seq_info[4]

        insertion_strand_te.append(te_strand)
        insertion_strand_genome.append(genome_strand)
                
        # bed中记录的TE方向都是正的，正负代表与序列的相对方向
        # insertion的方向完全取决于clip reads 比对到genome上的方向
        # 要提取的reads也是只跟比对到genome上的情况对应，与第一步的方向没有关系
        # 因为我要的是seq2上｜｜之间的seq，所以根据seq3的reads以及方向就有可以得到
        # raw reads
        # >>>>>>>>>>>>>>>>>>>///////////////>>>>>>>>>A>>>>>>>>>> seq1
        # >>>>>>>>>>>>>>>>>>>\\\\\\\\\\\\\\\>>>>>>>>>A>>>>>>>>>> 
        # map to TE:
        # >>>>>>>>>>>>>>>｜>>>>               >>>>？>>>>>A>>>>>>>>>> seq2    ///////////////   TE
        # <<<<<<T<<<<<<<<｜<<<<               <<<<｜<<<<<<<<<<<<<<<<         ///////////////

        # clip reads map to genome
        # <<<<<<T<<<<<<<<？<<<<               <<<<｜<<<<<<<<<<<<<<<< seq3     
        # >>>>>>>>>>>>>>>｜>>>>               >>>>｜>>>>>A>>>>>>>>>>

        cilp_mid = READS_SEQUENCE[read_id]
        if genome_strand == "+" :
            insertion_strand.append("+")
            READS_SEQ.append(read_id + "|" + "\033[34m" + clip_left + "\033[0m" + cilp_mid + "\033[34m" + clip_right + "\033[0m")
            align_seq_file.write( '>' + read_id + '\n' )
            align_seq_file.write( clip_left + cilp_mid + clip_right + '\n' )
        else :
            insertion_strand.append("-")
            READS_SEQ.append( read_id + "|" + "\033[34m" + comp_reverse(clip_right) + "\033[0m" + cilp_mid + "\033[34m" + comp_reverse(clip_left) + "\033[0m" )
            align_seq_file.write( '>' + read_id + '\n' )
            align_seq_file.write( comp_reverse(clip_right) + cilp_mid + comp_reverse(clip_left) + '\n' )

    align_seq_file.close()
    fre_cuoff = 0.8
    align_seq = consensus_seq(RESULT_PATH + PREFIX + ".seq.temp.fq", fre_cuoff)
    align_seq = ''
            
    # 标记检测，测试用的
    if len(Counter(insertion_strand).items()) == 1 and list(Counter(insertion_strand).items())[0][1] != 1 :
        insertion_strand.append('CAO')
    elif len(Counter(insertion_strand).items()) == 2:
        insertion_strand.append('GAN')
                

    strand_info = ['_'.join(insertion_strand_te), '_'.join(insertion_strand_genome), '_'.join(insertion_strand)] 
    # outfile.write(str(ref)+'\t'+str(insertion[0])+'\t'+str(insertion[1])+'\t'+te+'\t'+'|'.join(strand_info)+'\t'+ str(insertion[2])+'\t'+str(insertion[3]-insertion[2])+'\t'+str(float(insertion[2])/insertion[3])+'\t'+insertion[4]+'\t'+str(insertion[5]) + '\t' + align_seq + '\t'+'_'.join(READS_SEQ)+'\n')    
    insertion_final = [str(ref), str(insertion[0]), str(insertion[1]), te, '|'.join(strand_info), str(insertion[2]), str(insertion[3]-insertion[2]), str(float(insertion[2])/insertion[3]), insertion[4], str(insertion[5]), align_seq, '_'.join(READS_SEQ)]
    outfile.write('\t'.join(insertion_final) + '\n')


# 加载前一步的文件
if os.path.exists(outFile_pickle):
    with open(outFile_pickle, "rb") as f2:
        READS_CLUSTER = pickle.load(f2)
    print('get reads cluster from pkl file')
else:
    READS_CLUSTER = Cluster_raeds(READS, outFile_pickle)

with open(SEQ_INFO_FILE, "rb") as f3:
    READS_SEQUENCE = pickle.load(f3)


# 遍历insertion，把临近的insertion merge到一起，这时候还没有考虑不同方向，不同类型的insertion
# 遇到有特殊情况的需要标注开来

# 一级遍历 ： ref
#   二级遍历：ref - te
#       三级遍历：ref - te - insertion（sorted）
for ref in READS_CLUSTER:
    TE_CLUSTER = READS_CLUSTER[ref]

    for te in TE_CLUSTER:
        TE_CLUSTER[te] = sorted(TE_CLUSTER[te])

        if len(TE_CLUSTER[te])==0 :
            continue

        insertion_pre = TE_CLUSTER[te][0]
        supporting_reads = 1
        # j = 0
        flag = '1'
        for j in range(len(TE_CLUSTER[te][1:])):
            insertion = TE_CLUSTER[te][1:][j]
            # 以前一个insertion为参考，判断后一个insertion是否跟前一个有overlap
            # 这里判断的标准是看后一个insertion的中点是否在前一个insertion当中
            # 如果两个insertion有overlap，那么这两个insertion就合并在一起，区间合并在一起

            # if sum(insertion[0:2])/2 > insertion_pre[0] and sum(insertion[0:2])/2 < insertion_pre[1]:
            # insertion_pre[1] = insertion_pre[1] + 20    # 加入序列信息之后需要去掉 flag=delete

            if insertion[0] >= insertion_pre[0] and insertion[0] <= insertion_pre[1]+20:
                insertion_pre = merge_insertion(insertion_pre, insertion)
                supporting_reads = supporting_reads + 1
            else:
                #   mapped_reads是根据insertion的位置，再用pysam计算insertion的位置上有多少条reads比对到了这个位置
                mapped_reads = bamfile.count(contig=ref, start=insertion_pre[0], stop=insertion_pre[1]+500)
                # insertion_start, insertion_end, insertion_type, clip_te_len, cilp_seq
                insertion_record = [insertion_pre[0], insertion_pre[1], supporting_reads, mapped_reads, insertion_pre[2], insertion_pre[3], insertion_pre[4] ]
                record_insertion(ref, te, insertion_record)
                insertion_pre = insertion
                supporting_reads = 1
                if j == len(TE_CLUSTER[te][1:]):
                    flag='end'
        if flag !='end' :
            mapped_reads = bamfile.count(contig=ref, start=insertion_pre[0], stop=insertion_pre[1]+500)
            insertion_record = [insertion_pre[0], insertion_pre[1], supporting_reads, mapped_reads, insertion_pre[2], insertion_pre[3], insertion_pre[4]]
            record_insertion(ref, te, insertion_record)

outfile.close()        
print('>_< DONE >V<')
endtime = datetime.datetime.now()
print (endtime - starttime)
