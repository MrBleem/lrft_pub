

import os
import re
import sys
import pickle
import pysam
from operator import itemgetter, attrgetter
from collections import Counter


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

# READS = bamfile.fetch( "1", 10000, 1000000 )

outFile_pickle = RESULT_PATH + PREFIX + ".lrft.cluster.pkl"
SEQ_INFO_FILE = DATA_PATH + PREFIX.split('.')[0] + ".clip.reads.pkl"




# 对序列进行反向互补
def comp_reverse(seq):
    base = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    new_seq = ''
    for i in range(len(seq)):
        s = seq[i]
        new_seq = base[s] + new_seq
    return new_seq


 
    
# 判断这个insertion的reads上面是否有一个很大的deletion
# 如果有，就可能是一个germline insertion
# 如果没有，就可能是一个novel insertion   
def check_novel_insertion(cigar, clip_te_len):
    # ----------------------------------------   reference
    # --------------        ------------------   reads
    #              max deletion
    Type='novel'
    re_cigar = cigar
    s_cigar = list(set(re.findall('[0-9]*D', cigar)))
    if len(s_cigar) != 0 :
        max_deletion = max([ int(t[:-1]) for t in s_cigar ])
        # ? 如果两个deletion离的比较远怎么判断呢
        # 这种情况是不是应该只取最大的deletion周围的cigar来判断
        re_cigar2 = cigar.replace(str(max_deletion) + "D", str(max_deletion) + "M"  )
        s_cigar2 = list(set(re.findall('[0-9]*D',re_cigar2)))
        if len(s_cigar2) != 0 :
            max2_deletion = max([ int(t[:-1]) for t in  s_cigar2 ])
        else:
            max2_deletion = 0
        # 比较deletion的长度跟插入的转座子长度，如果相差较大，可能是比对时断点周围的序列也可以比对到这个缺口的地方，使大的insertion变成两段或多段
        if max_deletion >= clip_te_len/2 or max_deletion + max2_deletion >= clip_te_len/2 :
            if abs(max_deletion-clip_te_len) <= clip_te_len/2  :
                Type="germline"
                return [Type,cigar]
            else:
                re_cigar = cigar.replace(str(max_deletion)+"D", str(max_deletion)+"X"  )
                return check_novel_insertion(re_cigar, clip_te_len)
        else:
            Type="novel"
            return [Type,re_cigar]
    else:
        return [Type,cigar]
    
# 通过结果来看，这个地方好像是更针对于结构比较清晰的 insertion
def get_insertion_position_novel(genome_map_position_start, clip_te_len, clip_left_len, clip_right_len, cigar):
    s_cigar = re.findall('[0-9]*[A-Z]',cigar.replace('X', 'D'))
    s_cigar_len = len(s_cigar)
    re_type = "novel"
    # 需要先判断比对到基因组上时是否还有clip的字段
    # 因为是从左到右增加，所以只要看左边是否clip就好了
    if s_cigar[0][-1] == 'S'  or s_cigar[0][-1] == 'H':
        cigar_len_S = int(s_cigar[0][:-1])
    else:
        cigar_len_S = 0
    cigar_len = 0
    s=''
    # 因为找insertion是按照从比对起点开始往后找的，所以理应先与clip的左边的reads长度进行比较
    # 如果左边的clip 长度是0，可能是原reads只测到那里
    if clip_left_len == 0 or cigar_len_S >= clip_left_len:
        insertion_position_start = genome_map_position_start - clip_te_len
        if insertion_position_start <=0 :
            insertion_position_start = 1
        insertion_position_end = genome_map_position_start
        q_position = 0
        re_type = "germline"
        i=1
    else:
        # cigar_len_S <= clip_left_len:
        # 计算cigar的长度
        # 满足条件 去除Deletion的部分长度达到clip的长度停止计算

        # Deletion 不算在clip的长度增加中
        # insertion是算在clip的长度增加中的
        # (mis)match也是算在clip的长度增加中的

        # ---------   --------    ---｜---------     ---------   reference
        # ------   ---------  -------｜-        --------------   cliped reads
        
        cigar_len_dic = {'M':0, 'D':0, 'I':0, 'S':cigar_len_S}

        for i in range(s_cigar_len):
            s = s_cigar[i]
            cigar_len = cigar_len + int(s[:-1])

            cigar_len_dic[s[-1]] = cigar_len_dic[s[-1]] + int(s[:-1])
            cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']
            
            if cigar_len_ex_D >= clip_left_len:
                # 当延伸超过左边的长度时，需要判断下一个cigar是什么来判断insertion end在哪
                # 其实就是前面这段match的长度加deletions的长度
                if s[-1] == 'M':
                    insertion_position_start = genome_map_position_start + ( clip_left_len - cigar_len_dic['S'] + cigar_len_dic['D'] - cigar_len_dic['I'] ) - 1
                    insertion_position_end = insertion_position_start + 1
                    if i < s_cigar_len-1:
                        s_next = s_cigar[i+1]
                        if cigar_len_ex_D == clip_left_len and  s_next[-1] == 'D' : 
                            insertion_position_end = insertion_position_start + int(s_next[:-1])
                    q_position=cigar_len_ex_D

                elif s[-1] == 'I':
                    insertion_position_start = genome_map_position_start + ( cigar_len_dic['M'] + cigar_len_dic['D'] ) - 1
                    insertion_position_end = insertion_position_start + int(s[:-1])
                    q_position = cigar_len_ex_D - int(s[:-1])
                    
                else:
                    if s[-1] == 'H' or s[-1] == 'S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D']  - 1
                        insertion_position_end = insertion_position_start + clip_te_len + 1
                        q_position = cigar_len_ex_D
                break
            else:
                if i < s_cigar_len-1:
                    s_next = s_cigar[i+1]
                    if s_next[-1]=='H' or s_next[-1]=='S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        insertion_position_end = insertion_position_start + clip_te_len + 1
                        q_position = cigar_len_ex_D
                        re_type="germline"
                        break
    return [insertion_position_start, insertion_position_end, ''.join(s_cigar[i-1:i+2]), q_position, re_type]



def get_insertion_position_germline(genome_map_position_start, cigar):
    # 这两行是为了找到reads缺失最多的位置，通过观察，应该是前一步删除掉的转座子的部分
    s_cigar = list(set(re.findall('[0-9]*D', cigar.replace('X', 'D'))))
    max_deletion = max([ int(t[:-1]) for t in s_cigar ])
    # 根据cigar字段确定insertion位置与reads比对起始位置之间的距离
    pre_insertion_cigar = cigar.split(str(max_deletion)+'D')[0].split('S')[-1]
    pre_insertion_cigar_l = sum([ int(l[:-1])  for l in re.findall('[0-9]*[A-Z]', pre_insertion_cigar) ])
    # 确定insertion再基因组上的具体位置
    # pre_insertion_cigar_l_D = sum([ int(l[:-1])  for l in re.findall('[0-9]*[D]',pre_insertion_cigar) ])
    pre_insertion_cigar_l_I = sum([ int(l[:-1])  for l in re.findall('[0-9]*[I]',pre_insertion_cigar) ])
    insertion_position_start = genome_map_position_start + pre_insertion_cigar_l - pre_insertion_cigar_l_I
    insertion_position_end = insertion_position_start + max_deletion
    
    return [insertion_position_start, insertion_position_end]


print('ck1')


def Cluster_raeds(READS):
    print('?')
    READS_CLUSTER={}
    for read in READS:
        # pysam提取reads的基本信息
        # 再前面的处理中，reads的名字中包含部分关于比对到TE上的有关信息
        # eg. SRR11669560.sra.9661265.ALU:ALU1.1.2831.245.9196
        # eg. SRR11669560.sra.12503931.ALU:ALU1.+.11.6487.269.2381
        if read.is_secondary or read.is_supplementary:    # 竟然还留着这个，简直是智障 or read.flag==0: 
            continue
        name = read.qname.split(',') # reads name的分隔符
        ref = read.reference_name
        seq = read.query_sequence
        if ref not in READS_CLUSTER:
        # READS_CLUSTER[ref] =  {'ALU':[],'LINE1':[], 'SVA':[]}
            READS_CLUSTER[ref] =  {}

        # reads_id = name[2]
        reads_id = read.qname

        # 50a9e472-a887-4803-a02c-86f17a2da755.Ko.L1.22.-.1.11804.2797.1424
        # SRR11669560.sra.13891323. LINE1:LINE1. -.2269.5481.1501.13704
        
        clip_left_len = name[4]
        clip_te_len = int(name[5])
        clip_right_len = name[6]
        TE = name[1].split(':')[0]
        map_pos = read.pos + 1
        cigar = read.cigarstring
        TE_strand = name[2]

        
        
        # 判断reads比对的方向，方向会影响前一步cilp的左右方向
        flag=hex(read.flag)
        flag_strand = int(flag[2:])%100
        if flag_strand != 10 :
            genome_strand = '+'
        else :
            genome_strand = '-'
            t = clip_right_len
            clip_right_len = clip_left_len
            clip_left_len = int(t)
        
        # 根据cliped reads 的reads 找insertion，并记录到字典中，
        # 这一步应该记录到bed文件中，可以用bedtools来处理，进行insertion区域的合并
        # 相当于是用intersect 来取代merge insertion
        insertion_type = check_novel_insertion(cigar, int(clip_te_len))[0]
        insertion_cigar = check_novel_insertion(cigar, int(clip_te_len))[1]
        if insertion_type == "novel":
            insertion_position = get_insertion_position_novel(int(map_pos), int(clip_te_len), int(clip_left_len), int(clip_right_len), insertion_cigar)
            insertion_type = insertion_position[-1]
        else:
            insertion_position = get_insertion_position_germline(int(map_pos), insertion_cigar)

        if TE not in READS_CLUSTER[ref]:
            READS_CLUSTER[ref][TE]=[]
            
        # 一条reads断点两端的序列信息记录
        cilp_seq = reads_id + "|" + TE_strand + "|" + genome_strand + "|" + seq[int(clip_left_len) - 50 : int(clip_left_len)] + "|" + seq[int(clip_left_len):int(clip_left_len)+50] 
    
        with open(RESULT_PATH + "tmp/SRR11669560." + ref + "." + TE + ".cluster.bed", 'w') as cluster_bed:
            cluster_bed.write( "\t".join([ ref, str(insertion_position[0]), str(insertion_position[1]),  str(clip_te_len), TE, cilp_seq + "," + insertion_type ]) + "\n")

        # READS_CLUSTER[ref][TE].append([insertion_position[0], insertion_position[1], insertion_type, clip_te_len, cilp_seq])

    # with open(outFile_pickle, "wb") as f1:
    #     pickle.dump(READS_CLUSTER, f1)
    
    # return READS_CLUSTER


Cluster_raeds(READS)
print('>_< Cluster DONE >V<')

# if os.path.exists(outFile_pickle):
#     with open(outFile_pickle, "rb") as f2:
#         READS_CLUSTER = pickle.load(f2)
#     # print('get reads cluster from pkl file')
#     print('already cluster')
# else:
#     READS_CLUSTER = Cluster_raeds(READS)

  


