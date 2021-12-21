

import os
import re
import sys
import pickle
import pysam
from operator import itemgetter, attrgetter
from collections import Counter

from lrft_get_insertion_position import get_insertion_position

def Cluster_raeds(READS, outFile_pickle):
    READS_CLUSTER={}
    for read in READS:
        
        # print(re.findall('[0-9]*[A-Z]',read.cigarstring)[0])
        # print(read.is_supplementary, read.is_secondary)
        # pysam提取reads的基本信息
        # 再前面的处理中，reads的名字中包含部分关于比对到TE上的有关信息
        # eg. SRR11669560.sra.9661265.ALU:ALU1.1.2831.245.9196
        # eg. SRR11669560.sra.12503931.ALU:ALU1.+.11.6487.269.2381
        if  read.is_supplementary:
            print(read.qname + ":type:supplementary")
            continue
        
        if read.is_secondary:    # 竟然还留着这个，简直是智障 or read.flag==0: 
            print(read.qname + ":type:secondary")
            if read.seq == None or read.cigar == None:
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
        
        clip_left_len = int(name[4])
        clip_te_len = int(name[5])
        clip_right_len = int(name[6])
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
        
        # 如果remap后两边clip的长度太长超过了 左边或者右边截断的reads，就意味着这个insertion没有跨过这个insertion，那这个reads就不能采用
        # 再严格一点就要限制跨过insertion的长度要控制在 N bp 以内
        s_cigar = re.findall('[0-9]*[A-Z]',cigar)
        if s_cigar[0][-1] == 'S'  or s_cigar[0][-1] == 'H':
            cigar_len_S_f = int(s_cigar[0][:-1])
        else:
            cigar_len_S_f = 0
        if s_cigar[-1][-1] == 'S'  or s_cigar[-1][-1] == 'H':
            cigar_len_S_l = int(s_cigar[-1][:-1])
        else:
            cigar_len_S_l = 0
        
        

        if cigar_len_S_f > clip_left_len or cigar_len_S_l > clip_right_len:
            print(read.qname + ":type:clipped")
        else:
            print(read.qname + ":type:perfect")
        
        if cigar_len_S_f > clip_left_len+30 or cigar_len_S_l > clip_right_len+30:
            continue
            
        
        # 根据cliped reads 的reads 找insertion，并记录到字典中，
        # 这一步应该记录到bed文件中，可以用bedtools来处理，进行insertion区域的合并
        # 相当于是用intersect 来取代merge insertion

        # insertion_position = get_insertion_position(int(map_pos), int(clip_te_len), int(clip_left_len), int(clip_right_len), insertion_cigar)
        insertion_position = get_insertion_position(int(map_pos), int(clip_te_len), int(clip_left_len), int(clip_right_len), cigar)

        if TE not in READS_CLUSTER[ref]:
            READS_CLUSTER[ref][TE]=[]
            
        # 一条reads断点两端的序列信息记录
        query_position = insertion_position[3]
        # cilp_seq = reads_id + "|" + TE_strand + "|" + genome_strand + "|" + seq[int(clip_left_len) - 50 : int(clip_left_len)] + "|" + seq[int(clip_left_len):int(clip_left_len)+50] 
        cilp_seq = reads_id + "|" + \
                    TE_strand + "|" + \
                    genome_strand + "|" + \
                    seq[ int(query_position[0]) - 50 : int(query_position[0]) ] + "|" + \
                    seq[ int(query_position[0]) : int(query_position[1]) ]  + "|" + \
                    seq[ int(query_position[2]) : int(query_position[3]) ]  + "|" + \
                    seq[ int(query_position[3]) : int(query_position[3] + 50) ]
                    # seq[ int(query_position[1]) : int(query_position[2]) ]  + "|" + \
    
        # print(reads_id + "\t" + str(insertion_position[0]) +"\t" + str(insertion_position[1]) + "\t" + insertion_type)
        if insertion_position[1] - insertion_position[0] >= 20 : 
            insertion_type = "germline"
        else:
            insertion_type = "novel"
        READS_CLUSTER[ref][TE].append([insertion_position[0], insertion_position[1], insertion_type, clip_te_len, cilp_seq])

    with open(outFile_pickle, "wb") as f1:
        pickle.dump(READS_CLUSTER, f1)
    
    return READS_CLUSTER


