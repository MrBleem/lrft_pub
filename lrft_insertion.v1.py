
import datetime
starttime = datetime.datetime.now()
#long running

import os
import re
import sys
import pickle
import subprocess
import pysam
from operator import itemgetter, attrgetter
from collections import Counter


# imput
PROJECT_PATH = sys.argv[1]    # "/data/tusers/boxu/lrft/result/nanopore/human"
BAMFILE = sys.argv[2]        # "SRR11669560.mapped.sorted.tcg.bam"
PREFIX = sys.argv[3]        # "SRR11669560"


DATA_PATH = PROJECT_PATH + "/map/"
RESULT_PATH = PROJECT_PATH + "/validation/"
bamfile = pysam.AlignmentFile( DATA_PATH + BAMFILE, "rb" )
# PREFIX = bamFile.split('.')[0]

# filter_reads
# if sys.argv[4] :
    # CHROMS = sys.argv[4]
    # chrom = CHROMS.strip().split(',')
    # READS = bamfile.fetch( chrom[0], 10000, 1000000 )
# else:
# READS = bamfile

# READS = bamfile.fetch( "1" )
READS = bamfile.fetch( "1")
outFile = RESULT_PATH + PREFIX + ".lrft.table.chr1.txt"
outfile = open(outFile, 'w')
outfile.write('chr\tinsertion_start\tinsertion_end\tTE\tstrand\tsurpporting_reads\tun_surpporting_reads\tfrequency\tinsertion_type\tisnertion_te_length\tinsertion_te_sequence\n')


outFile_pickle = RESULT_PATH + PREFIX + ".lrft.chr1.cluster.pkl"
SEQ_INFO_FILE = DATA_PATH + PREFIX.split('.')[0] + ".clip.reads.pkl"

if os.path.exists(outFile_pickle):
    print('?')
    rm_args = ['rm', outFile_pickle]
    subprocess.Popen(rm_args)


# 多序列比对
def mafft(align_seq_file):
    ''' use MAFFT to create MSA '''

    out_fn = '.'.join(align_seq_file.split('.')[:-1]) + '.msa.fa'
    
    # 根据序列的长度进行参数的选择
    #     if 
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
        print(str(max_deletion)+"\t"+str(max2_deletion)+"\t"+str(clip_te_len) )
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
    cigar_temp = cigar.replace('H', 'S')
    s_cigar = re.findall('[0-9]*[A-Z]', cigar_temp.replace('X', 'D'))
    s_cigar_len = len(s_cigar)
    re_type = "novel"
    # 需要先判断比对到基因组上时是否还有clip的字段
    # 因为是从左到右增加，所以只要看左边是否clip就好了
    # print(s_cigar[0])
    if s_cigar[0][-1] == 'S' :
        cigar_len_S = int(s_cigar[0][:-1])
    else:
        cigar_len_S = 0
    cigar_len = 0
    s=''
    # 因为找insertion是按照从比对起点开始往后找的，所以理应先与clip的左边的reads长度进行比较
    # 如果左边的clip 长度是0，可能是原reads只测到那里
    if clip_left_len == 0 or cigar_len_S >= clip_left_len:
        # insertion_position_start = genome_map_position_start - clip_te_len
        # print('test')
        insertion_position_start = genome_map_position_start - 1
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
        
        cigar_len_dic = {'M':0, 'D':0, 'I':0, 'S':0}

        for i in range(s_cigar_len):
            s = s_cigar[i]
            cigar_len = cigar_len + int(s[:-1])

            cigar_len_dic[s[-1]] = cigar_len_dic[s[-1]] + int(s[:-1])
            cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']
            
            if cigar_len_ex_D >= clip_left_len:
                # 当延伸超过左边的长度时，需要判断下一个cigar是什么来判断insertion end在哪
                # 其实就是前面这段match的长度加deletions的长度
                if s[-1] == 'M':
                    insertion_position_start = genome_map_position_start + ( clip_left_len - cigar_len_dic['S'] + cigar_len_dic['D'] - cigar_len_dic['I'] ) -1
                    insertion_position_end = insertion_position_start + 1
                    if i < s_cigar_len-1:
                        
                        s_next = s_cigar[i+1]
                        print(s_next)
                        if abs(cigar_len_ex_D - clip_left_len) <= 10 and  s_next[-1] == 'D' : 
                            print('ck')
                            insertion_position_end = insertion_position_start + int(s_next[:-1])
                    q_position=cigar_len_ex_D

                elif s[-1] == 'I':
                    insertion_position_start = genome_map_position_start + ( cigar_len_dic['M'] + cigar_len_dic['D'] ) -1
                    # insertion_position_end = insertion_position_start + int(s[:-1])
                    insertion_position_end = insertion_position_start + 1
                    
                    q_position = cigar_len_ex_D - int(s[:-1])

                # elif s[-1] == 'D':
                # 不存在D的情况，因为碱基的累积不会停在缺失的地方
                    
                else:
                    if s[-1] == 'H' or s[-1] == 'S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        # insertion_position_end = insertion_position_start + clip_te_len + 1
                        insertion_position_end = insertion_position_start + 1
                        q_position = cigar_len_ex_D
                break
            else:
                if i < s_cigar_len-1:
                    s_next = s_cigar[i+1]
                    if s_next[-1]=='H' or s_next[-1]=='S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        # insertion_position_end = insertion_position_start + clip_te_len + 1
                        insertion_position_end = insertion_position_start + 1
                        q_position = cigar_len_ex_D
                        re_type="germline"
                        break
    return [insertion_position_start, insertion_position_end, ''.join(s_cigar[i-1:i+2]), q_position, re_type]



def get_insertion_position_germline(genome_map_position_start, cigar, clip_left_len):
    s_cigar = list(set(re.findall('[0-9]*D', cigar.replace('X', 'D'))))
    #   这两行是为了找到reads缺失最多的位置，通过观察，应该是前一步删除掉的转座子的部分
    s_cigar_S = list(set(re.findall('[0-9]*[HS]', cigar[:-1] )))
    if  s_cigar_S == "S" or s_cigar_S == "H":
        s_cigar_S = int(s_cigar_S[:-1])
    else:
        s_cigar_S = 0
    max_deletion = max([ int(t[:-1]) for t in s_cigar ])
    #   根据cigar字段确定insertion位置与reads比对起始位置之间的距离
    pre_insertion_cigar = cigar.replace('H', 'S').split(str(max_deletion)+'D')[0].split('S')[-1]
    pre_insertion_cigar_l = sum([ int(l[:-1])  for l in re.findall('[0-9]*[A-Z]', pre_insertion_cigar) ])
    pre_insertion_cigar_l_I = sum([ int(l[:-1])  for l in re.findall('[0-9]*[I]',pre_insertion_cigar) ])

    re_cigar2 = cigar.replace(str(max_deletion) + "D", str(max_deletion) + "M"  )
    s_cigar2 = list(set(re.findall('[0-9]*D',re_cigar2)))
    if len(s_cigar2) != 0 :
        max2_deletion = max([ int(t[:-1]) for t in  s_cigar2 ])
    else:
        max2_deletion = 0    
    pre_insertion_cigar2 = re_cigar2.replace('H', 'S').split(str(max2_deletion)+'D')[0].split('S')[-1]
    pre_insertion_cigar_l2 = sum([ int(l[:-1])  for l in re.findall('[0-9]*[A-Z]', pre_insertion_cigar2) ])
    pre_insertion_cigar_l_I2 = sum([ int(l[:-1])  for l in re.findall('[0-9]*[I]',pre_insertion_cigar2) ])

    
    # 确定insertion再基因组上的具体位置
    # pre_insertion_cigar_l_D = sum([ int(l[:-1])  for l in re.findall('[0-9]*[D]',pre_insertion_cigar) ])
    
    if abs( pre_insertion_cigar_l - ( clip_left_len - s_cigar_S ) ) > 10000:
        print('test')
        return get_insertion_position_germline(genome_map_position_start, re_cigar2, clip_left_len)
    else:
        insertion_position_start = genome_map_position_start + pre_insertion_cigar_l - pre_insertion_cigar_l_I
        insertion_position_end = insertion_position_start + max_deletion
        x1 = pre_insertion_cigar_l - pre_insertion_cigar_l_I
        x2 = pre_insertion_cigar_l2 - pre_insertion_cigar_l_I2
        if pre_insertion_cigar_l2 > pre_insertion_cigar_l and x2 - x1 - max_deletion  < 15 :
            # insertion_position_end = genome_map_position_start + pre_insertion_cigar_l - pre_insertion_cigar_l_I + max2_deletion
            insertion_position_end = genome_map_position_start + pre_insertion_cigar_l2 - pre_insertion_cigar_l_I2 + max2_deletion
            
        elif pre_insertion_cigar_l2 < pre_insertion_cigar_l  and x1 - x2 - max2_deletion  < 15:
            insertion_position_start = genome_map_position_start + pre_insertion_cigar_l2 - pre_insertion_cigar_l_I2
        
        return [insertion_position_start, insertion_position_end]


print('ck1')


def Cluster_raeds(READS):
    print('?')
    READS_CLUSTER={}
    for read in READS:
        # print(read.qname)
        # print(re.findall('[0-9]*[A-Z]',read.cigarstring)[0])
        # print(read.is_supplementary, read.is_secondary)
        # pysam提取reads的基本信息
        # 再前面的处理中，reads的名字中包含部分关于比对到TE上的有关信息
        # eg. SRR11669560.sra.9661265.ALU:ALU1.1.2831.245.9196
        # eg. SRR11669560.sra.12503931.ALU:ALU1.+.11.6487.269.2381
        if  read.is_supplementary or read.is_secondary:    # 竟然还留着这个，简直是智障 or read.flag==0: 
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
        
        if cigar_len_S_f >= clip_left_len or cigar_len_S_l >= clip_right_len:
            continue
        
        # 根据cliped reads 的reads 找insertion，并记录到字典中，
        # 这一步应该记录到bed文件中，可以用bedtools来处理，进行insertion区域的合并
        # 相当于是用intersect 来取代merge insertion
        cni =  check_novel_insertion(cigar, int(clip_te_len))
        insertion_type =cni[0]
        insertion_cigar = cni[1]
        insertion_position = get_insertion_position_novel(int(map_pos), int(clip_te_len), int(clip_left_len), int(clip_right_len), insertion_cigar)

        # if insertion_type == "novel":
        #     insertion_position = get_insertion_position_novel(int(map_pos), int(clip_te_len), int(clip_left_len), int(clip_right_len), insertion_cigar)
        #     insertion_type = insertion_position[-1]
        # else:
        #     insertion_position = get_insertion_position_germline(int(map_pos), insertion_cigar, int(clip_left_len))

        if TE not in READS_CLUSTER[ref]:
            READS_CLUSTER[ref][TE]=[]
            
        # 一条reads断点两端的序列信息记录
        cilp_seq = reads_id + "|" + TE_strand + "|" + genome_strand + "|" + seq[int(clip_left_len) - 50 : int(clip_left_len)] + "|" + seq[int(clip_left_len):int(clip_left_len)+50] 
    
        print(reads_id + "\t" + str(insertion_position[0]) +"\t" + str(insertion_position[1]) + "\t" + insertion_type)
        READS_CLUSTER[ref][TE].append([insertion_position[0], insertion_position[1], insertion_type, clip_te_len, cilp_seq])

    with open(outFile_pickle, "wb") as f1:
        pickle.dump(READS_CLUSTER, f1)
    
    return READS_CLUSTER


if os.path.exists(outFile_pickle):
    with open(outFile_pickle, "rb") as f2:
        READS_CLUSTER = pickle.load(f2)
    print('get reads cluster from pkl file')
else:
    READS_CLUSTER = Cluster_raeds(READS)


print('ck2')
TE_INSERTION_CLUSTER = {}
for ref in READS_CLUSTER:
    TE_CLUSTER=READS_CLUSTER[ref]
#     print(TE_CLUSTER[1])
    if ref not in TE_INSERTION_CLUSTER:
        TE_INSERTION_CLUSTER[ref] =  {}
    for te in TE_CLUSTER:
        TE_CLUSTER[te] = sorted(TE_CLUSTER[te])
        if len(TE_CLUSTER[te])==0 :
            continue
        TE_INSERTION_CLUSTER[ref][te] = []
        insertion_pre = TE_CLUSTER[te][0]
        supporting_reads = 1
        j=0
        flag = '1'
        for j in range(len(TE_CLUSTER[te][1:])):
            insertion = TE_CLUSTER[te][1:][j]
#             for insertion in TE_CLUSTER[te][1:]:
            # 以前一个insertion为参考，判断后一个insertion是否跟前一个有overlap
            # 这里判断的标准是看后一个insertion的中点是否在前一个insertion当中
            # 如果两个insertion有overlap，那么这两个insertion就合并在一起，区间合并在一起
#             if sum(insertion[0:2])/2 > insertion_pre[0] and sum(insertion[0:2])/2 < insertion_pre[1]:
            # insertion_pre[1] = insertion_pre[1] + 20    # 加入序列信息之后需要去掉 flag=delete
            if insertion[0] >= insertion_pre[0] and insertion[0] <= insertion_pre[1]+20:
                insertion_pre = merge_insertion(insertion_pre, insertion)
                supporting_reads = supporting_reads + 1
            else:
#                mapped_reads是根据insertion的位置，再用pysam计算insertion的位置上有多少条reads比对到了这个位置
                mapped_reads = bamfile.count(contig=ref, start=insertion_pre[0], stop=insertion_pre[1]+500)
                # insertion_start, insertion_end, insertion_type, clip_te_len, cilp_seq
                TE_INSERTION_CLUSTER[ref][te].append([insertion_pre[0], insertion_pre[1], supporting_reads, mapped_reads, insertion_pre[2], insertion_pre[3], insertion_pre[4] ])
#                 TE_INSERTION_CLUSTER[te].append(insertion_pre)
                insertion_pre = insertion
                supporting_reads = 1
                if j == len(TE_CLUSTER[te][1:]):
                    flag='end'
        if flag !='end' :
            mapped_reads = bamfile.count(contig=ref, start=insertion_pre[0], stop=insertion_pre[1]+500)
            TE_INSERTION_CLUSTER[ref][te].append([insertion_pre[0], insertion_pre[1], supporting_reads, mapped_reads, insertion_pre[2], insertion_pre[3], insertion_pre[4]])


with open(SEQ_INFO_FILE, "rb") as f3:
    READS_SEQUENCE = pickle.load(f3)
    
    
print('ck3')
for ref in TE_INSERTION_CLUSTER:
    for te in TE_INSERTION_CLUSTER[ref]:
        for insertion in TE_INSERTION_CLUSTER[ref][te]:
    #         chr  insertion_start  insertion_end  TE  TE_start  TE_end  surpporting_reads  un_surpporting_reads  fre  
            # 仍然存在在genome中找到对应insertion位置的
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
            # multi sequence alignment
            align_sequences = {}
            id = None
            seq = ''
            for line in open(mafft( RESULT_PATH + PREFIX + ".seq.temp.fq" ), 'r'):
                if line.startswith('>'):
                    if id != None:
                        align_sequences[id] = seq.upper()
                    id = line.strip().split('>')[1]
                    seq = ''
                else:
                    seq = seq + line.strip()
            align_sequences[id] = seq.upper()
            align_seq = '\t|\t'.join(get_consensus_seq(align_sequences))
            
            # 标记检测，测试用的
            if len(Counter(insertion_strand).items()) == 1 and list(Counter(insertion_strand).items())[0][1] != 1 :
                insertion_strand.append('CAO')
            elif len(Counter(insertion_strand).items()) == 2:
                insertion_strand.append('GAN')
            
            strand_info = ['_'.join(insertion_strand_te), '_'.join(insertion_strand_genome), '_'.join(insertion_strand)] 
            # outfile.write(str(ref)+'\t'+str(insertion[0])+'\t'+str(insertion[1])+'\t'+te+'\t'+'|'.join(strand_info)+'\t'+ str(insertion[2])+'\t'+str(insertion[3]-insertion[2])+'\t'+str(float(insertion[2])/insertion[3])+'\t'+insertion[4]+'\t'+str(insertion[5]) + '\t' + align_seq + '\t'+'_'.join(READS_SEQ)+'\n')    
            insertion_final = [str(ref), str(insertion[0]), str(insertion[1]), te, '|'.join(strand_info), str(insertion[2]), str(insertion[3]-insertion[2]), str(float(insertion[2])/insertion[3]), insertion[4], str(insertion[5]), align_seq, '_'.join(READS_SEQ)]
            outfile.write('\t'.join(insertion_final) + '\n')    

        
outfile.close()        
print('>_< DONE >V<')
endtime = datetime.datetime.now()
print (endtime - starttime)

