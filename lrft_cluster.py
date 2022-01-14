



import re
import pickle
import pysam

from lrft_get_insertion_position import get_insertion_position

def Cluster_reads(READS, READS_SEQUENCE, outFile_pickle):
    READS_CLUSTER={}
    for read in READS:
        # pysam提取reads的基本信息
        # 再前面的处理中，reads的名字中包含部分关于比对到TE上的有关信息
        # eg. SRR11669560.sra.9661265.ALU:ALU1.1.2831.245.9196
        # eg. SRR11669560.sra.12503931.ALU:ALU1.+.11.6487.269.2381
        if  read.is_supplementary:
            # print(read.qname + ":type:supplementary")
            continue
        
        if read.is_secondary:   
            # print(read.qname + ":type:secondary")
            if read.query_sequence == None or read.cigar == None:
                continue  
        name = read.qname.split(',') # reads name的分隔符
        ref = read.reference_name
        seq = read.query_sequence
        if ref not in READS_CLUSTER:
            READS_CLUSTER[ref] =  {}
        reads_id = read.qname
        
        clip_left_len = int(name[4])
        clip_te_len = int(name[5])
        clip_right_len = int(name[6])
        TE = name[1].split(':')[0]
        map_pos = read.pos + 1
        cigar = read.cigarstring
        TE_strand = name[2]
        
        # 判断reads比对的方向，方向会影响前一步cilp的左右方向
        flag = hex(read.flag)
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

        cigar_list = re.findall('[0-9]*[A-Z]',cigar)
        if cigar_list[0][-1] == 'S'  or cigar_list[0][-1] == 'H':
            cigar_len_S_start  = int(cigar_list[0][:-1])
        else:
            cigar_len_S_start = 0

        if cigar_list[-1][-1] == 'S'  or cigar_list[-1][-1] == 'H':
            cigar_len_S_end = int(cigar_list[-1][:-1])
        else:
            cigar_len_S_end = 0
        
        if cigar_len_S_start > clip_left_len + 60 or cigar_len_S_end > clip_right_len + 60:
            print(reads_id, cigar_len_S_start, cigar_len_S_end)
            continue

        # reads 要跨过insertion position 30 bp
        # spanding_length
        if  clip_left_len - cigar_len_S_start <= 30 or   clip_right_len - cigar_len_S_end <= 30:
            insertion_tag = "truncked"
        else:
            insertion_tag = 'spaned'
            
        
            
        
        # 根据cliped reads 的reads 找insertion，并记录到字典中，
        # 这一步应该记录到bed文件中，可以用bedtools来处理，进行insertion区域的合并
        # 相当于是用intersect 来取代merge insertion

        insertion_position = get_insertion_position(int(map_pos), int(clip_te_len), int(clip_left_len), int(clip_right_len), cigar)
        # print(reads_id, insertion_position[0], insertion_position[1])
        if TE not in READS_CLUSTER[ref]:
            READS_CLUSTER[ref][TE] = []
            
        # 一条reads断点两端的序列信息记录
        # 通过断点周围的信息，
        query_position = insertion_position[3]
        cilp_te = READS_SEQUENCE[reads_id]   # clipped TE sequence

        # ref_seq = genome_fa_seq.fetch(ref,int(insertion_position[0])-15 , int(insertion_position[1])+14 )
        
        # if int(insertion_position[1]) - int(insertion_position[0]) <= 15:
        #     candinate_TSD_genome = ref_seq
        # else:
        #     candinate_TSD_genome = "".join([ref_seq[0:20], ref_seq[-20:]])
        
        # query_position[0] 是左边clip reads比对到基因组的那个点，并不一定是clip 的那个点
        # 对应的query_position[3] 是右边的那个点

        if query_position[0] == query_position[1] :
            candinate_TSD_left = ''
            candinate_TSD_right = seq[ int(query_position[2]) - 15 : int(query_position[2]) ] + seq[ int(query_position[2]) : int(query_position[2] + 15) ]
        elif query_position[1] == query_position[2]:
            candinate_TSD_left = seq[ int(query_position[0]) - 15 : int(query_position[0]) ] + seq[ int(query_position[0]) : int(query_position[0]) + 15 ]
            candinate_TSD_right = ''
        else:
            if query_position[1] - query_position[0] < 15:
                candinate_TSD_left = seq[ int(query_position[0]) - 15 : int(query_position[0]) ] + seq[ int(query_position[0]) : int(query_position[1]) ]
            else:
                candinate_TSD_left = seq[ int(query_position[0]) - 15 : int(query_position[0]) ] + seq[ int(query_position[0]) : int(query_position[0]) + 15 ]

            if query_position[2] - query_position[1] < 15:
                candinate_TSD_right = seq[ int(query_position[1]) : int(query_position[2]) ] + seq[ int(query_position[2]) : int(query_position[2] + 15) ]
            else:
                candinate_TSD_right = seq[ int(query_position[2]) - 15 : int(query_position[2]) ] + seq[ int(query_position[2]) : int(query_position[2] + 15) ]




        # cilp_seq = reads_id + "|" + TE_strand + "|" + genome_strand + "|" + seq[int(clip_left_len) - 50 : int(clip_left_len)] + "|" + seq[int(clip_left_len):int(clip_left_len)+50] 
        # 这里记录的信息是
        # -------------------------------------------------------------------------    genome
        #                  ---------+---------/---------/--------+-------              reads
        #                       1       2          3         4        5       

        cilp_seq =  reads_id  + "|" + \
                    TE_strand  + "|" + \
                    genome_strand + "|" + \
                    seq[ int(query_position[0]) - 50 : int(query_position[0]) ] + "|" + \
                    seq[ int(query_position[0]) : int(query_position[1]) ]  + "|" + \
                    cilp_te   + "|" + \
                    seq[ int(query_position[1]) : int(query_position[2]) ]  + "|" + \
                    seq[ int(query_position[2]) : int(query_position[2] + 50) ] + "|" + \
                    candinate_TSD_left + "|" + \
                    candinate_TSD_right

        if insertion_position[1] - insertion_position[0] >= 20 : 
            insertion_type = "germline"
        else:
            insertion_type = "novel"


        READS_CLUSTER[ref][TE].append([insertion_position[0], insertion_position[1], genome_strand, insertion_type, clip_te_len, cilp_seq, insertion_tag])

    with open(outFile_pickle, "wb") as f1:
        pickle.dump(READS_CLUSTER, f1)
    
    return READS_CLUSTER


