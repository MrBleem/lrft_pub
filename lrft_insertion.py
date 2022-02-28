
####

# 比v1添加了insertion类型的判断，暂时将insertion的长度作为判断标准 

# 比v2添加了trunked的判断，就是截断reads的insertion的判断

# 比v3:尝试全部模块化

# 比v4:添加了TSD的模块

# 比v5:Tcandadate TSD筛选多添加了一些

# 比v6:需要再完成的事情：
        # 只有一端的clip reads
        # 按染色体位置进行遍历
        # germlien 和 novel insertion的区别

# 想办法把complement的insertion找出来，这可能可以是一个亮点吧



####



import datetime
starttime = datetime.datetime.now()
#long running

import os
import sys
import pickle
import subprocess
import pysam
from collections import Counter

from lrft_find_tsd import get_tsd_in_genome
from lrft_consensus import get_consensus_seq
from lrft_cluster import Cluster_reads_by_TE

# help文档
# if __name__ == '__main__':
#     helpdoc = '''
# Descrtption:
#     long reads for transposon

# Usage
#     python lrft_insertion -i path_to_bam_file -g path_to_genome_fasta -o out_path -p prefix

# Parameter
#     -h/--help
#         Print helpdoc
#     -i/--infile
#         Input file
#     -g/--genomeFa
#         Genome fasta
#     -o/--outpath
#         Out path
#     -p/--prefix
#         Result prefix
#     '''
#     try:
#         opts,args = getopt.getopt(sys.argv[1:], "hi:g:o:p:", ["infile=", "genomeFa=", "outpath=", "prefix="])
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
#         elif opt in ("-g","--genomeFa"):
#             GENOME_FA = arg
#         elif opt in ("-o","--outpath"):
#             OUTPATH = arg
#         elif opt in ("-p", "--prefix"):
#             PREFIX = arg



# 这个版本只用一种insertion的方法， max_deletion的方法暂时不用
# 因为想在大部分物种中使用，主要依靠一种思想比较靠谱

# imput
PROJECT_PATH = sys.argv[1]     # "/data/tusers/boxu/lrft/result/nanopore/human"
BAMFILE = sys.argv[2]         # "SRR11669560.mapped.sorted.tcg.bam"
GENOME_FA = sys.argv[3]      # "/data/tusers/boxu/annotation/hs37d5/hs37d5.fa"
PREFIX = sys.argv[4]        # "SRR11669560"


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
 
# region="2:40,613,659-40,613,752"
# region = "1:23,981,124-23,981,132" # AACAGATGGGGCATC
# region = '1:77,764,983-77,764,992'
# region = '1:112,991,995-112,992,009' # CAGACACGTTATTT
# region = '1:234,805,709-234,805,732' # TGATTATATGTTTT
# region = '1:246,592,328-246,592,345' # AGTCTCACTTT
# region = '10:9,474,847-9,474,855' # TTTTT
# region = '10:26,852,961-26,852,978' # AAATAAAAAGAAAA
# region = '7:107,829,249-107,829,298'
# region = '1:41,593,610-41,602,609'
# region= '1:12,790,468-12,791,007'

# chr = region.split(':')[0]
# region_start = region.split(':')[1].split('-')[0].replace(',', '')
# region_end = region.split(':')[1].split('-')[1].replace(',', '')
# READS = bamfile.fetch( chr, int(region_start) , int(region_end))

outFile = RESULT_PATH + PREFIX + ".lrft.table.txt"
outfile = open(outFile, 'w')
outfile.write('chr\tinsertion_start\tinsertion_end\tTE\tstrand\tsurpporting_reads\tun_surpporting_reads\tfrequency\tinsertion_type\tisnertion_te_length\tinsertion_te_sequence\n')


outFile_pickle = RESULT_PATH + PREFIX + ".lrft.cluster.pkl"
SEQ_INFO_FILE = DATA_PATH + PREFIX.split('.')[0] + ".clip.reads.pkl"
# GENOME_FA = "/data/tusers/boxu/annotation/hs37d5/hs37d5.fa"
genome_fa_seq = pysam.Fastafile(GENOME_FA)


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
# 取并集的时候还是需要看一下方向的问题
def merge_insertion(insertion1, insertion2):
    # merge insertion应该只去merge方向相同的，所以方向判断应该在前一步

    # --------      insertion1
    #    --------   insertion2
    # -----------   merged_insertion
    # ++++++++
    #    ++++++++
    # 
    new_insertion_start = insertion1[0]
    new_insertion_end = insertion1[1]
    # if insertion2[0] < new_insertion_s:
    #     new_insertion_s = insertion2[0]
    if insertion2[1] > new_insertion_end:
        new_insertion_end = insertion2[1] 
    
    if insertion1[6] == "truncked" and insertion2[6] == "truncked":
        insertion_tag = "truncked"
    else:
        insertion_tag = "spanning"
        
    return [new_insertion_start, new_insertion_end, insertion1[2],  insertion1[3] + ";" + insertion2[3], sum([insertion1[4],insertion2[4]])/2,  insertion1[5]+";"+insertion2[5], insertion_tag ]


def record_insertion(ref, te, insertion):
    # 将insertion记录到文件中
    # 这一步包括：
    #       · 判断insertion的方向 
    #       · insertion的consensus sequence
    
    insertion_start = insertion[0]
    insertion_end = insertion[1]

    supporting_reads_num = insertion[2]
    spanding_reads_num = insertion[3]
    if spanding_reads_num == 0 :
        spanding_reads_num = supporting_reads_num

    insertion_strand = []
    insertion_strand_te = []
    insertion_strand_genome = []

    insertion_type = insertion[4]
    mapped_te_position = insertion[5]
    READS_SEQ_INFO = insertion[6]
    READS_SEQ = []

    insertion_tag = insertion[7]

    

    # align_seq_file_name = RESULT_PATH + PREFIX + ".seq.temp." + str(test_index) + ".fq"
    # align_tsd_file_name = RESULT_PATH + PREFIX +".tsd.temp." + str(test_index) + ".fq"

    align_seq_file_name = RESULT_PATH + PREFIX + ".seq.temp.fq"

    align_tsd_file_name = RESULT_PATH + PREFIX + ".tsd.temp.fq"

    align_tsd_left_file_name = RESULT_PATH + PREFIX + ".tsd.temp.1.fq"
    align_tsd_right_file_name = RESULT_PATH + PREFIX + ".tsd.temp.2.fq"

    align_seq_file = open( align_seq_file_name, 'w' )
    align_tsd_file = open( align_tsd_file_name, 'w' )

    align_tsd_left_file = open( align_tsd_left_file_name, 'w' )
    align_tsd_right_file = open( align_tsd_right_file_name, 'w' )


    ref_seq = genome_fa_seq.fetch(ref, int(insertion_start) - 20 , int(insertion_end) + 19 )

    if int(insertion_end) - int(insertion_start) <= 20:
        tsd_genome = ref_seq
        print(tsd_genome)
        print(insertion_start, insertion_end)
    else:
        tsd_genome = "".join([ref_seq[0:20], ref_seq[-20:]])

    left_tsd_num = 0
    right_tsd_num = 0

    for read_seq in READS_SEQ_INFO.split(';'):           
        read_seq_info = read_seq.split('|')
                
        read_id = read_seq_info[0]
        te_strand = read_seq_info[1]
        genome_strand = read_seq_info[2]

        insertion_sequence = read_seq_info[3:]

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
        
        # 以reference 正向为导向
        if genome_strand == "+":
            insertion_strand.append("+")
        else:
            insertion_strand.append("-")
            insertion_sequence[2] = comp_reverse(insertion_sequence[2])
        
        READS_SEQ.append(read_id + "|" + \
                        insertion_sequence[0] + \
                        "\033[33m" + insertion_sequence[1] + "\033[0m" + \
                        "\033[34m" + insertion_sequence[2] + "\033[0m" + \
                        "\033[33m" + insertion_sequence[3] + "\033[0m" + \
                        insertion_sequence[4] )

        # tsd
        # align_tsd_file.write( '>' + read_id + ".left" + '\n' )
        # align_tsd_file.write( insertion_sequence[5] + '\n' )

        # align_tsd_file.write( '>' + read_id + ".right" + '\n' )
        # align_tsd_file.write( insertion_sequence[6] + '\n' )

        # tsd left 和 right 分别先做consensus
        align_tsd_left_file.write( '>' + read_id + ".left" + '\n' )
        align_tsd_left_file.write( insertion_sequence[5] + '\n' )

        align_tsd_right_file.write( '>' + read_id + ".right" + '\n' )
        align_tsd_right_file.write( insertion_sequence[6] + '\n' )

        if len(insertion_sequence[5]) >= 15:
            left_tsd_num += 1
        if len(insertion_sequence[6]) >= 15:
            right_tsd_num += 1

        # consensus seq
        align_seq_file.write( '>' + read_id + '\n' )
        # align_seq_file.write( clip_left + cilp_mid + clip_right + '\n' )
        align_seq_file.write( "".join(insertion_sequence[0:4]) + '\n' )


    # print("TSD_left_right_propotion:" + str(left_tsd_num) +":"+ str(right_tsd_num) )
    align_seq_file.close()
    # align_tsd_file.close()

    align_tsd_right_file.close()
    align_tsd_left_file.close()

    fre_cutoff = 0.75


    if insertion_end - insertion_start >= 50 or insertion_tag == 'truncked': 
        insertion_type = "germline"
        tsd_seq = "NA"
    else:
        insertion_type = "novel"
        tsd_left = get_consensus_seq(align_tsd_left_file_name, fre_cutoff, 'seq')
        tsd_right = get_consensus_seq(align_tsd_right_file_name, fre_cutoff, 'seq')

        align_tsd_file.write( '>left\n' )
        align_tsd_file.write( tsd_left + '\n' )

        align_tsd_file.write( '>right\n' )
        align_tsd_file.write( tsd_right + '\n' )
        align_tsd_file.close()

        tsd_seq  = get_consensus_seq(align_tsd_file_name, fre_cutoff, 'tsd')
    
    # print("type:" + insertion_type + ":TSD_left_right_propotion:" + str(left_tsd_num) +":"+ str(right_tsd_num) )
    tsd_num_ck = str(left_tsd_num) +":"+ str(right_tsd_num) 

    if tsd_seq == 'NA':
        tsd_in_genome = ['NA','NA', 'NA']
    else:
        tsd_in_genome = get_tsd_in_genome(tsd_seq, tsd_genome)

        insertion_start = int(insertion_start)-20 + int(tsd_in_genome[0])
        insertion_end = insertion_start + len(tsd_in_genome[1]) -1 


    # align_reads = get_consensus_seq(align_seq_file_name, fre_cutoff, 'seq')
    align_reads = ''
    align_seq = align_reads + "<<<< >>>>" + tsd_in_genome[1]

    
            
    # 标记检测，测试用的
    if len(Counter(insertion_strand).items()) == 1 and list(Counter(insertion_strand).items())[0][1] != 1 :
        insertion_strand.append('CAO')
    elif len(Counter(insertion_strand).items()) == 2:
        insertion_strand.append('GAN')
                

    strand_info = ['_'.join(insertion_strand_te), '_'.join(insertion_strand_genome), '_'.join(insertion_strand)] 
    insertion_final = [str(ref), str(insertion_start), str(insertion_end), te, '|'.join(strand_info), str(supporting_reads_num) + "|" + tsd_num_ck, str(spanding_reads_num-supporting_reads_num), str(float(supporting_reads_num)/spanding_reads_num), insertion_type, str(mapped_te_position), "\033[31m" + tsd_in_genome[1] + "\033[0m", align_seq, '_'.join(READS_SEQ)]
    
    
    outfile.write('\t'.join(insertion_final) + '\n')

with open(SEQ_INFO_FILE, "rb") as f3:
    READS_SEQUENCE = pickle.load(f3)

# 加载前一步的文件
if os.path.exists(outFile_pickle):
    with open(outFile_pickle, "rb") as f2:
        READS_CLUSTER = pickle.load(f2)
    print('get reads cluster from pkl file')
else:
    READS_CLUSTER = Cluster_reads_by_TE(READS, READS_SEQUENCE, outFile_pickle)




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
        # cluster insertion by position
        for j in range(len(TE_CLUSTER[te][1:])):
            insertion = TE_CLUSTER[te][1:][j]
            # 以前一个insertion为参考，判断后一个insertion是否跟前一个有overlap
            # 这里判断的标准是看后一个insertion的中点是否在前一个insertion当中
            # 如果两个insertion有overlap，那么这两个insertion就合并在一起，区间合并在一起

            # if sum(insertion[0:2])/2 > insertion_pre[0] and sum(insertion[0:2])/2 < insertion_pre[1]:
            # insertion_pre[1] = insertion_pre[1] + 20    # 加入序列信息之后需要去掉 flag=delete

            # if insertion[0] >= insertion_pre[0] and insertion[0] <= insertion_pre[1]+20 and insertion[2] == insertion_pre[2]: ????? insertion[2] == insertion_pre[2]这个是为了判断insertion type 用的
            if insertion[0] >= insertion_pre[0] and insertion[0] <= insertion_pre[1]+20 and insertion[2] == insertion_pre[2] : # 判断方向
                # print(insertion[2],insertion_pre[2])
                insertion_pre = merge_insertion(insertion_pre, insertion) 
                supporting_reads = supporting_reads + 1
            else:
                #   mapped_reads是根据insertion的位置，再用pysam计算insertion的位置上有多少条reads比对到了这个位置
                mapped_reads = bamfile.count(contig=ref, start=insertion_pre[0], stop=insertion_pre[1]+500)
                # insertion_start, insertion_end, insertion_type, clip_te_len, cilp_seq, insertion_tag
                insertion_record = [insertion_pre[0], insertion_pre[1], supporting_reads, mapped_reads,  insertion_pre[3], insertion_pre[4], insertion_pre[5], insertion_pre[6]  ]
                record_insertion(ref, te, insertion_record)
                insertion_pre = insertion
                supporting_reads = 1
                if j == len(TE_CLUSTER[te][1:]):
                    flag='end'
                    
        if flag !='end' :
            mapped_reads = bamfile.count(contig=ref, start=insertion_pre[0], stop=insertion_pre[1]+500)
            # insertion_start, insertion_end, insertion_type, clip_te_len, cilp_seq, insertion_tag
            insertion_record = [insertion_pre[0], insertion_pre[1], supporting_reads, mapped_reads,  insertion_pre[3], insertion_pre[4], insertion_pre[5], insertion_pre[6] ]
            record_insertion(ref, te, insertion_record)

outfile.close()
print('>_< DONE >V<')
endtime = datetime.datetime.now()
print (endtime - starttime)
rm_args = ['rm', outFile_pickle]
subprocess.Popen(rm_args)



