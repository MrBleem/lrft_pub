
from operator import itemgetter
from collections import Counter
import pickle
import re

score_matrix = {'-':-2, 'N':2 }

# m1
# 以15为窗口，遍历这个tsd，找到分数最高的那个
# window_size = 15


# m2
# 找到一个只有一个gap的字符串
# 没有长度限制


# tsd_seq = '---AAC-GATGGG----T---'
# tsd_seq = 'A-A--ACA--AACAGATG--'
# tsd_seq = '-------A--ACA--AACAGATG-GGGCATC--------'

# tsd_seq = '-------A--ACA--AACAGATG-GGGCATC--------'
# tsd_seq = '---AAC-GATGGG----T---'
# tsd_genome = 'AAAAAACAAAAACAGATGGG'
# tsd_genome = 'AAAAAACAAAAACAGATGGG'



def get_candinate_tsd_with_seq(tsd):
    candinate_tsd = []
    tsd_split = tsd.split('-')
    if len(tsd_split) == 1:
        candinate_tsd.append([0, tsd, 0])
        return candinate_tsd[0]
    for i in range(len(tsd.split('-'))):
        if tsd_split[i] != '':
            tsd_tmp = ''
            tsd_tmp_score = 0
            gap_flag = 0

            start_index = len( '-'.join(tsd_split[0:i]) ) + 1

            for j in range(start_index, len(tsd)):
                if tsd[j] != '-':
                    tsd_tmp = tsd_tmp + tsd[j]
                    tsd_tmp_score = tsd_tmp_score + 2
                    if j == len(tsd) and len(tsd_tmp) >= 4:
                        candinate_tsd.append([i, tsd_tmp, tsd_tmp_score])
                else:
                    gap_flag = gap_flag + 1 
                    if gap_flag <= 3:
                        if len(tsd_tmp) >= 4:
                            if len( re.findall('T{4,}', tsd_tmp) ) or len( re.findall('A{4,}', tsd_tmp) )>= 1 :
                                candinate_tsd.append([i, tsd_tmp, tsd_tmp_score-10])
                            else:
                                candinate_tsd.append([i, tsd_tmp, tsd_tmp_score])
                        tsd_tmp = tsd_tmp + tsd[j]
                        tsd_tmp_score = tsd_tmp_score - 4
                    else:
                        break
    if len(candinate_tsd) <= 0:
        return ['NA', 'NA', 'NA']
    TSD_sorted = sorted( candinate_tsd, key = itemgetter(2), reverse=1 )
    # print(TSD_sorted)
    TSD = TSD_sorted[0]

    return TSD


def score_tsd_genome(seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score = score + 2
        else:
            score = score - 2
    return score

def get_tsd_in_genome(tsd_seq, tsd_genome):
    # 把得到的candidate tsd与genome的序列进行比较，看tsd在genome上的位置
    TSD_with_seq = get_candinate_tsd_with_seq(tsd_seq)
    # print('>>>>')
    # print(TSD_with_seq)
    # TSD_with_seq = [20,'CAGAC-CG-TATTT']
    # print(tsd_genome)
    candinate_tsd_with_genome = []
    for i in range(len(tsd_genome)):
        window_genome = tsd_genome[i:i+len(TSD_with_seq[1])]
        candinate_tsd_with_genome.append([i, window_genome, score_tsd_genome(window_genome, TSD_with_seq[1])])
    # print(sorted( candinate_tsd_with_genome, key = itemgetter(2) ))
    tsd_in_genome = sorted( candinate_tsd_with_genome, key = itemgetter(2) )[-1]
    # print(tsd_in_genome)
    # print('<<<<<<')
    return tsd_in_genome





