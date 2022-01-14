# %%

from operator import itemgetter

score_matrix = {'-':-2, }
window_size = 10

# m1
# 以10为窗口，遍历这个tsd，找到分数最高的那个

# m2
# 找到一个只有一个gap的字符串
# 没有长度限制

tsd_seq = '---AAC-GATGGG----T---'
tsd_seq = 'A-A--ACA--AACAGATG--'
tsd_seq = '-------A--ACA--AACAGATG-GGGCATC--------'

'ACACCACAGCAACAGATGCA'
'AAAAAACAAAAACAGATGGG'
'A-A--ACA--AACAGATG--'

# tsd_left =   'AACATCAAGTAAATGAAAGAGATAATGCCA'
# tsd_right =  'TGCCAAAAAAAAAAAAAAAAAGAGATAATG'
# tsd_genome = 'AACATCAAGTAAATGAAAAAAGAGATAATG'
#            '-----CAA--AAA-------AAAGAGATAATG--- '



####
# methods2
####

tsd_seq = '-------A--ACA--AACAGATG-GGGCATC--------'
tsd_seq = '---AAC-GATGGG----T---'
tsd_genome = 'AAAAAACAAAAACAGATGGG'


def get_candinate_tsd_with_seq(tsd):
    candinate_tsd = []
    tsd_split = tsd.split('-')
    for i in range(len(tsd.split('-'))):
        if tsd_split[i] != '':
            tsd_tmp = ''
            tsd_tmp_score = 0
            gap_flag = 0

            start_index = len( '-'.join(tsd_split[0:i]) ) + 1
            
            for j in range(start_index, len(tsd)):
                # print(tsd[j])
                if tsd[j] != '-':
                    tsd_tmp = tsd_tmp + tsd[j]
                    tsd_tmp_score = tsd_tmp_score + 2
                else:
                    gap_flag = gap_flag + 1 
                    if gap_flag < 2:
                        tsd_tmp = tsd_tmp + tsd[j]
                        tsd_tmp_score = tsd_tmp_score + 0
                    else:
                        break
            candinate_tsd.append([i, tsd_tmp, tsd_tmp_score])
    TSD = sorted( candinate_tsd, key = itemgetter(2) )[-1]
    if len(TSD[1]) < 5:
        TSD = 'NA'
    return TSD

# 把得到的candidate tsd与genome的序列进行比较，看tsd在genome上的位置
tsd_genome = 'AAAAAACAAAAACAGATGGG'


def score_tsd_genome(seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score = score + 2
        else:
            score = score - 2
    return score

def get_tsd_in_genome(tsd_seq, tsd_genome):
    TSD_with_seq = get_candinate_tsd_with_seq(tsd_seq)
    candinate_tsd_with_genome = []
    for i in range(len(tsd_genome)):
        window_genome = tsd_genome[i:i+len(TSD_with_seq[1])]
        candinate_tsd_with_genome.append([i, window_genome, score_tsd_genome(window_genome, TSD_with_seq[1])])

    tsd_in_genome = sorted( candinate_tsd_with_genome, key = itemgetter(2) )[-1]
    return tsd_in_genome




print( get_tsd_in_genome(tsd_seq, tsd_genome) )
# 这样子可以得到TSD在consensus seq上的位置
# 但是并不知道分别在两端的位置具体是怎么样的

# %%
####
# methods1
####

# candinate_tsd = []
# for i in range(len(tsd) - window_size):
#     if tsd[i] != "-":
#         tsd_window = tsd[i: i + window_size]
#         tsd_window_score = 0
#         for j in range(window_size):
#             if tsd_window[j] == "-":
#                 tsd_window_score = tsd_window_score - 2
#             else:
#                 tsd_window_score = tsd_window_score + 2
#         candinate_tsd.append( [ i, tsd_window, tsd_window_score ] )


