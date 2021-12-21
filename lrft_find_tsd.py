

from lrft_consensus_tsd import consensus_seq

# 获得的序列是break points(左右两个点)
# 如果用consensus sequence的话，不太能将TSD放入到insertion position的修正当中
# 因为不同的supporting reads TSD比对的情况都不相同，不好修正 
# 所以需要逐个考虑


# 目前选取clip位置左右500bp的reads
# 想起来需要知道clip位置与break point之间差了几个bp
# 因为在找clip left的长度的时候，cigar的延伸不一定刚好到clip left的长度



supporting_reads = { '9396933': 'ATGGCCTCACAAATCTTAACTCGCCACTTCAAGTTAAAACCATCAGAGATGGCTGGGCATGGTGCTCACACCTGTAATCCCAGCACTTTGGGAGGCCAAGGCGTGGCGGATCACGTAGGTCAGAGATCGCACCATCCTGGCTGGGTCACAGTGAAACCCTGTCTCTACTAAAATACAAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAGTCCCAGTAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCAGACCTGGTGAGTGAGCCAAGATCACACCACTGCACTCCAGCCTGATGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAACCATCAGAAGTTCCACTATTGTAAATTT'
                    }
ALU_seq = 'GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTC'
mafft_tmp = open('./tmp.fa', 'w')



for read in supporting_reads:
    read_seq = supporting_reads[read]
    read_seq_left = read_seq[:50]
    read_seq_right = read_seq[-50:]
    mafft_tmp.write('>left\n')
    mafft_tmp.write(read_seq_left+'\n')
    # mafft_tmp.write(left_test + '\n')
    mafft_tmp.write('>right\n')
    mafft_tmp.write(read_seq_right)
    # mafft_tmp.write(right_test)
    mafft_tmp.close()
    fre_cuoff = 1

    substitude_tsd = consensus_seq('./tmp.fa', fre_cuoff)
    print(substitude_tsd)




