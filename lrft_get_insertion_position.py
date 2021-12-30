
import re




# 修正：需要添加clip left和insertion cigar之间的距离有多长
# query_position 包含里几个特定位置，用以确定具体比对情况
#       ----------------       isnertion     ------------      reference genome
#       ----------------_____insertion seq___------------      map to TE to extract TE 
#       ----------------_____             ___------------      map to genome after clip
#                      |     |           |   |                 query_position
#                   ----_____             ___-----             alignment to find TSD
#                    ---                  ___                  TSD


# 通过结果来看，这个地方好像是更针对于结构比较清晰的 insertion
def get_insertion_position(genome_map_position_start, clip_te_len, clip_left_len, clip_right_len, cigar):
    cigar_temp = cigar.replace('H', 'S')
    # s_cigar = re.findall('[0-9]*[A-Z]', cigar_temp.replace('X', 'D'))
    s_cigar = re.findall('[0-9]*[A-Z]', cigar_temp)
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
        insertion_position_start = genome_map_position_start - 1
        if insertion_position_start <=0 :
            insertion_position_start = 1
        insertion_position_end = genome_map_position_start

        insertion_position_to_clip_left = clip_left_len
        insertion_position_to_clip_right = cigar_len_S - clip_left_len
        
        query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len + clip_te_len, clip_left_len + clip_te_len + insertion_position_to_clip_right ]
        # print('ck0')
        # print(query_position)
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
            
            if cigar_len_ex_D < clip_left_len:
                # deletion位置在clip前面的结果
                if s[-1] == 'D' and int(s[:-1]) >= 50: 
                    if clip_left_len - cigar_len_ex_D < 100:
                        insertion_position_start = genome_map_position_start + ( cigar_len_dic['M'] +  cigar_len_dic['D'] - int(s[:-1])  ) -1
                        insertion_position_end = insertion_position_start + int(s[:-1])
                        
                        insertion_position_to_clip_left = clip_left_len - cigar_len_ex_D
                        insertion_position_to_clip_right = clip_te_len + insertion_position_to_clip_left - int(s[:-1])

                        # print(insertion_position_to_clip_right)

                        # insertion_position_to_clip_right = cigar_len_ex_D + int(s[:-1]) - ( insertion_position_to_clip_left + clip_te_len )

                        if insertion_position_to_clip_right > 0:
                            # query_position = [cigar_len_ex_D, clip_left_len, cigar_len_ex_D + int(s[:-1]), cigar_len_ex_D + int(s[:-1]) + insertion_position_to_clip_right ]
                            # query_position = [ cigar_len_ex_D, clip_left_len, cigar_len_ex_D + int(s[:-1]), cigar_len_ex_D + int(s[:-1]) ]                            
                            query_position = [ cigar_len_ex_D, clip_left_len, clip_left_len, clip_left_len ]
                            # print('ck1')
                            # print(query_position)
                        else:
                            query_position = [ cigar_len_ex_D, clip_left_len, clip_left_len , clip_left_len - insertion_position_to_clip_right ]
                            # print('ck2')
                            # print(query_position)

                        break
                if i < s_cigar_len - 1:
                    s_next = s_cigar[i+1]
                    if s_next[-1]=='H' or s_next[-1]=='S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        # insertion_position_end = insertion_position_start + clip_te_len + 1
                        insertion_position_end = insertion_position_start + 1
                        
                        insertion_position_to_clip_left = clip_left_len - cigar_len_ex_D
                        insertion_position_to_clip_right = 0

                        query_position = cigar_len_ex_D
                        query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len + clip_te_len, clip_left_len + clip_te_len + insertion_position_to_clip_right ]
                        # print('ck3')
                        # print(query_position)

                        re_type="germline"
                        break
            else:
            # if cigar_len_ex_D >= clip_left_len:
                # 当延伸超过左边的长度时，需要判断下一个cigar是什么来判断insertion end在哪
                # 其实就是前面这段match的长度加deletions的长度
                if s[-1] == 'M':
                    insertion_position_start = genome_map_position_start + ( clip_left_len - cigar_len_dic['S'] + cigar_len_dic['D'] - cigar_len_dic['I'] ) -1
                    insertion_position_end = insertion_position_start + 1
                    
                    insertion_position_to_clip_left = 0
                    insertion_position_to_clip_right = 0

                    query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len , clip_left_len + insertion_position_to_clip_right ]
                    
                    # print('ck4')
                    # print(query_position)
                    # deletion在clip位置后面的情况
                    # 在insertion位置两边延伸一个大概TSD的长度
                    # 
                    for j in range(14):
                        if i+j+1 < s_cigar_len-1:
                            s_next = s_cigar[i+j+1]
                            # print(s_next)
                            cigar_len = cigar_len + int(s_next[:-1])
                            cigar_len_dic[s_next[-1]] = cigar_len_dic[s_next[-1]] + int(s_next[:-1])
                            cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']
                            # print( cigar_len_ex_D, clip_left_len )

                            if abs(cigar_len_ex_D - clip_left_len) <= 40 and  s_next[-1] == 'D'  and int(s_next[:-1]) >= 50: 
                                # print('ck')
                                # insertion_position_start = genome_map_position_start + ( clip_left_len - cigar_len_dic['S'] + cigar_len_dic['D'] - cigar_len_dic['I'] ) -1
                                insertion_position_end = genome_map_position_start + ( cigar_len_dic['D'] +  cigar_len_dic['M'] ) -1
                                # insertion_position_start = insertion_position_end - int(s_next[:-1])
                                
                                insertion_position_to_clip_left = 0
                                insertion_position_to_clip_right = cigar_len_ex_D - clip_left_len

                                query_position = [ clip_left_len - insertion_position_to_clip_left, clip_left_len  , clip_left_len , clip_left_len  + insertion_position_to_clip_right ]
                                # print('ck5')
                                # print(query_position)
                                break
                    # if i < s_cigar_len-1:
                    #     s_next = s_cigar[i+1]
                    #     print(s_next)
                    #     if abs(cigar_len_ex_D - clip_left_len) <= 10 and  s_next[-1] == 'D' : 
                    #         print('ck')
                    #         insertion_position_end = insertion_position_start + int(s_next[:-1])
                    

                elif s[-1] == 'I':
                    insertion_position_start = genome_map_position_start + ( cigar_len_dic['M'] + cigar_len_dic['D'] ) -1
                    # insertion_position_end = insertion_position_start + int(s[:-1])
                    insertion_position_end = insertion_position_start + 1

                    insertion_position_to_clip_left = clip_left_len - ( cigar_len_ex_D - int(s[:-1]) ) 
                    insertion_position_to_clip_right = int(s[:-1]) - insertion_position_to_clip_left
                                       
                    query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len , clip_left_len + insertion_position_to_clip_right ]
                    # print('ck6')
                    # print(query_position)
                # elif s[-1] == 'D':
                # 不存在D的情况，因为碱基的累积不会停在缺失的地方
                    
                else:
                    if s[-1] == 'H' or s[-1] == 'S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        # insertion_position_end = insertion_position_start + clip_te_len + 1
                        insertion_position_end = insertion_position_start + 1

                        insertion_position_to_clip_left =  clip_left_len - ( cigar_len_ex_D - int(s[:-1]) ) 
                        insertion_position_to_clip_right = 0
                        
                        query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len , clip_left_len  + insertion_position_to_clip_right ]
                        # print('ck7')
                        # print(query_position)
                break
                
    # return [insertion_position_start, insertion_position_end, ''.join(s_cigar[i-1:i+2]), query_position, re_type]
    return [insertion_position_start, insertion_position_end, ''.join(s_cigar[i-1:i+2]), query_position, re_type ]

