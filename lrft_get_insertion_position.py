
import re

# 修正：需要添加clip left和insertion cigar之间的距离有多长
# query_position 包含里几个特定位置，用以确定具体比对情况
#       ----------------       isnertion     ------------      reference genome
#       ----------------_____insertion seq___------------      map to TE to extract TE 
#       ----------------_____             ___------------      map to genome after clip
#                      |     |           |   |                 query_position
#                   ----_____             ___-----             alignment to find TSD
#                    ---                  ___                  TSD


# 在ciagr的延伸到达clip_left_len的时候，仍然往后判断是否有毛糙的比对
# 如果是一个比较大的deletion的时候，要把insertion修正到这个地方
extend_cigar_len = 20


# 通过结果来看，这个地方好像是更针对于结构比较清晰的 insertion
def get_insertion_position(genome_map_position_start, clip_te_len, clip_left_len, clip_right_len, cigar):
    # 得到的insertion结果包括 genome上insertion的位置
    # 包括query position：reads上clip的位置，以及reads上比对到genome上碱基的位置
    cigar_temp = cigar.replace('H', 'S')
    # cigar_list = re.findall('[0-9]*[A-Z]', cigar_temp.replace('X', 'D'))
    cigar_list = re.findall('[0-9]*[A-Z]', cigar_temp)
    cigar_list_len = len(cigar_list)
    re_type = "novel"
    # 需要先判断比对到基因组上时是否还有clip的字段
    # 因为是从左到右增加，所以只要看左边是否clip就好了
    # print(cigar_list[0])
    if cigar_list[0][-1] == 'S' :
        cigar_len_S = int(cigar_list[0][:-1])
    else:
        cigar_len_S = 0
    cigar_len = 0

    # 因为找insertion是按照从比对起点开始往后找的，所以理应先与clip的左边的reads长度进行比较
    # 如果左边的clip 长度是0，可能是原reads只测到那里
    if clip_left_len == 0 or cigar_len_S >= clip_left_len:
        insertion_position_start = genome_map_position_start - 1
        if insertion_position_start <= 0 :
            insertion_position_start = 1
        insertion_position_end = genome_map_position_start

        insertion_position_to_clip_left = clip_left_len
        insertion_position_to_clip_right = cigar_len_S - clip_left_len
        
        query_position = [ clip_left_len - insertion_position_to_clip_left, clip_left_len , clip_left_len + insertion_position_to_clip_right ]
        re_type = "germline"
        i = 1
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

        for i in range(cigar_list_len):
            cigar_item = cigar_list[i]
            cigar_len = cigar_len + int(cigar_item[:-1])

            cigar_len_dic[cigar_item[-1]] = cigar_len_dic[cigar_item[-1]] + int(cigar_item[:-1])
            cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']
            
            if cigar_len_ex_D < clip_left_len:
                # deletion位置在clip前面的结果
                if cigar_item[-1] == 'D' and int(cigar_item[:-1]) >= 50: 
                    if clip_left_len - cigar_len_ex_D < 200:
                        insertion_position_start = genome_map_position_start + ( cigar_len_dic['M'] +  cigar_len_dic['D'] - int(cigar_item[:-1])  ) -1
                        insertion_position_end = insertion_position_start + int(cigar_item[:-1])
                        
                        insertion_position_to_clip_left = clip_left_len - cigar_len_ex_D
                        insertion_position_to_clip_right = clip_te_len + insertion_position_to_clip_left - int(cigar_item[:-1])
                        for j in range(extend_cigar_len):
                            if i+j+1 < cigar_list_len-1:
                                cigar_item_next = cigar_list[i+j+1]
                                cigar_len = cigar_len + int(cigar_item_next[:-1])
                                cigar_len_dic[cigar_item_next[-1]] = cigar_len_dic[cigar_item_next[-1]] + int(cigar_item_next[:-1])
                                cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']

                                if abs(cigar_len_ex_D - clip_left_len) <= 40 and  cigar_item_next[-1] == 'D'  and int(cigar_item_next[:-1]) >= 50: 
                                    insertion_position_end = genome_map_position_start + ( cigar_len_dic['D'] +  cigar_len_dic['M'] ) -1                                
                                    insertion_position_to_clip_left = 0
                                    insertion_position_to_clip_right = cigar_len_ex_D - clip_left_len

                                    query_position = [ clip_left_len - insertion_position_to_clip_left,  clip_left_len , clip_left_len  + insertion_position_to_clip_right ]
                                    break
                        if insertion_position_to_clip_right > 0:
                            query_position = [ cigar_len_ex_D, clip_left_len, clip_left_len ]
                        else:
                            query_position = [ cigar_len_ex_D, clip_left_len, clip_left_len - insertion_position_to_clip_right ]
                        break
                
                if i < cigar_list_len - 1:
                    cigar_item_next = cigar_list[i+1]
                    if cigar_item_next[-1]=='H' or cigar_item_next[-1]=='S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        insertion_position_end = insertion_position_start + 1
                        
                        insertion_position_to_clip_left = clip_left_len - cigar_len_ex_D
                        insertion_position_to_clip_right = 0

                        query_position = cigar_len_ex_D
                        query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len + insertion_position_to_clip_right ]

                        re_type="germline"
                        break
                else:
                    print('here')
            else:
            # if cigar_len_ex_D >= clip_left_len:
                # 当延伸超过左边的长度时，需要判断下一个cigar是什么来判断insertion end在哪
                # 其实就是前面这段match的长度加deletions的长度
                if cigar_item[-1] == 'M':
                    insertion_position_start = genome_map_position_start + ( clip_left_len - cigar_len_dic['S'] + cigar_len_dic['D'] - cigar_len_dic['I'] ) -1
                    insertion_position_end = insertion_position_start + 1
                    
                    insertion_position_to_clip_left = 0
                    insertion_position_to_clip_right = 0

                    query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len,  clip_left_len + insertion_position_to_clip_right ]
                    
                    # print('ck4')
                    # print(query_position)
                    # deletion在clip位置后面的情况
                    # 在insertion位置两边延伸一个大概TSD的长度
                    # 
                    for j in range(extend_cigar_len):
                        if i+j+1 < cigar_list_len-1:
                            cigar_item_next = cigar_list[i+j+1]
                            cigar_len = cigar_len + int(cigar_item_next[:-1])
                            cigar_len_dic[cigar_item_next[-1]] = cigar_len_dic[cigar_item_next[-1]] + int(cigar_item_next[:-1])
                            cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']

                            if abs(cigar_len_ex_D - clip_left_len) <= 40 and  cigar_item_next[-1] == 'D'  and int(cigar_item_next[:-1]) >= 50: 
                                insertion_position_end = genome_map_position_start + ( cigar_len_dic['D'] +  cigar_len_dic['M'] ) -1                                
                                insertion_position_to_clip_left = 0
                                insertion_position_to_clip_right = cigar_len_ex_D - clip_left_len

                                query_position = [ clip_left_len - insertion_position_to_clip_left,  clip_left_len , clip_left_len  + insertion_position_to_clip_right ]
                                break

                    

                elif cigar_item[-1] == 'I':
                    print('ckkkk')
                    insertion_position_start = genome_map_position_start + ( cigar_len_dic['M'] + cigar_len_dic['D'] ) -1
                    insertion_position_end = insertion_position_start + 1

                    insertion_position_to_clip_left = clip_left_len - ( cigar_len_ex_D - int(cigar_item[:-1]) ) 
                    insertion_position_to_clip_right = int(cigar_item[:-1]) - insertion_position_to_clip_left
                                       
                    query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len + insertion_position_to_clip_right ]
                    for j in range(extend_cigar_len):
                        if i+j+1 < cigar_list_len-1:
                            cigar_item_next = cigar_list[i+j+1]
                            cigar_len = cigar_len + int(cigar_item_next[:-1])
                            cigar_len_dic[cigar_item_next[-1]] = cigar_len_dic[cigar_item_next[-1]] + int(cigar_item_next[:-1])
                            cigar_len_ex_D = cigar_len_dic['S'] + cigar_len_dic['M'] + cigar_len_dic['I']

                            if abs(cigar_len_ex_D - clip_left_len) <= 40 and  cigar_item_next[-1] == 'D'  and int(cigar_item_next[:-1]) >= 50: 
                                insertion_position_end = genome_map_position_start + ( cigar_len_dic['D'] +  cigar_len_dic['M'] ) -1                                
                                insertion_position_to_clip_left = 0
                                insertion_position_to_clip_right = cigar_len_ex_D - clip_left_len

                                query_position = [ clip_left_len - insertion_position_to_clip_left,  clip_left_len , clip_left_len  + insertion_position_to_clip_right ]
                                break

                # elif s[-1] == 'D':
                # 不存在D的情况，因为碱基的累积不会停在缺失的地方
                    
                else:
                    if cigar_item[-1] == 'H' or cigar_item[-1] == 'S':
                        insertion_position_start = genome_map_position_start + cigar_len_dic['M'] + cigar_len_dic['D'] - 1
                        insertion_position_end = insertion_position_start + 1

                        insertion_position_to_clip_left =  clip_left_len - ( cigar_len_ex_D - int(cigar_item[:-1]) ) 
                        insertion_position_to_clip_right = 0
                        
                        query_position = [clip_left_len - insertion_position_to_clip_left, clip_left_len, clip_left_len  + insertion_position_to_clip_right ]
                break
                
    
    return [insertion_position_start, insertion_position_end, ''.join(cigar_list[i-1:i+2]), query_position, re_type ]

