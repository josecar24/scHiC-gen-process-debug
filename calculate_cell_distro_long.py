import os
import sys
import re
from collections import Counter
import random

# 打开输入文件
fhi = open(sys.argv[1])  # Open filehandle of associations bedpe
cnt = Counter()  # Set up a Counter object for counting the occurrence of barcode pairs
frag_occurrence = {}
bc = {}  # Dictionary for barcode associations
bc_long = {}  # Barcode associations for those above a certain distance class
bc_random = {}  # Dictionary for random barcode associations
bc_random_long = {}  # Dictionary for long-distance random associations
bcs = []  # Store all the barcode labels to shuffle them later

# 第一遍遍历，统计条形码对的出现次数，并生成条形码组合列表
for line in fhi:
    split = line.split()
    cnt["%s-%s" % (split[11], split[12])] += 1
    bcs.append("%s-%s" % (split[11], split[12]))

# 打乱条形码组合，供随机使用
random.shuffle(bcs)

# 初始化条形码关联的字典
for combo in cnt:
    frag_occurrence[combo] = Counter()
    bc[combo] = {"human": 0}
    bc_random[combo] = {"human": 0}
    bc_long[combo] = {"human": 0}
    bc_random_long[combo] = {"human": 0}
fhi.close()

# 重新打开文件以进行第二遍遍历
fhi = open(sys.argv[1])
i = 0  # Counter for random barcodes
for line in fhi:
    split = line.split()
    barcode = "%s-%s" % (split[11], split[12])
    fragid1 = split[13]
    fragid2 = split[15]
    strand1, strand2 = split[9:11]
    dist1, dist2 = split[14], split[16]

    # 根据片段信息更新frag_occurrence字典
    if fragid1 == fragid2:
        frag_occurrence[barcode][(fragid1, strand1, dist1)] += 1
    elif split[1] == split[4]:
        frag_occurrence[barcode][(fragid1, strand1, dist1)] += 1
    else:
        frag_occurrence[barcode][(fragid1, strand1, dist1)] += 1
        frag_occurrence[barcode][(fragid2, strand2, dist2)] += 1

    # 更新随机条形码的统计
    random_bc = bcs[i] # 从打乱后的条形码列表中取一个
    if re.search("human", split[0]):
        if re.search("human", split[3]):
            # 同属于人类物种且为同一条染色体
            bc[barcode]["human"] += 1
            bc_random[random_bc]["human"] += 1
            # 判断是否为长距离片段（大于1000个碱基）
            if int(split[4]) - int(split[1]) > 1000:
                bc_long[barcode]["human"] += 1
                bc_random_long[random_bc]["human"] += 1
        else:
            # 不同染色体或不同片段但仍然属于人类物种
            bc[barcode]["human"] += 1
            bc_random[random_bc]["human"] += 1
            bc_long[barcode]["human"] += 1
            bc_random_long[random_bc]["human"] += 1
    i += 1
fhi.close()

# 计算并输出统计信息
for combo in cnt:
    total = bc[combo]["human"]  #对于contact次数列的统计
    total_random = bc_random[combo]["human"]
    if total == 0:
        continue
    # 计算人类物种的片段比例
    frac_human = float(bc[combo]["human"]) / float(total)
    # 如果total_random为0，随机比例为0，否则计算随机的比例
    frac_human_random = 0 if total_random == 0 else float(bc_random[combo]["human"]) / float(total_random)

    # 将barcode对拆分为bc1和bc2
    bc1, bc2 = combo.split("-")
    ones, twos, threes, fours = 0, 0, 0, 0
    # 遍历每个片段信息，统计片段数量
    for fragid in frag_occurrence[combo]:
        if frag_occurrence[combo][fragid] == 1:
            ones += 1
        elif frag_occurrence[combo][fragid] == 2:
            twos += 1
            sys.stderr.write(f"{fragid[0]}\t{fragid[1]}\t{fragid[2]}\t{frag_occurrence[combo][fragid]}\t{combo}\n")
        elif frag_occurrence[combo][fragid] == 3:
            threes += 1
            sys.stderr.write(f"{fragid[0]}\t{fragid[1]}\t{fragid[2]}\t{frag_occurrence[combo][fragid]}\t{combo}\n")
        elif frag_occurrence[combo][fragid] >= 4:
            fours += 1
            sys.stderr.write(f"{fragid[0]}\t{fragid[1]}\t{fragid[2]}\t{frag_occurrence[combo][fragid]}\t{combo}\n")

    # 打印出真实和随机条形码对的统计信息（所有长度）
    print(f"{frac_human}\t{frac_mouse}\t{bc[combo]['human']}\t{bc[combo]['mouse']}\t{total}\t{cnt[combo]}\t{bc1}\t{bc2}\tTrue\tAll\t{ones}\t{twos}\t{threes}\t{fours}")
    print(f"{frac_human_random}\t{frac_mouse_random}\t{bc_random[combo]['human']}\t{bc_random[combo]['mouse']}\t{total_random}\t{cnt[combo]}\t{bc1}\t{bc2}\tRandomized\tAll\t{ones}\t{twos}\t{threes}\t{fours}")

    # 再计算长距离片段的统计信息
    total = bc_long[combo]["human"] 
    if total == 0:
        continue
    frac_human = float(bc_long[combo]["human"]) / float(total)
    total_random = bc_random_long[combo]["human"]
    frac_human_random = 0 if total_random == 0 else float(bc_random_long[combo]["human"]) / float(total_random)

    # 打印出长距离片段的统计信息
    print(f"{frac_human}\t{frac_mouse}\t{bc_long[combo]['human']}\t{bc_long[combo]['mouse']}\t{total}\t{cnt[combo]}\t{bc1}\t{bc2}\tTrue\tLong\t{ones}\t{twos}\t{threes}\t{fours}")
    print(f"{frac_human_random}\t{frac_mouse_random}\t{bc_random_long[combo]['human']}\t{bc_random_long[combo]['mouse']}\t{total_random}\t{cnt[combo]}\t{bc1}\t{bc2}\tRandomized\tLong\t{ones}\t{twos}\t{threes}\t{fours}")
