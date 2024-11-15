import os
import sys
import re
from tqdm import tqdm  # 引入 tqdm 以添加进度条

def uncouple_reads(fhi, barcodes, fho, stats):
    '''
    输入一个排序的BED文件处理的文件接口，和一个名称和标识码的关联列表，並输出匹配结果和统计信息
    :param fhi: 输入 BED 文件的文件接口
    :param barcodes: 标识码列表的文件接口
    :param fho: 输出的结果文件接口
    :param stats: 输出匹配的统计结果
    '''
    matched = 0
    unmatched = 0
    linkage = {}  # 连接关联名称和标识码状态的字典
    unlinked = {}  # 进行匹配失败的字典

    # 处理标识码文件
    for line in barcodes:
        name, bc1, bc2, bcterm = line.split()
        name_correct = name.split("/")[0]
        if bc1 == bc2:  # 检查标识码的两半是否匹配
            barcode = bc1  # 如果匹配则记录标识码
            matched += 1
            linkage[name_correct] = [barcode, bcterm]
        else:
            barcode = "%s-%s" % (bc1, bc2)
            unlinked[name_correct] = [barcode, bcterm]
            unmatched += 1

    stats.write(f"<H3>Read Pairs with Matched Barcodes: {matched}</H3>\n<H3>Read Pairs with Unmatched Barcodes: {unmatched}</H3>\n")
    print(f"Total matched barcodes: {matched}, Total unmatched barcodes: {unmatched}")

    # 读取输入文件的总行数，用于显示进度条
    total_lines = sum(1 for _ in fhi)
    fhi.seek(0)  # 将文件指针重置回开头

    read_lb = [1, 2, 3, 4, 5, 6]  # 初始化值，用于记录后续进程

    # 使用 tqdm 添加进度条
    for line in tqdm(fhi, total=total_lines, desc="Processing BED records", mininterval=10):
        chrom, fend, rend, name, mapq, strand, chrom_frag, fend_frag, rend_frag, dist = line.split()
        name_correct = name.split("/")[0]
    
        # 判断该条记录是否在标识码字典中（匹配成功或失败）
        if name_correct in linkage:
            # 如果当前记录与之前一条记录的名称相同
            if read_lb[3] == name_correct:
                # 如果染色体相同，进入匹配逻辑
                if chrom == read_lb[0]:
                    # 根据 fend 和 rend 的条件选择不同的写入逻辑
                    # 情况 1: 当前记录的 fend 或 rend 与上一条记录的 fend 或 rend 相等，说明这两条记录的位置在染色体上非常接近，可能是互补配对的情况，在这种情况下，将当前记录和前一条记录拼接起来，标记配对(+/-)
                    if int(fend) == int(read_lb[1]) or int(rend) == int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t+\t-\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{dist}\n")
                    elif int(fend) < int(read_lb[2]):
                    # 情况 2: 当前记录的 fend 小于上一条记录的 rend，表明当前记录在前一条记录的区域内，因此它们之间存在一定的重叠，这种情况下，标记配对的方向，并保留原始的方向信息
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{strand}\t{read_lb[5]}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{dist}\n")

                    else: # 情况 3: 当前记录的 fend 大于前一条记录的 rend
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{dist}\n")
                else:
                    # 如果染色体不同
                    if chrom < read_lb[0]:
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{strand}\t{read_lb[5]}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{dist}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{dist}\n")

        elif name_correct in unlinked:
            # 处理没有匹配成功的情况
            if read_lb[3] == name_correct:
                if chrom == read_lb[0]:
                    if int(fend) == int(read_lb[1]) or int(rend) == int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t+\t-\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{dist}\n")
                    elif int(fend) < int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{strand}\t{read_lb[5]}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{dist}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{dist}\n")
                else:
                    if chrom < read_lb[0]:
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{strand}\t{read_lb[5]}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{dist}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{dist}\n")

        # 更新 read_lb 为当前行的值
        read_lb = [chrom, fend, rend, name_correct, mapq, strand, chrom_frag, fend_frag, rend_frag, dist]

def main():
    fhi = open(sys.argv[1])  # 输入 BED 文件
    barcodes = open(sys.argv[2])  # 输入标识码文件
    mismatch_out = open(sys.argv[3], 'w')  # 输出不匹配的文件
    outfile = open(sys.argv[4], 'w')  # 输出的文件
    uncouple_reads(fhi, barcodes, mismatch_out, outfile)
    fhi.close()
    barcodes.close()
    mismatch_out.close()
    outfile.close()

if __name__ == "__main__":
    main()


