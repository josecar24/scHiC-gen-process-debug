import os, sys, re

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

    for line in barcodes:  # 迭代标识码文件
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

    read_lb = [1, 2, 3, 4, 5, 6]  # 初始化值，用于记录后续进程

    for line in fhi:
        chrom, fend, rend, name, mapq, strand, chrom_frag, fend_frag, rend_frag, frag_id, space, frag_strand, dist = line.split()
        name_correct = name.split("/")[0]
        # 处理匹配成功的情况
        if name_correct in linkage:
            if read_lb[3] == name_correct:
                if chrom == read_lb[0]:
                    if int(fend) == int(read_lb[1]) or int(rend) == int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{read_lb[4]}\t+\t-\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{frag_id}\t{dist}\t{read_lb[9]}\t{read_lb[12]}\n")
                    elif int(fend) < int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{read_lb[4]}\t{strand}\t{read_lb[5]}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{frag_id}\t{dist}\t{read_lb[9]}\t{read_lb[12]}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{read_lb[9]}\t{read_lb[12]}\t{frag_id}\t{dist}\n")
                else:
                    if chrom < read_lb[0]:
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{read_lb[4]}\t{strand}\t{read_lb[5]}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{frag_id}\t{dist}\t{read_lb[9]}\t{read_lb[12]}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{linkage[name_correct][0]}\t{linkage[name_correct][1]}\t{read_lb[9]}\t{read_lb[12]}\t{frag_id}\t{dist}\n")
        # 处理没有匹配成功的情况
        elif name_correct in unlinked:
            if read_lb[3] == name_correct:
                if chrom == read_lb[0]:
                    if int(fend) == int(read_lb[1]) or int(rend) == int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{read_lb[4]}\t+\t-\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{frag_id}\t{dist}\t{read_lb[9]}\t{read_lb[12]}\n")
                    elif int(fend) < int(read_lb[2]):
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{read_lb[4]}\t{strand}\t{read_lb[5]}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{frag_id}\t{dist}\t{read_lb[9]}\t{read_lb[12]}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{read_lb[9]}\t{read_lb[12]}\t{frag_id}\t{dist}\n")
                else:
                    if chrom < read_lb[0]:
                        fho.write(f"{chrom}\t{fend}\t{rend}\t{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{name}\t{mapq}\t{read_lb[4]}\t{strand}\t{read_lb[5]}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{frag_id}\t{dist}\t{read_lb[9]}\t{read_lb[12]}\n")
                    else:
                        fho.write(f"{read_lb[0]}\t{read_lb[1]}\t{read_lb[2]}\t{chrom}\t{fend}\t{rend}\t{name}\t{read_lb[4]}\t{mapq}\t{read_lb[5]}\t{strand}\t{unlinked[name_correct][0]}\t{unlinked[name_correct][1]}\t{read_lb[9]}\t{read_lb[12]}\t{frag_id}\t{dist}\n")
        read_lb = [chrom, fend, rend, name_correct, mapq, strand, chrom_frag, fend_frag, rend_frag, frag_id, space, frag_strand, dist]

def main():
    fhi = open(sys.argv[1])  # 输入 BED 文件
    barcodes = open(sys.argv[2])  # 输入标识码文件
    mismatch_out = open(sys.argv[3], 'w')  # 输出统计错误的文件
    outfile = open(sys.argv[4], 'w')  # 输出的文件
    uncouple_reads(fhi, barcodes, mismatch_out, outfile)
    fhi.close()
    barcodes.close()
    mismatch_out.close()
    outfile.close()

if __name__ == "__main__":
    main()
