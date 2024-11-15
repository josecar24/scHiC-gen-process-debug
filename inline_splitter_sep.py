import os, sys, re
import itertools as it
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Levenshtein import hamming
import gzip as gz
from tqdm import tqdm

# 模块和函数实现，用于应用变程和参数的处理

def checkHamming(barcodes, barcode):
    '''
    根据标识码列表，检查目标标识码是否近似于列表中的标识码（编辑距离小于或等于2）
    :param barcodes: 用来比对的标识码列表
    :param barcode: 待检查的标识码
    :return: 返回是否匹配和匹配的标识码
    '''
    for bc in barcodes:
        match = False
        hd = hamming(barcode, bc)
        if hd <= 2:
            match = True
            barcode = bc
            break
    return (match, barcode)

# 功能：给定配对的快速库数据，根据标识码进行分割

def split_fastqs(r1, r2, r1_o, r2_o, barcodes, total_reads):
    '''
    根据配对的子带标识码，分割或过滤配对组对进行读取。
    :param r1: 输入的带有标识码的对应 R1 文件
    :param r2: 输入的带有标识码的对应 R2 文件
    :param r1_o: 输出的对应 R1 文件
    :param r2_o: 输出的对应 R2 文件
    :param barcodes: 标识码列表
    :param total_reads: 总读次数
    '''
    mismatch = 0
    not_found_R1 = 0
    not_found_R2 = 0
    reads = 0
    fqr1 = FastqGeneralIterator(r1)
    fqr2 = FastqGeneralIterator(r2)  # 将两个迭代器绑定在一起，便于简化进程
    seqzip = it.zip_longest(fqr1, fqr2)  # 将两个迭代器编组为对象
    for pairs in tqdm(seqzip, total=total_reads, desc="Processing reads"):
        title1, seq1, qual1 = pairs[0]
        title2, seq2, qual2 = pairs[1]
        barcode1 = seq1[:8]  # 只检测 R1 序列的前 8 位
        barcode2 = seq2[:8]  # 同时检测 R2 序列的前 8 位
        test1 = checkHamming(barcodes, barcode1)
        test2 = checkHamming(barcodes, barcode2)
        
        if test1[0]:
            print(f"Matched barcode in R1: {barcode1}")
        if test2[0]:
            print(f"Matched barcode in R2: {barcode2}")
        
        if test1[0]:
            if test2[0]:
                # 如果两个序列的标识码均匹配，输出被修剪的序列
                if test1[1] == test2[1]:
                    r1_o.write("@%s&%s\n%s\n+\n%s\n" % (title1, test1[1], seq1[11:], qual1[11:]))
                    r2_o.write("@%s&%s\n%s\n+\n%s\n" % (title2, test2[1], seq2[11:], qual2[11:]))
                else:
                    mismatch += 1  # 如果 R1 和 R2 的标识码不匹配，进行记录
            else:  # R2 的标识码不匹配
                not_found_R2 += 1
        elif test2[0]:
            not_found_R1 += 1
        else:
            not_found_R1 += 1
            not_found_R2 += 1
        reads += 1
    # 返回结果统计信息
    out_error0 = "<H3>Total number of reads:%s</H3>" % reads
    out_error1 = "<H3>Total number of barcode mismatches:%s</H3>" % mismatch
    out_error2 = "<H3>Total number of missed R1 barcodes:%s</H3>" % not_found_R1
    out_error3 = "<H3>Total number of missed R2 barcodes:%s</H3>" % not_found_R2
    sys.stderr.write(out_error0 + '\n' + out_error1 + '\n' + out_error2 + '\n' + out_error3 + '\n')

# 主函数：进行模块的运行

def main():
    r1 = gz.open(sys.argv[1], 'rt')  # R1 文件的输入处理
    r2 = gz.open(sys.argv[2], 'rt')  # R2 文件的输入处理
    barcode_fhi = open(sys.argv[3])  # 标识码文件的输入

    # 读取标识码文件，并根据细胞系进行分类
    barcodes_gm12878 = []
    barcodes_hff = []
    for line in barcode_fhi:
        barcode, cell_type = line.split() # 默认分割
        if cell_type == "GM12878":
            barcodes_gm12878.append(barcode)
        elif cell_type == "HFFc6":
            barcodes_hff.append(barcode)
            
    # 打印条形码加载情况
    print(f"GM12878 barcodes loaded: {len(barcodes_gm12878)}")
    print(f"HFFc6 barcodes loaded: {len(barcodes_hff)}")

    # 创建压缩的输出文件
    r1_o_gm = gz.open(sys.argv[4], 'wt')  # GM12878 R1 输出（压缩）
    r2_o_gm = gz.open(sys.argv[5], 'wt')  # GM12878 R2 输出（压缩）
    r1_o_hff = gz.open(sys.argv[6], 'wt')  # HFFc6 R1 输出（压缩）
    r2_o_hff = gz.open(sys.argv[7], 'wt')  # HFFc6 R2 输出（压缩）

    # 计算总读次数
    total_reads = sum(1 for _ in FastqGeneralIterator(r1))
    r1.seek(0)  # 重置文件读取位置
    r2.seek(0)

    # 分别进行拆分
    split_fastqs(r1, r2, r1_o_gm, r2_o_gm, barcodes_gm12878, total_reads)  # 拆分 GM12878 数据
    split_fastqs(r1, r2, r1_o_hff, r2_o_hff, barcodes_hff, total_reads)  # 拆分 HFFc6 数据

    # 关闭文件
    r1.close()
    r2.close()
    r1_o_gm.close()
    r2_o_gm.close()
    r1_o_hff.close()
    r2_o_hff.close()
    barcode_fhi.close()

if __name__ == "__main__":
    main()
