import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# 打开两个文件
with open('./1018_tmp/4DNFI7HX4GMM.fastq.gz.bc_clipped', 'r') as r1, open('./1018_tmp/4DNFIF88P2F1.fastq.gz.bc_clipped', 'r') as r2:
    # 使用 FastqGeneralIterator 逐个读取序列
    fqr1 = FastqGeneralIterator(r1)
    fqr2 = FastqGeneralIterator(r2)
    
    mismatch_count = 0
    reads = 0

    for (title1, seq1, qual1), (title2, seq2, qual2) in zip(fqr1, fqr2):
        barcode1 = title1.split('&')[1]
        barcode2 = title2.split('&')[1]
        if barcode1 != barcode2:
            mismatch_count += 1
            print(f"Mismatch in read pair {reads + 1}: R1 barcode {barcode1} != R2 barcode {barcode2}")
        reads += 1

    print(f"Total reads checked: {reads}")
    print(f"Total mismatches found: {mismatch_count}")