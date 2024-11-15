import os
import sys
import re
from collections import Counter

def cell_sort(percentages):
    '''Iterate through the percentages file generated by the scHi-C pipeline, and return a list
       of cell barcodes where the cell has >1000 unique reads and is 0.95 one species'''
    cells = {}
    ids = {}
    for line in percentages:
        # 对每个条形码组合进行处理
        hup, huc, val, tot, bc1, bc2, sig, length, ones, twos, threes, fours = line.split()

        # 忽略长距离和随机化行
        if length == "Long": 
            continue
        if sig == "Randomized": 
            continue

        # 筛选读取数大于 1000 的条形码组合
        if int(huc) >= 1000: 
            identity = ("human", huc)
            barcode = f"{bc1}-{bc2}"
            cells[barcode] = Counter()  # 将条形码组合存储在细胞字典中
            ids[barcode] = identity
    
    return (cells, ids)

def single_cell_contacts(fhi, cell_contacts):
    # 对每个条形码组合的交互信息进行分类
    for line in fhi:
        split = line.split()
        bc1 = split[11]
        bc2 = split[12]
        barcode = f"{bc1}-{bc2}"

        # 如果条形码在筛选后的单细胞条形码列表中，则处理该条形码
        if barcode in cell_contacts:
            fcoord1 = int(split[1])
            rcoord2 = int(split[5])
            chr1 = split[0]
            chr2 = split[3]

            # 只处理相同物种的交互（在本例中只有人类）
            if chr1 != chr2:
                cell_contacts[barcode]["Interchromosomal"] += 1
            else:
                dist = abs(fcoord1 - rcoord2)
                if dist > 20000:
                    cell_contacts[barcode][">20kb"] += 1
                elif 1000 <= dist <= 20000:
                    cell_contacts[barcode]["1kb-20kb"] += 1
                else:
                    cell_contacts[barcode]["<1kb"] += 1

    return cell_contacts

def main():
    # 使用上下文管理器打开文件，保证文件可以自动关闭
    with open(sys.argv[1]) as percentages:
        cells = cell_sort(percentages)

    with open(sys.argv[2]) as bedpe:
        cell_contacts = single_cell_contacts(bedpe, cells[0])

    # 逐行输出条形码组合的交互信息
    for barcode in cell_contacts:
        for entry in cell_contacts[barcode]:
            print(f"{sys.argv[3]}\t{barcode}\t{entry}\t{cell_contacts[barcode][entry]}\t{sum(cell_contacts[barcode].values())}\t{cells[1][barcode][0]}")

if __name__ == "__main__":
    main()

