import os
import sys
import re
from collections import Counter

# 使用上下文管理器打开文件，确保其在处理完成后正确关闭
with open(sys.argv[1]) as fhi:
    barcodes = {}

    for line in fhi:
        split = line.split('\t')
        barcode = split[1]

        # 如果条形码不在字典中，则初始化它
        if barcode not in barcodes:
            barcodes[barcode] = Counter()

        # 对不同类型的距离进行处理
        if split[2] in [">20kb", "<1kb", "1kb-20kb"]:
            barcodes[barcode][">20kb"] += float(split[3]) / float(split[4])
        
        if split[2] == "Interchromosomal":
            barcodes[barcode]["Interchromosomal"] = float(split[3]) / float(split[4])

# 计算并输出每个条形码的 cistrans 比例
for barcode in barcodes:
    if barcodes[barcode]["Interchromosomal"] != 0:  # 确保分母不为 0
        cistrans_ratio = barcodes[barcode][">20kb"] / barcodes[barcode]["Interchromosomal"]
    else:
        cistrans_ratio = "undefined"  # 避免除零错误，返回未定义的比例

    print(f"{barcode}\t{cistrans_ratio}")
