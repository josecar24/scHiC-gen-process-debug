import os
import sys
import json


# 使用上下文管理器打开第一个文件：cis-trans
with open(sys.argv[1]) as fhi:
    ratios = {}

    # 读取并存储每个条形码的比值
    for line in fhi:
        split = line.split()
        barcode = split[0]
        ratio_value = split[1]

        # 确保数据被正确读取
        print(f"[DEBUG] Read from cis-trans: Barcode: {barcode}, Ratio: {ratio_value}", file=sys.stderr)

        ratios[barcode] = ratio_value
        
# 转储 ratios 字典到一个调试文件中
#with open("ratios_dump.json", "w") as dump_file:
#    json.dump(ratios, dump_file, indent=4)
#    print("[DEBUG] Ratios have been dumped to ratios_dump.json", file=sys.stderr)

# 初始化计数器
found_count = 0
zero_count = 0
num_ex = 0
total_count = 0  # 总计数器，用于统计百分比文件中的总行数

# 使用上下文管理器打开第二个文件：percentages
with open(sys.argv[2]) as fhi:
    # 处理 percentages 文件中的每一行

    for line in fhi:
        total_count += 1  # 每处理一行总数加一
        split = line.split()
        barcode = f"{split[4]}-{split[5]}"

        # 确保百分比文件中的条形码也被正确读取
        print(f"[DEBUG] Processing percentages file: Barcode: {barcode}", file=sys.stderr)
        # 判断 ratios 中是否已经存在对应条形码
        if barcode in ratios and ratios[barcode] != 0:
            # 如果条形码存在且其值不是 0
            print(f"[DEBUG] Found barcode in cis-trans, using ratio: {ratios[barcode]}", file=sys.stderr)
            found_count += 1
        else:
            # 如果 ratios 中没有对应的条形码，设为 0
            print(f"[DEBUG] Barcode not found in cis-trans, setting ratio to 0: {barcode}", file=sys.stderr)
            ratios[barcode] = 0
            zero_count += 1  # 统计被设为 0 的条形码数量

        # 打印每一行的数据，增加 cis-trans 比例
        print(f"{line.rstrip()}\t{ratios[barcode]}")

# 输出总计数结果
print(f"[SUMMARY] Total barcodes found in cis-trans: {found_count}", file=sys.stderr)
print(f"[SUMMARY] Total barcodes set to 0: {zero_count}", file=sys.stderr)
print(f"[SUMMARY] Total barcodes processed from percentages: {total_count}", file=sys.stderr)
print(f"[SUMMARY] Exs: {num_ex}", file=sys.stderr)


