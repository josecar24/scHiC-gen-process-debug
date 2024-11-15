import csv
from tqdm import tqdm  # 导入 tqdm 模块用于进度条

input_file = "1018_out_dir/br1_tr1_pair1.concat.sorted.bed"
output_file = "1018_out_dir/br1_tr1_pair1.concat.sorted.modified.bed"

# 统计总行数用于进度条
with open(input_file, 'r') as f_in:
    total_lines = sum(1 for _ in f_in)

# 打开原始文件，添加三列占位符信息
with open(input_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
    reader = csv.reader(f_in, delimiter='\t')
    writer = csv.writer(f_out, delimiter='\t')
    
    # 使用 tqdm 添加进度条
    for row in tqdm(reader, total=total_lines, desc="Processing Rows", mininterval=1.0):
        if len(row) == 10:  # 修改条件为 10 列
            # 在第 9 列和第 10 列之间插入三列信息
            row = row[:9] + ['placeholder1', 'placeholder2', 'placeholder3'] + row[9:]
        writer.writerow(row)

print(f"File {output_file} has been created with 13 columns.")

