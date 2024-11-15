import sys
import random

def sample_bed(input_bed, output_bed, sample_rate=0.01):
    """
    以指定的采样率从输入BED文件中抽取子样本，确保成对的Read保持在一起。
    
    :param input_bed: 原始BED文件路径
    :param output_bed: 输出子样本BED文件路径
    :param sample_rate: 采样率（默认1%）
    """
    with open(input_bed, 'r') as fin, open(output_bed, 'w') as fout:
        pair = []
        line_count = 0
        sampled_pairs = 0
        total_pairs = 0
        for line in fin:
            pair.append(line)
            if len(pair) == 2:
                total_pairs += 1
                if random.random() < sample_rate:
                    fout.writelines(pair)
                    sampled_pairs += 1
                pair = []
        # 处理文件中行数为奇数的情况
        if pair:
            print("警告：输入文件中存在未成对的行，将被忽略。")
        
        print(f"总对数: {total_pairs}, 采样对数: {sampled_pairs}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python sample_bed.py <输入BED文件> <输出子样本BED文件>")
        sys.exit(1)
    
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    sample_bed(input_bed, output_bed)
