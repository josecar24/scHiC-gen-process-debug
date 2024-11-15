# Goal: Given a "percentages" file which specifies the
# species breakdown and coverage within each barcode combination,
# generate sparse matrix files at specified resolutions for each single cell
import os
import sys
import re
from collections import Counter
from tqdm import tqdm
from math import sqrt

def normalize_matrix(matrix):
    '''Given a matrix, root normalize the values in the matrix'''
    cov = Counter()  # Store coverage of each bin in a Counter object
    normed = {}  # Store the normalized values in a dict
    for i in matrix:
        # Iterate through the matrix and keep track of the counts of each bin
        bin1, bin2, chrom1, chrom2 = i
        count = matrix[i]
        cov[bin1] += count
        cov[bin2] += count
    for i in matrix:
        # Walk through the matrix once again and calculate the normalized values
        # based on the computed coverages
        bin1, bin2, chrom1, chrom2 = i
        normed[i] = float(matrix[i]) / (sqrt(cov[bin1]) * sqrt(cov[bin2]))
    return normed

def define_bins(chromsizes, resolutions):
    '''Define the bins that cover the genome at the desired resolution(s)'''
    bins = {}  # Store the bins
    valid_chroms = {}  # Store valid chromosome names
    lines = chromsizes.readlines()  # Read chromosome sizes file as a list

    for resolution in resolutions:
        bins[resolution] = {}  # Each resolution has its own dictionary

    for resolution in resolutions:
        hindex = 0  # Index for human chromosomes
        for line in lines:
            # Split the line and ensure we only use the first two columns: chromname and length
            split_line = line.split()
            chromname, length = split_line[0], split_line[1]
            valid_chroms[chromname] = True
            for i in range(0, int(length), resolution):
                bins[resolution][(chromname, i)] = hindex
                hindex += 1

    return bins, valid_chroms

def cell_sort(percentages):
    '''Return a list of cell barcodes where the cell has >=1000 unique reads, 
       passes cistrans_cutoff and is >= 95% from one species'''
    cells = {}
    valid_count = 0
    for line in percentages:
        # Destructure the line using human-only format
        hup, huc, val, tot, bc1, bc2, sig, length, ones, twos, threes, fours, ratio = line.split()

        if length == "Long":
            continue  # Ignore "Long" rows
        if sig == "Randomized":
            continue  # Ignore "Randomized" rows
        if float(hup) >= 0.95 and int(huc) >= 1000 and float(ratio) >= 1:
            identity = ("human", huc)
            barcode = f"{bc1}-{bc2}"
            cells[barcode] = identity  # Store species type and coverage with the cell barcode
            valid_count += 1
            
    # Debug: Print the number of valid cells
    print(f"[DEBUG] Number of valid cells: {valid_count}", file=sys.stderr)

    return cells

#print(f"[DEBUG] Number of valid cells: {len(cell_list)}", file=sys.stderr)


def bedpe_walk(bedpe, cell_list, resolutions, bins, valid_chroms):
    '''Walk through a bedpe file, generate sparse matrices'''
    cell_matrices = {resolution: {} for resolution in resolutions}
    
    valid_entries = 0
    
    total_lines = 0
    processed_lines = 0

    for line in bedpe:
        total_lines += 1
        n1, f1, r1, n2, f2, r2, name, q1, a2, s1, s2, bc1, bc2, frag1, dist1, frag2, dist2, dupcount = line.split()

        species1 = n1.split("_")[0]
        species2 = n2.split("_")[0]
        barcode = f"{bc1}-{bc2}"

        if barcode in cell_list:
            if species1 != species2:
                continue  # Skip mismatched or contaminating species reads
            if n1 not in valid_chroms or n2 not in valid_chroms:
                continue  # Skip invalid chromosomes
            
            valid_entries += 1

            for resolution in resolutions:
                if barcode not in cell_matrices[resolution]: # 如果对于当前条形码（即 barcode）和分辨率还没有交互记录，那么需要为这个细胞的交互矩阵建立一个空的计数器，以便后续记录每个 bin 对之间的交互次数。
                    cell_matrices[resolution][barcode] = Counter() 

                # Calculate bin positions
                pos1_reduce = (int(f1) + int(r1)) // 2 // resolution * resolution  # Use the fragment midpoint 这行代码计算了第一个 read 的位置，目的是找出该 read 所在的 bin
                pos2_reduce = (int(f2) + int(r2)) // 2 // resolution * resolution  # Use the fragment midpoint 它是对第二个 read 进行计算，用于找出第二个 read 所在的 bin
                bin1 = bins[resolution][(n1, pos1_reduce)] # 这行代码将第一个 read 的染色体编号 n1 和位置 pos1_reduce 映射到 bins 字典中，以找到它对应的 bin ID
                bin2 = bins[resolution][(n2, pos2_reduce)]

                # Only report bins where bin1 <= bin2 to save space
                if bin1 <= bin2:
                    key = (bin1, bin2, n1, n2)
                else:
                    key = (bin2, bin1, n2, n1)

                cell_matrices[resolution][barcode][key] += 1
            processed_lines += 1
                
    print(f"[DEBUG] Number of valid bedpe entries processed: {valid_entries}", file=sys.stderr)

    return cell_matrices

def main():
    # 使用上下文管理器打开文件
    with open(sys.argv[1]) as genome_file, open(sys.argv[2]) as percentages, open(sys.argv[3]) as bedpe:
        # User defined resolutions
        resolutions = [500000]

        # 定义染色体的 bins
        bins, valid_chroms = define_bins(genome_file, resolutions)

        # 筛选出符合条件的单细胞条形码
        cell_list = cell_sort(percentages)
        
        # 在这里添加调试语句来查看 cell_list 的数量
        print(f"[DEBUG] Number of valid cells: {len(cell_list)}", file=sys.stderr)

        # 处理 BEDPE 文件生成矩阵
        cell_matrices = bedpe_walk(bedpe, cell_list, resolutions, bins, valid_chroms)

        # 对每个条形码和分辨率输出矩阵
        #for resolution in resolutions:
            #for barcode in cell_matrices[resolution]:
                #fho_name = f"{cell_list[barcode][0]}_{cell_list[barcode][1]}_{barcode}_{resolution}.matrix"
                #with open(fho_name, 'w') as fho:
                    #norm = normalize_matrix(cell_matrices[resolution][barcode])
                    #for i in norm:
                        #fho.write(f"{i[0]}\t{i[1]}\t{cell_matrices[resolution][barcode][i]}\t{norm[i]}\t{i[2]}\t{i[3]}\n")
                        
                # 输出矩阵到指定的目录
        output_dir = "1018_out_dir/br1_tr1_pair1.matrices"
        os.makedirs(output_dir, exist_ok=True)

        valid_cells_list_path = os.path.join(output_dir, "valid_cells.list")
        with open(valid_cells_list_path, 'w') as valid_cells_list:
            # 对每个条形码和分辨率输出矩阵
            for resolution in resolutions:
                for barcode in tqdm(cell_matrices[resolution], desc=f"Generating matrices at resolution {resolution}"):
                    fho_name = f"{cell_list[barcode][0]}_{cell_list[barcode][1]}_{barcode}_{resolution}.matrix"
                    fho_path = os.path.join(output_dir, fho_name)

                    # 将矩阵文件路径写入到 valid_cells.list
                    valid_cells_list.write(f"{fho_name}\n")

                    # 输出矩阵文件
                    with open(fho_path, 'w') as fho:
                        norm = normalize_matrix(cell_matrices[resolution][barcode])
                        for i in norm:
                            fho.write(f"{i[0]}\t{i[1]}\t{cell_matrices[resolution][barcode][i]}\t{norm[i]}\t{i[2]}\t{i[3]}\n")


if __name__ == "__main__":
    main()

