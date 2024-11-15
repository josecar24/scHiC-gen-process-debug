import os
import sys
import re

def define_valid_chroms(genome_file):
    # Create a dictionary of valid chromosomes from the genome file
    valid_chroms = {}
    for line in genome_file:
        split = line.split()
        valid_chroms[split[0]] = True #使用列表中的第一个元素 split[0] 作为键，把它添加到 valid_chroms 字典中
    return valid_chroms

def filter_bed(bedpe, valid_chroms):
    # Filter the BEDPE file, printing only entries with valid chromosomes
    for line in bedpe:
        split = line.split()
        chrid1 = split[0]
        chrid2 = split[3]
        if chrid1 in valid_chroms and chrid2 in valid_chroms:
            print(line.rstrip('\n'))

def main():
    # Open genome and BEDPE files passed as arguments
    with open(sys.argv[1]) as genome_file:
        with open(sys.argv[2]) as bedpe:
            valid_chroms = define_valid_chroms(genome_file)
            filter_bed(bedpe, valid_chroms)

if __name__ == "__main__":
    main()

