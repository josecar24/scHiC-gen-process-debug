#Script for analyzing sciHi-C fastqs.
import os,sys,re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools as it
import gzip as gz
from tqdm import tqdm
from Levenshtein import hamming

# 添加通配符
def addWildcards(string):
	foo= ""
	for i in range(len(string)):
        	foo+=string[:i] + '.' + string[i+1:] + "|"
	return foo[:-1]

# 实现反向互补的方法
def rc(seq):
	rcd = ""
	rcs = {"A":"T","C":"G","T":"A","G":"C", "N":"N"}
	for char in seq[::-1]:
		rcd+=rcs[char]
	return rcd

# 检查给定条形码的编辑距离是否小于或等于2
def checkHamming(barcodes,barcode):
	for bc in barcodes:
		match = False
		hd = hamming(barcode, bc)
		if hd <= 2:
			match = True
			barcode = bc
			break
	return (match, barcode)

# 关联条形码并根据内部条形码处理输入的fastq数据
#@profile
def associateBarcodes(fhi_r1, fhi_r2, fho_r1, fho_r2, internal_bc):
	'''Takes a Fastq filehandle in and filehandle out for each read; writes barcode-clipped fastqs to\
	filehandle out and returns a dictionary with barcode pairs associated with given read names'''
	seq_r1 = FastqGeneralIterator(fhi_r1) # 读取R1文件的序列
	seq_r2 = FastqGeneralIterator(fhi_r2)
	seqzip = zip(seq_r1, seq_r2) #zip up the two fastq iterators to increase performance
	bc_seq = {}
	bc1s = [] #Container to store internal barcode IDs
	bridge = ''
	#Create a list of all possible internal barcodes (necessary
	#for hamming distance calculations)
	for line in internal_bc:
		bc = line.rstrip('\n')
		strip = 'GATC' + bc + 'GAATTC' # 构造包含内部条形码的序列
		bc1s.append(bc)
    # 正则表达式匹配包含内部适配子的序列
	bridge_re = re.compile('GATC[ACGT]{8}GAATTC')
    #We want to find all reads where
    #the internal adaptor is present
    #in either R1 or R2. Accepting that we're
    #throwing out some stuff because of mismatches,
    ##let's just use the re module to do this.
    
    # 设置进度条
	total_reads = sum(1 for _ in seq_r1)
	fhi_r1.seek(0)  # 重置 R1 文件指针以重新读取
	fhi_r2.seek(0)  # 重置 R2 文件指针以重新读取
	seq_r1 = FastqGeneralIterator(fhi_r1)
	seq_r2 = FastqGeneralIterator(fhi_r2)
	seqzip = zip(seq_r1, seq_r2)
    
    # 处理每对序列
	for i in tqdm(seqzip, desc="Processing reads", total=total_reads, mininterval=100.0):
		title1, seq1, qual1 = i[0] # 解析R1的标题、序列和质量分数
		title2, seq2, qual2 = i[1] # 解析R2的标题、序列和质量分数
		new_title = title1.split()[0] # 从标题中提取基本信息
		m1 = bridge_re.search(seq1) # 查找R1中的适配子
		m2 = bridge_re.search(seq2) # 查找R2中的适配子
		terminal_bc = title1.split('&')[1] # 从标题中获取terminal条形码
		if m1: # 如果在R1中找到了适配子
			if m2: # 如果在R2中也找到了适配子
				pos1 = m1.start() # 获取R1中适配子的起始位置
				pos2 = m2.start() # 获取R2中适配子的起始位置
				trim_seq1 = seq1[:pos1]
				trim_seq2 = seq2[:pos2]
				trim_qual1 = qual1[:pos1]
				trim_qual2 = qual2[:pos2]
				#Check for length of the read  如果裁剪后的序列太短，跳过
				if len(trim_seq1) < 30 or len(trim_seq2) < 30:
					continue 
 
                # 获取序列中条形码位置的序列
				barcode_seq1 = seq1[pos1 + 4:pos1 + 12]
				barcode_seq2 = seq2[pos2 + 4:pos2 + 12]
                # 检查内部条形码的编辑距离是否小于等于2
				ham_bc1 = checkHamming(bc1s, barcode_seq1)
				ham_bc2 = checkHamming(bc1s, barcode_seq2)
    			# 如果两个条形码均匹配
				if ham_bc1[0]:
					if ham_bc2[0]:
						print(f"@{new_title}\n{trim_seq1}\n+\n{trim_qual1}", file=fho_r1)
						print(f"@{new_title}\n{trim_seq2}\n+\n{trim_qual2}", file=fho_r2)
						bc1 = ham_bc1[1]
						bc2 = ham_bc2[1]
						print(f"{new_title}\t{bc1}\t{bc2}\t{terminal_bc}")
			else: # 如果只在R1中找到了适配子
				pos1 = m1.start()
				pos1_end = m1.end()
				trim_seq1 = seq1[:pos1]
				trim_qual1 = qual1[:pos1]
				#Check for length of the rea
				if len(trim_seq1) < 30:
					continue #and continue if either one is too short
				if len(seq1[pos1_end:]) < 12:
					continue
				# 反向互补获取条形码
				barcode_seq2 = rc(seq1[pos1_end:pos1_end + 8])
				barcode_seq1 = seq1[pos1 + 4:pos1 + 12]
				ham_bc1 = checkHamming(bc1s, barcode_seq1)
				ham_bc2 = checkHamming(bc1s, barcode_seq2)
				if ham_bc1[0]:
					if ham_bc2[0]:
						print(f"@{new_title}\n{trim_seq1}\n+\n{trim_qual1}", file=fho_r1)
						print(f"@{new_title}\n{seq2}\n+\n{qual2}", file=fho_r2)
						bc1 = ham_bc1[1]
						bc2 = ham_bc2[1]
						print(f"{new_title}\t{bc1}\t{bc2}\t{terminal_bc}")
		else: # 如果只在R2中找到了适配子
			if m2:
				pos2 = m2.start()
				pos2_end = m2.end()
				trim_seq2 = seq2[:pos2]
				trim_qual2 = qual2[:pos2]
				#Check for length of the read
				#print seq2[pos2_end:]
				if len(trim_seq2) < 30:
					continue #and continue if either one is too short
				if len(seq2[pos2_end:]) < 12:
					continue #and continue if either one is too short
				# 反向互补获取条形码
				barcode_seq1 = rc(seq2[pos2_end:pos2_end + 8])
				barcode_seq2 = seq2[pos2 + 4:pos2 + 12]
				ham_bc1 = checkHamming(bc1s, barcode_seq1)
				ham_bc2 = checkHamming(bc1s, barcode_seq2)
    
				if ham_bc1[0]:
					if ham_bc2[0]:
						print(f"@{new_title}\n{seq1}\n+\n{qual1}", file=fho_r1)
						print(f"@{new_title}\n{trim_seq2}\n+\n{trim_qual2}", file=fho_r2)
						bc1 = ham_bc1[1]
						bc2 = ham_bc2[1]
						print(f"{new_title}\t{bc1}\t{bc2}\t{terminal_bc}")

def main():
	fhi_bc = open(sys.argv[1])
	fhi_r1 = open(sys.argv[2])
	fhi_r2 = open(sys.argv[3])
	fho_r1 = open(sys.argv[4], 'w')
	fho_r2 = open(sys.argv[5], 'w')
	assoc = associateBarcodes(fhi_r1, fhi_r2, fho_r1, fho_r2, fhi_bc)
	fhi_bc.close()
	fhi_r1.close()
	fhi_r2.close()
	fho_r1.close()
	fho_r2.close()

if __name__ == "__main__":
	main()
