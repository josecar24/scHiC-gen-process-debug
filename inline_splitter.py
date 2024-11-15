import os,sys,re
import itertools as it
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Levenshtein import hamming
from tqdm import tqdm
import gzip as gz


def checkHamming(barcodes,barcode):
	'''Given a list of barcodes, check that the given barcode is within edit\
		distance 2 to any of the list of barcodes
    给定一个标识码列表和一个目标标识码，检查此标识码是否在给定标识码列表中近似（编辑距离小于或等于2）
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

def split_fastqs(r1, r2, r1_o, r2_o, barcodes):
	'''Given paired-end fastq data, split reads based off of an inline 10 (or 11) mer barcode
    根据配对的子带标识码，分割或过滤配对的组对进行读取。
    :param r1: 输入的带有标识码的对应R1文件
    :param r2: 输入的带有标识码的对应R2文件
    :param r1_o: 输出的对应R1文件
    :param r2_o: 输出的对应R2文件
    :param barcodes: 标识码列表
    '''
	mismatch = 0
	not_found_R1 = 0
	not_found_R2 = 0
	reads = 0
	fqr1 = FastqGeneralIterator(r1)
	fqr2 = FastqGeneralIterator(r2) # 将两个进行器绑定在一起便于简化进程
	seqzip = zip(fqr1, fqr2) #Zip up the two iterators for expediency
 
	for pairs in tqdm(seqzip, desc="Processing reads", mininterval=100.0):
		title1, seq1, qual1 = pairs[0]
		title2, seq2, qual2 = pairs[1]
		barcode1 = seq1[:8] #Just look in read 1--barcodes SHOULD be the same on both ends of the molecule)  # 只测试R1序列的前8位
		barcode2 = seq2[:8] #Check barcode 2 as well; print out how many times they disagree # 同时测试R2序列的前8位；指出它们的不匹配情况
		test1 = checkHamming(barcodes, barcode1)
		test2 = checkHamming(barcodes, barcode2)
		if test1[0]:
			if test2[0]:
			# 如果两个序列的标识码均匹配，输出被修剪的序列
			#if the barcodes match, print out the trimmed / split reads to new files
				if test1[1] == test2[1]:
					print(f"@{title1}&{test1[1]}\n{seq1[11:]}\n+\n{qual1[11:]}", file=r1_o)
					print(f"@{title2}&{test2[1]}\n{seq2[11:]}\n+\n{qual2[11:]}", file=r2_o)
				else:
					mismatch += 1 # 如果R1和R2的标识码不匹配，进行记录
			else: #If there isn't a match in R1
				not_found_R2 += 1
		elif test2[0]:
			not_found_R1 += 1
		else:
			not_found_R1 += 1
			not_found_R2 += 1
		reads += 1
	out_error0 = "<H3>Total number of reads:%s</H3>" % reads
	out_error1 = "<H3>Total number of barcode mismatches:%s</H3>" % mismatch
	out_error2 = "<H3>Total number of missed R1 barcodes:%s</H3>" % not_found_R1
	out_error3 = "<H3>Total number of missed R2 barcodes:%s</H3>" % not_found_R2
	sys.stderr.write(out_error0 + '\n' + out_error1 + '\n' + out_error2 + '\n' + out_error3 + '\n')

def main():
	r1 = gz.open(sys.argv[1], 'rt') #R1 filehandle
	r2 = gz.open(sys.argv[2], 'rt') #R2 filehandle
	barcode_fhi = open(sys.argv[3],'r') #Barcode filehandle
	#Barcodes
	barcodes = []
	for line in barcode_fhi:
		barcodes.append(line.split()[0])
	r1_o = open(sys.argv[4], 'w') #R1 filehandle out
	r2_o = open(sys.argv[5], 'w') #R2 filehandle out
	split_fastqs(r1, r2, r1_o, r2_o, barcodes) #Split some fastqs, yo
	#Close the filehandles
	r1.close()
	r2.close()
	r1_o.close()
	r2_o.close()
	barcode_fhi.close()

if __name__ == "__main__":
	main()
