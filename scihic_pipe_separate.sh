#全局变量
# 输入参数（根据您的命令）
bc1="inner_barcodes.txt"                # Inner barcode file
fq_r1="4DNFI7HX4GMM.fastq.gz"           # R1 Fastq
fq_r2="4DNFIF88P2F1.fastq.gz"           # R2 Fastq
barcodes="outer_barcodes.txt"           # Outer barcode file
bc_assoc="br1_tr1_pair1"                # Outfile prefix
outdir="1018_out_dir"                   # Output directory

# 临时目录
tmp_dir="./1018_tmp/"
mkdir -p $tmp_dir
mkdir -p $outdir

#1
# 使用SeqPrep进行适配子剪切
#SeqPrep -A AGATCGGAAGAGCGATCGG -B AGATCGGAAGAGCGTCGTG \  #这里要用fastqc找到的正确seq来代替,虽然Adapter Content没有直接给序列，但是看图应该没有多少序列污染
SeqPrep -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -f ./1018_tmp/4DNFI7HX4GMM.fastq.gz \
    -r ./1018_tmp/4DNFIF88P2F1.fastq.gz \
    -1 ./1018_tmp/4DNFI7HX4GMM.fastq.gz.clipped_seqp2 \
    -2 ./1018_tmp/4DNFIF88P2F1.fastq.gz.clipped_seqp2 \
    > 1018_out_dir/seq_prep.txt 2>> 1018_out_dir/adaptor_clipping_stats

#由于SeqPrep这个工具处理得非常慢且不能显示进度条不能并行，我是不是可以采取其他替代品？cutadapt可以多核心并行
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -o ./1018_tmp/4DNFI7HX4GMM.fastq.gz.clipped \
         -p ./1018_tmp/4DNFIF88P2F1.fastq.gz.clipped \
         ./1018_tmp/4DNFI7HX4GMM.fastq.gz \
         ./1018_tmp/4DNFIF88P2F1.fastq.gz \
         -j 8 > 1018_out_dir/seq_prep_cutadapt.txt 2>> 1018_out_dir/adaptor_clipping_stats


#2
#从GM12878 or HFFc6混合样本中分离 GM12878 细胞系
python inline_splitter_sep.py \
    ./1018_tmp/4DNFI7HX4GMM.fastq.gz.clipped \
    ./1018_tmp/4DNFIF88P2F1.fastq.gz.clipped \
    ./barcode/4DNFIFOSOLDM.txt \
    ./1018_tmp/GM12878_R1.fastq \
    ./1018_tmp/GM12878_R2.fastq \
    ./1018_tmp/HFFc6_R1.fastq \
    ./1018_tmp/HFFc6_R2.fastq \
    2> 1018_out_dir/splitting_stats.html

python inline_splitter_sep.py \ 
    ./1018_tmp/4DNFI7HX4GMM.fastq.gz \
    ./1018_tmp/4DNFIF88P2F1.fastq.gz \
    ./barcode/4DNFIFOSOLDM.txt \
    ./1018_tmp/GM12878_R1.fastq.gz \
    ./1018_tmp/GM12878_R2.fastq.gz \
    ./1018_tmp/HFFc6_R1.fastq.gz \
    ./1018_tmp/HFFc6_R2.fastq.gz \
    2> 1018_out_dir/splitting_stats.html
#目前来看上面那段代码(#2)没有作用

#3.1
# 使用inline_splitter.py进行按外部条形码拆分
python inline_splitter.py \
    ./1018_tmp/4DNFI7HX4GMM.fastq.gz \
    ./1018_tmp/4DNFIF88P2F1.fastq.gz \
    ./barcode/4DNFINE9JUAI.txt \
    ./1018_tmp/4DNFI7HX4GMM.fastq.gz.split \
    ./1018_tmp/4DNFIF88P2F1.fastq.gz.split \
    2> 1018_out_dir/splitting_stats.html

#3.2
# 使用analyze_scDHC_V2design.py处理内部条形码
python analyze_scDHC_V2design.py \
    ./barcode/4DNFIUSK42KP.txt \
    ./1018_tmp/4DNFI7HX4GMM.fastq.gz.split \
    ./1018_tmp/4DNFIF88P2F1.fastq.gz.split \
    ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bc_clipped \
    ./1018_tmp/4DNFIF88P2F1.fastq.gz.bc_clipped \
    > 1018_out_dir/br1_tr1_pair1.analyze_scDHC.log

#4
# 比对R1
bowtie2 -x combo_hg38 -p 8 -U ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bc_clipped \
    2> 1018_out_dir/4DNFI7HX4GMM.fastq.gz.mapping_stats \
    | samtools view -bS - > 1018_tmp/4DNFI7HX4GMM.fastq.gz.bam

# 比对R2
bowtie2 -x combo_hg38 -p 8 -U ./1018_tmp/4DNFIF88P2F1.fastq.gz.bc_clipped \
    2> 1018_out_dir/4DNFIF88P2F1.fastq.gz.mapping_stats \
    | samtools view -bS - > 1018_tmp/4DNFIF88P2F1.fastq.gz.bam

#建议的检查代码
samtools flagstat 1018_tmp/4DNFI7HX4GMM.fastq.gz.bam
samtools flagstat 1018_tmp/4DNFIF88P2F1.fastq.gz.bam

#(原脚本中的参考工序，现在应该不用管了)
wc -l ./1018_tmp/4DNFIF88P2F1.fastq.gz.bc_clipped > ./1018_tmp/br1_tr1_pair1.ph1
read reads_r1 junk < ./1018_tmp/br1_tr1_pair1.ph1
total_reads=`expr $reads_r1 / 4`

#5
# 将BAM文件转换为BED格式
bedtools bamtobed -i ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bam > ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bed
bedtools bamtobed -i ./1018_tmp/4DNFIF88P2F1.fastq.gz.bam > ./1018_tmp/4DNFIF88P2F1.fastq.gz.bed

#6
# 过滤映射质量（MAPQ > 30）的reads
awk '$5 > 30' ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bed > 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30
awk '$5 > 30' ./1018_tmp/4DNFIF88P2F1.fastq.gz.bed > 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30

#7
# 对BED文件进行排序
sort -k1,1 -k2,2n 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30 > 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted
sort -k1,1 -k2,2n 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30 > 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted

# 提取匹配的 reads ID
comm -12 <(cut -f4 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted | sort) <(cut -f4 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted | sort) > common_ids.txt

# 使用这些 ID 过滤原始的 BED 文件，只保留匹配的 reads
awk 'NR==FNR{a[$1];next} $4 in a' common_ids.txt 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted > 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered
awk 'NR==FNR{a[$1];next} $4 in a' common_ids.txt 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted > 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered

#再排序一次
# 使用 sort -k1,1 -k2,2n 进行排序，标识为 ".sorted.by_position"
sort -k1,1 -k2,2n 4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered > 4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered.by_position
sort -k1,1 -k2,2n 4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered > 4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered.by_position

# 使用 sort -k4,4 进行排序，标识为 ".sorted.by_readID"
sort -k4,4 4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered > 4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered.by_readID
sort -k4,4 4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered > 4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered.by_readID

# 检查使用 sort -k1,1 -k2,2n 排序的文件（后续用这个）
head -n 10 4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered.by_position
head -n 10 4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered.by_position

# 检查使用 sort -k4,4 排序的文件(在这里主要是看一看，起到检查作用)
head -n 10 4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered.by_readID
head -n 10 4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered.by_readID


#在此之前先用reference.fa和find_dpnii.py来生成combo_hg38.dpnii.bed (1105补丁：改用源码digest_genome.py来做)
python digest_genome.py reference.fa -r dpnii -o combo_hg38_ref.dpnii.bed

#再去掉杂质
grep -v "^KI" combo_hg38_ref.dpnii.bed | grep -v "^GL" > combo_hg38.filtered.dpnii.bed

#8 使用bedtools closest找到最近的DpnII位点
bedtools closest -t first -d \
    -a 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.sorted.filtered.by_position \
    -b combo_hg38.filtered.dpnii.bed \
    > 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.dre

bedtools closest -t first -d \
    -a 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.sorted.filtered.by_position \
    -b combo_hg38.filtered.dpnii.bed \
    > 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.dre



#9
# 合并并排序
cat 1018_out_dir/4DNFI7HX4GMM.fastq.gz.bed.mapq30.dre 1018_out_dir/4DNFIF88P2F1.fastq.gz.bed.mapq30.dre > 1018_out_dir/br1_tr1_pair1.concat.bed
sort -s -k 4 1018_out_dir/br1_tr1_pair1.concat.bed > 1018_out_dir/br1_tr1_pair1.concat.sorted.bed

#在转换为bedpe之前，先提取小样本测试一下代码?
#shuf -n 700000 1018_out_dir/br1_tr1_pair1.concat.sorted.bed > 1018_out_dir/br1_tr1_pair1.concat.sorted.sample.bed
python sample.py 1018_out_dir/br1_tr1_pair1.concat.sorted.bed 1018_out_dir/br1_tr1_pair1.concat.sorted.m.samp.bed


#把目前10列版本的文档填充成13列
#python fill_bed.py

# 将BED文件转换为BEDPE格式
python catbed2bedpe_2.py \
    1018_out_dir/br1_tr1_pair1.concat.sorted.bed \
    1018_out_dir/br1_tr1_pair1.analyze_scDHC.log \
    1018_out_dir/br1_tr1_pair1.all.bedpe \
    1018_out_dir/br1_tr1_pair1.matching.html \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq0 \
    1018_out_dir/br1_tr1_pair1.bedpe.mapunmatch \
    > 1018_out_dir/br1_tr1_pair1.bedpe.mapq_nan 
# 注意这里生成的文件是抽样版的，后面也是；调通之后记得覆盖

# 将BED文件转换为BEDPE格式
#python catbed2bedpe.py \
#    1018_out_dir/br1_tr1_pair1.concat.sorted.bed \
#    1018_out_dir/br1_tr1_pair1.analyze_scDHC.log \
#    1018_out_dir/br1_tr1_pair1.unmatched.bedpe \
#    1018_out_dir/br1_tr1_pair1.matching.html \
#    > 1018_out_dir/br1_tr1_pair1.bedpe.mapq30

# 对BEDPE文件进行排序 
sort -k12,13 -k1,1 -k2,2n -k3,3 -k4,4n 1018_out_dir/br1_tr1_pair1.bedpe.mapq0 > 1018_out_dir/br1_tr1_pair1.bedpe.mapq0.sorted
sort -k7 1018_out_dir/br1_tr1_pair1.bedpe.mapq0 > 1018_out_dir/br1_tr1_pair1.bedpe.mapq0.sorted

# 去重
python dedupe_scDHC.py \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq0.sorted \
    > 1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped

#10
# 过滤未定位的contig和卫星序列
python filter_bedpe.py \
    combo_hg38.genomesize \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped \
    > 1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped.filtered

#11
# 计算总的barcode reads数量
wc -l ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bc_clipped > 1018_out_dir/br1_tr1_pair1.ph1
read reads_r1 junk < 1018_out_dir/br1_tr1_pair1.ph1
total_reads=$((reads_r1 / 4))
print('total_reads')
#rm 1018_out_dir/br1_tr1_pair1.ph1

# 统计映射的reads数量
wc -l 1018_out_dir/br1_tr1_pair1.bedpe.mapq0 > 1018_out_dir/br1_tr1_pair1.ph1
read mapped_reads_mapq junk < 1018_out_dir/br1_tr1_pair1.ph1
#rm 1018_out_dir/br1_tr1_pair1.ph1

# 统计去重后的reads数量
wc -l 1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped > 1018_out_dir/br1_tr1_pair1.ph1
read associated_reads_mapq_dedupe junk < 1018_out_dir/br1_tr1_pair1.ph1
rm 1018_out_dir/br1_tr1_pair1.ph1

# 生成HTML报告
cat 1018_out_dir/splitting_stats.html >> 1018_out_dir/br1_tr1_pair1.baseline_stats.html
echo "<H3>Total Reads With Barcode Found: $total_reads</H3>" >> 1018_out_dir/br1_tr1_pair1.baseline_stats.html
cat 1018_out_dir/br1_tr1_pair1.matching.html >> 1018_out_dir/br1_tr1_pair1.baseline_stats.html
echo "<H3>Total Mapped Reads MAPQ > 30: $mapped_reads_mapq</H3>" >> 1018_out_dir/br1_tr1_pair1.baseline_stats.html
echo "<H3>Total Deduped Mapped Reads MAPQ > 30: $associated_reads_mapq_dedupe</H3>" >> 1018_out_dir/br1_tr1_pair1.baseline_stats.html

#漏执行部分
#python sort_strandedness.py $bc_assoc.bedpe.mapq0 $bc_assoc.deduped.stats.percent_plots > $bc_assoc.deduped.stats.html&
#python sort_strandedness.py $bc_assoc.bedpe.mapq0.deduped.filtered $bc_assoc.associated.stats.percent_plots > $bc_assoc.associated.stats.html&



#12
# 分析物种特异性条形码
python calculate_cell_distro_long.py \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq0 \
    > 1018_out_dir/br1_tr1_pair1.percentages \
    2> 1018_out_dir/br1_tr1_pair1.REoccurrences

python calculate_cell_distro_long.py \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped.filtered \
    > 1018_out_dir/br1_tr1_pair1.deduped.percentages \
    2> 1018_out_dir/br1_tr1_pair1.deduped.REoccurrences  

#13
# 计算单细胞片段长度分布 
python single_cell_length_distros.py \
    1018_out_dir/br1_tr1_pair1.deduped.percentages \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped.filtered \
    br1_tr1_pair1 \
    > 1018_out_dir/br1_tr1_pair1.single_cell_lengths

# 计算Cis-Trans比率
python calculate_median_cistrans_ratio.py \
    1018_out_dir/br1_tr1_pair1.single_cell_lengths \
    > 1018_out_dir/br1_tr1_pair1.cistrans.txt

# 将Cis-Trans比率整合到百分比文件中
python incorporate_cistrans_genotype.py \
    1018_out_dir/br1_tr1_pair1.cistrans.txt \
    1018_out_dir/br1_tr1_pair1.deduped.percentages \
    > 1018_out_dir/br1_tr1_pair1.deduped.percentages.filterable

# (1111处理到这一步)
#14
# 生成单细胞接触矩阵
python bin_schic.py \
    combo_hg38.genomesize \
    1018_out_dir/br1_tr1_pair1.deduped.percentages.filterable \
    1018_out_dir/br1_tr1_pair1.bedpe.mapq30.deduped.filtered > foo

# 将矩阵文件移动到单独的目录
mkdir 1018_out_dir/br1_tr1_pair1.matrices
mv ./*.matrix 1018_out_dir/br1_tr1_pair1.matrices

# 生成有效的细胞矩阵列表
cd 1018_out_dir/br1_tr1_pair1.matrices
ls -1 *.matrix > valid_cells.list

#(optional)
# 删除临时目录
#rm -rf $tmp_dir
