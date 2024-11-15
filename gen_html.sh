# 计算总的barcode reads数量
wc -l ./1018_tmp/4DNFI7HX4GMM.fastq.gz.bc_clipped > 1018_out_dir/br1_tr1_pair1.ph1
read reads_r1 junk < 1018_out_dir/br1_tr1_pair1.ph1
total_reads=$((reads_r1 / 4))
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