#! /bin/bash

## 质控
fastp -i DBYX1FB_1.fq.gz -o DBYX1FB_clean.1.fastq.gz -I DBYX1FB_2.fq.gz -O DBYX1FB_clean.2.fastq.gz -w 8 -h DBYX1FB.html -j DBYX1FB.json

## 去宿主
bowtie2 -x /root/database/kneaddata_human/hg37dec_v0.1 -1 DBYX1FB_clean.1.fastq.gz -2 DBYX1FB_clean.2.fastq.gz -S DBYX1FB.sam  2>DBYX1FB.mapping.log
# 提取未比对到宿主的序列, -f 4表示：过滤掉未比对到参考基因组的reads（Flag=4）
samtools fastq -@ 8 -f 4 DBYX1FB.sam -1 DBYX1FB.unmap.1.fastq.gz -2 DBYX1FB.unmap.2.fastq.gz -s DBYX1FB.unmap.single.fastq.gz

## 基于metaphlan物种注释
zcat DBYX1FB.unmap.1.fastq.gz DBYX1FB.unmap.2.fastq.gz|metaphlan --input_type fastq --bowtie2out DBYX1FB_bowtie2.bz2 --output_file DBYX1FB_metaphlan.tsv --nproc 8
# 单个样本合并物种信息
merge_metaphlan_tables.py metaphlan.tsv > merged_abundance_table.txt
# 多个样本合并物种信息
merge_metaphlan_tables.py *.tsv > merged_abundance_table.txt
# 按照界门纲目科属种来生成物种组成表
grep -E '(p__)|(clade_name)' merged_abundance_table.txt |grep -v 'c__'|sed 's/|/;/g' > merged_abundance_table_phylum.txt
grep -E '(c__)|(clade_name)' merged_abundance_table.txt |grep -v 'o__'|sed 's/|/;/g' > merged_abundance_table_class.txt
grep -E '(o__)|(clade_name)' merged_abundance_table.txt |grep -v 'f__'|sed 's/|/;/g' > merged_abundance_table_order.txt
grep -E '(f__)|(clade_name)' merged_abundance_table.txt |grep -v 'g__'|sed 's/|/;/g' > merged_abundance_table_family.txt
grep -E '(g__)|(clade_name)' merged_abundance_table.txt |grep -v 's__'|sed 's/|/;/g' > merged_abundance_table_genus.txt
grep -E '(s__)|(clade_name)' merged_abundance_table.txt |grep -v 't__'|sed 's/|/;/g' > merged_abundance_table_species.txt

## megahit组装
megahit -1 DBYX1FB.unmap.1.fastq.gz -2 DBYX1FB.unmap.2.fastq.gz -o DBYX1FB_megahit --out-prefix DBYX1FB -t 8
# 过滤500bp以下的组装序列,修改序列名称
seqkit seq -m 500 DBYX1FB_megahit/DBYX1FB.contigs.fa --remove-gaps > DBYX1FB.contigs_500.fa
sed -i 's/>/>DBYX1FB_/g' DBYX1FB.contigs_500.fa

## prodigal基因预测
prodigal -p meta -a DBYX1FB_prot.faa -m -d DBYX1FB_nucl.fna -o DBYX1FB_genes.gff -f gff -s DBYX1FB.stat -i DBYX1FB.contigs_500.fa

## 假设这里已经是所有样本都处理完了,分别合并了所有样本的fna和faa
mkdir total
cp DBYX1FB_prot.faa total/prot.faa
cp DBYX1FB_nucl.fna total/nucl.fna

## cd-hit去冗余,先用蛋白序列去冗余
# 由于不同的核酸序列翻译后可能产生相同的蛋白质序列，因此对核酸序列的去冗余需要基于蛋白质序列
# 基于蛋白质序列文件与核酸序列文件中序列号的对应关系，根据去冗余后的蛋白质序列文件的序列号筛选核酸序列
cd total
cd-hit -i prot.faa -o prot_nonerude.faa -c 0.95 -T 8 -n 5 -d 0 -aS 0.9 -g 1 -sc 1 -sf 1 -M 0
less prot_nonerude.faa|grep '>'|awk -F ' ' '{print $1}'|sed 's/>//g' > prot_nonerude.list
seqtk subseq nucl.fna prot_nonerude.list > nucl_nonerude.fna
# 非冗余基因集构建bwa索引
bwa index nucl_nonerude.fna -p geneset_bwa
# 计算非冗余基因集中各个基因的长度，为了后面可能计算基因的RPKM等信息，需要基因的长度信息
(echo -e "gene\tlength"; bioawk -c fastx '{print $name, length($seq)}' nucl_nonerude.fna) > geneset_length.txt

## 非冗余基因集prot_nonerude.faa利用eggnog-mapper功能注释, --cpu 0表示用所有cpu, --usemem表示断点续存
emapper.py -i prot_nonerude.faa -o eggnog --cpu 0 --usemem
# 生成非冗余基因集的KO列表
less -S eggnog.emapper.annotations|cut -f1,12|grep -v "^#"|sed 's/ko://g'|sed '1i gene\tko'|awk '$2 !~ /-/ {print}' > KEGG_KO.txt
# 生成非冗余基因集的PATHWAY列表
less -S eggnog.emapper.annotations|cut -f1,13|grep -v "^#"|sed '1i gene\tpathway'|awk '$2 !~ /-/ {print}' > KEGG_PATHWAY.txt

## 样本基因丰度计算
cd ../
bwa mem -t 4 total/geneset_bwa DBYX1FB.unmap.1.fastq.gz DBYX1FB.unmap.2.fastq.gz | samtools view -bS - | samtools sort - > DBYX1FB_mapping_geneset.bam
# 过滤掉未比对到基因集的reads(Flag=4),过滤掉二次比对的reads(Flag=256),过滤掉嵌合reads(Flag=2048),过滤掉比对数小于2的
samtools view -F 4 -F 256 -F 2048 DBYX1FB_mapping_geneset.bam|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\tDBYX1FB' > DBYX1FB.count
# 不过滤比对数目小于2的
samtools view -F 4 -F 256 -F 2048 DBYX1FB_mapping_geneset.bam|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {print $1"\t"$2}'|sed '1i gene\tDBYX1FB' > DBYX1FB.count

## 合并多个样本的基因丰度表
# DBYX1FB.count.test1、DBYX1FB.count.test2是测试有多个样本的基因丰度表，和DBYX1FB.count一起合并
# 最终生成的merged_file.txt就是原始(prodigal预测后的基因)的样本基因丰度表，后续比对各功能数据库时，直接用id信息匹配获取对应功能数据库的基因丰度表
# R实现
'''R
files <- c("DBYX1FB.count", "DBYX1FB.count.test1", "DBYX1FB.count.test2")
df_list <- lapply(files, function(x) read.table(x, sep = '\t', header = TRUE))
df_merged <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), df_list)
df_merged[is.na(df_merged)] <- 0
write.table(df_merged, "merged_file.txt", sep = '\t', row.names = FALSE, quote = F)
'''
# python实现
'''python
import pandas as pd
files = ['DBYX1FB.count', 'DBYX1FB.count.test1', 'DBYX1FB.count.test2']
df_list = [pd.read_csv(f, sep='\t') for f in files]
df_merged = pd.concat(df_list).groupby('gene', as_index=False).sum()
df_merged.fillna(0, inplace=True)
df_merged.to_csv('merged_file.txt', sep='\t', index=False)
'''
