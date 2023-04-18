#! /bin/bash

bowtie2 -x $1 -1 $2 -2 $3 -S $4  2>$5
samtools fastq -@ 8 -f 4 $4 -1 $6 -2 $7 -s $8

## demo
# 去宿主
# bowtie2 -x /root/database/kneaddata_human/hg37dec_v0.1 -1 DBYX1FB_clean.1.fastq.gz -2 DBYX1FB_clean.2.fastq.gz -S DBYX1FB.sam  2>DBYX1FB.mapping.log
# 提取未比对到宿主的序列, -f 4表示：过滤掉未比对到参考基因组的reads（Flag=4）
# samtools fastq -@ 8 -f 4 DBYX1FB.sam -1 DBYX1FB.unmap.1.fastq.gz -2 DBYX1FB.unmap.2.fastq.gz -s DBYX1FB.unmap.single.fastq.gz
