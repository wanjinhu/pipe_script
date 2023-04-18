# ------------------------ 
# kneaddata质控
# ------------------------
kneaddata -i1 DBYX1FB_FDSW210039579-1r_1.fq.gz -i2 DBYX1FB_FDSW210039579-1r_2.fq.gz -o test --cat-final-output --output-prefix test --serial -v -t 24 -db /root/database/kneaddata_human --run-fastqc-start --run-fastqc-end --remove-intermediate-output
kneaddata_read_count_table --input test --output kneaddata_read_count_table.tsv

kneaddata --input1 DBYX1FB_FDSW210039579-1r_1.fq.gz --input2 DBYX1FB_FDSW210039579-1r_2.fq.gz -o DBYX1FB --output-prefix DBYX1FB -t 16 -db /root/database/kneaddata_human/hg37dec_v0.1 --remove-intermediate-output --bowtie2-options="--very-sensitive" --trimmomatic-options="MINLEN:50 ILLUMINACLIP:/usr/lib/python2.7/site-packages/kneaddata/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20 MINLEN:50"

# ------------------------
# metaphlan物种注释
# ------------------------
conda activate metaphlan4
metaphlan test/test.fastq --input_type fastq --bowtie2out metagenome.bowtie2.bz2 --output_file metaphlan.tsv --nproc 24 --bowtie2db /root/database/metaphlan_databases --index mpa_vOct22_CHOCOPhlAnSGB_202212
# 单个样本合并物种信息
merge_metaphlan_tables.py metaphlan.tsv > merged_abundance_table.txt
# 多个样本合并物种信息
merge_metaphlan_tables.py *.txt > merged_abundance_table.txt
# 按照界门纲目科属种来生成物种组成表
grep -E '(p__)|(clade_name)' merged_abundance_table.txt |grep -v 'c__'|sed 's/|/;/g' > merged_abundance_table_phylum.txt
grep -E '(c__)|(clade_name)' merged_abundance_table.txt |grep -v 'o__'|sed 's/|/;/g' > merged_abundance_table_class.txt
grep -E '(o__)|(clade_name)' merged_abundance_table.txt |grep -v 'f__'|sed 's/|/;/g' > merged_abundance_table_order.txt
grep -E '(f__)|(clade_name)' merged_abundance_table.txt |grep -v 'g__'|sed 's/|/;/g' > merged_abundance_table_family.txt
grep -E '(g__)|(clade_name)' merged_abundance_table.txt |grep -v 's__'|sed 's/|/;/g' > merged_abundance_table_genus.txt
grep -E '(s__)|(clade_name)' merged_abundance_table.txt |grep -v 't__'|sed 's/|/;/g' > merged_abundance_table_species.txt

# 其他可选处理
# onvert a SGB-based MetaPhlAn 4 output into a GTDB-taxonomy-based profile
sgb_to_gtdb_profile.py -i metaphlan.tsv -o metaphlan_output_gtdb.txt -d /root/database/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl

# -------------------------
# 功能注释(利用拼接的方式)
# -------------------------

