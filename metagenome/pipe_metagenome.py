#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   pipe_metagenome.py
@Time    :   2023/04/14 10:39:41
@Author  :   Wanjin.Hu 
@Version :   1.0
@Contact :   wanjin.hu@diprobio.com
@Description : pipeline of metagenome
'''

import argparse
import os
import pandas as pd

class PipeMetagenome(object):
    """
    @功能: 宏基因组分析流程
    """

    def __init__(self):
        super(PipeMetagenome).__init__()
        self.bash_env = "/usr/bin/bash"
        self.py_env = "/root/miniconda3/bin/python"
        self.fastp_cmd = "/root/pipe_script/metagenome/pipe_script/fastp.sh"
        self.ref_remove_cmd = "/root/pipe_script/metagenome/pipe_script/ref_remove.sh"
        self.metaphlan_cmd = "/root/pipe_script/metagenome/pipe_script/metaphlan.sh"
        self.metaphlan_merge_cmd = "/root/pipe_script/metagenome/pipe_script/metaphlan_merge.sh"
        self.metaphlan_merge = "/root/miniconda3/bin/merge_metaphlan_tables.py"
        self.megahit_cmd = "/root/pipe_script/metagenome/pipe_script/megahit.sh"
        self.prodigal_cmd = "/root/pipe_script/metagenome/pipe_script/prodigal.sh"
        self.cdhit_cmd = "/root/pipe_script/metagenome/pipe_script/cd-hit.sh"
        self.emapper_cmd = "/root/pipe_script/metagenome/pipe_script/emapper.sh"
        self.count_cmd = "/root/pipe_script/metagenome/pipe_script/gene_count.sh"
        self.kegg_cmd = "/root/pipe_script/metagenome/pipe_script/kegg/kegg.py"
        if os.path.exists(args.outDir):
            print("{}输出文件夹已存在,进行覆盖操作!".format(args.outDir))
            self.log = open((args.outDir + "/run.log"),mode="a",encoding="utf-8")
            # quit()
        if not os.path.exists(args.outDir):
            os.mkdir(args.outDir)
            self.log = open((args.outDir + "/run.log"),mode="a",encoding="utf-8")
        self.fastp_path = os.path.join(args.outDir,"01-fastp_trim")
        if os.path.exists(self.fastp_path):
            pass
        else:
            os.makedirs(self.fastp_path, 0o755)
        self.ref_remove_path = os.path.join(args.outDir,"02-ref_remove")
        if os.path.exists(self.ref_remove_path):
            pass
        else:
            os.makedirs(self.ref_remove_path, 0o755)
        self.metaphlan_path = os.path.join(args.outDir,"03-metaphlan")
        if os.path.exists(self.metaphlan_path):
            pass
        else:
            os.makedirs(self.metaphlan_path, 0o755)
        self.megahit_path = os.path.join(args.outDir,"04-megahit")
        if os.path.exists(self.megahit_path):
            pass
        else:
            os.makedirs(self.megahit_path, 0o755)
        self.prodigal_path = os.path.join(args.outDir,"05-prodigal")
        if os.path.exists(self.prodigal_path):
            pass
        else:
            os.makedirs(self.prodigal_path, 0o755)
        self.cdhit_path = os.path.join(args.outDir,"06-cdhit")
        if os.path.exists(self.cdhit_path):
            pass
        else:
            os.makedirs(self.cdhit_path, 0o755)
        self.emapper_path = os.path.join(args.outDir,"07-emapper")
        if os.path.exists(self.emapper_path):
            pass
        else:
            os.makedirs(self.emapper_path, 0o755)
        self.samcount_path = os.path.join(args.outDir,"08-sam_count")
        if os.path.exists(self.samcount_path):
            pass
        else:
            os.makedirs(self.samcount_path, 0o755)
        self.kegg_path = os.path.join(args.outDir,"09-emapper_kegg")
        if os.path.exists(self.kegg_path):
            pass
        else:
            os.makedirs(self.kegg_path, 0o755)
        self.raw_dict = {}
        with open(args.fqList) as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    line = line.strip().split("\t")
                    self.raw_dict[line[0]] = line[1:]
        
    def fastp_trim(self):
        """
        @功能: fastp质控
        """
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_r1 = v[0]
            sam_r2 = v[1]
            sam_fastp_path = os.path.join(self.fastp_path,sam_name)
            if os.path.exists(sam_fastp_path):
                continue
            else:
                os.makedirs(sam_fastp_path, 0o755)
                sam_r1_clean = sam_fastp_path + "/{}_clean.1.fastq.gz".format(sam_name)
                sam_r2_clean = sam_fastp_path + "/{}_clean.2.fastq.gz".format(sam_name)
                sam_html = sam_fastp_path + "/{}.html".format(sam_name)
                sam_json = sam_fastp_path + "/{}.json".format(sam_name)
                cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.fastp_cmd,
                                                    sam_r1,
                                                    sam_r1_clean,
                                                    sam_r2,
                                                    sam_r2_clean,
                                                    sam_html,
                                                    sam_json)
                os.system(command=cmd)
                print(cmd,file=self.log)
                print("Fastp finished: {}, check result dir: {}".format(sam_name,sam_fastp_path))
        print("All fastp finished, check dir: {}".format(self.fastp_path))

    def ref_remove(self):
        """
        @功能: 去宿主
        """
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_ref_remove_path = os.path.join(self.ref_remove_path,sam_name)
            if os.path.exists(sam_ref_remove_path):
                continue
            else:
                os.makedirs(sam_ref_remove_path,0o755)
                sam_r1_clean = self.fastp_path + "/{}/{}_clean.1.fastq.gz".format(sam_name,sam_name)
                sam_r2_clean = self.fastp_path + "/{}/{}_clean.2.fastq.gz".format(sam_name,sam_name)
                sam_ref_sam = sam_ref_remove_path + "/{}.sam".format(sam_name)
                sam_ref_sam_log = sam_ref_remove_path + "/{}.mapping.log".format(sam_name)
                sam_r1_unmap = sam_ref_remove_path + "/{}.unmap.1.fastq.gz".format(sam_name) 
                sam_r2_unmap = sam_ref_remove_path + "/{}.unmap.2.fastq.gz".format(sam_name)
                sam_unmap_single = sam_ref_remove_path + "/{}.unmap.single.fastq.gz".format(sam_name)
                cmd = "{} {} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                             self.ref_remove_cmd,
                                                             args.ref,
                                                             sam_r1_clean,
                                                             sam_r2_clean,
                                                             sam_ref_sam,
                                                             sam_ref_sam_log,
                                                             sam_r1_unmap,
                                                             sam_r2_unmap,
                                                             sam_unmap_single)
                os.system(command=cmd)
                print(cmd,file=self.log)
                print("Ref sequence remove finished: {}, check result dir: {}".format(sam_name,sam_ref_remove_path))
        print("All ref sequence remove finished, check dir: {}".format(self.ref_remove_path))

    def metaphlan(self):
        """
        @功能: metaphlan物种注释
        """
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_metaphlan_path = os.path.join(self.metaphlan_path,sam_name)
            if os.path.exists(sam_metaphlan_path):
                continue
            else:
                os.makedirs(sam_metaphlan_path,0o755)
                sam_r1_unmap = self.ref_remove_path + "/{}/{}.unmap.1.fastq.gz".format(sam_name,sam_name) 
                sam_r2_unmap = self.ref_remove_path + "/{}/{}.unmap.2.fastq.gz".format(sam_name,sam_name)
                sam_metaphlan_bz2 = sam_metaphlan_path + "/{}_bowtie2.bz2".format(sam_name)
                sam_metaphlan_out = sam_metaphlan_path + "/{}.tsv".format(sam_name) # 输出文件名为sample.tsv, 后续合并表格就是样本名
                cmd = "{} {} {} {} {} {}".format(self.bash_env,
                                                 self.metaphlan_cmd,
                                                 sam_r1_unmap,
                                                 sam_r2_unmap,
                                                 sam_metaphlan_bz2,
                                                 sam_metaphlan_out)
                os.system(command=cmd)
                print(cmd,file=self.log)
                print("Metaphlan finished: {}, check result dir: {}".format(sam_name,sam_metaphlan_path))
        metaphlan_merge = self.metaphlan_path + "/00_merged_abundance_table.txt"
        metaphlan_phylum = self.metaphlan_path + "/01_metaphlan_phylum.txt"
        metaphlan_class = self.metaphlan_path + "/02_metaphlan_class.txt"
        metaphlan_order = self.metaphlan_path + "/03_metaphlan_order.txt"
        metaphlan_family = self.metaphlan_path + "/04_metaphlan_family.txt"
        metaphlan_genus = self.metaphlan_path + "/05_metaphlan_genus.txt"
        metaphlan_species = self.metaphlan_path + "/06_metaphlan_species.txt"
        if os.path.exists(metaphlan_merge):
            pass
        else:  
            cmd1 = "{} {}/*/*.tsv > {}".format(self.metaphlan_merge,self.metaphlan_path,metaphlan_merge)
            os.system(command=cmd1)
            print(cmd1,file=self.log)
            cmd2 = "{} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.metaphlan_merge_cmd,
                                                    metaphlan_merge,
                                                    metaphlan_phylum,
                                                    metaphlan_class,
                                                    metaphlan_order,
                                                    metaphlan_family,
                                                    metaphlan_genus,
                                                    metaphlan_species)
            os.system(command=cmd2)
            print(cmd2,file=self.log)
        print("Metaphlan merge finished, check dir: {}".format(self.metaphlan_path))
    
    def megahit(self):
        """
        @功能: 组装
        """
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_megahit_path = os.path.join(self.megahit_path,sam_name)
            if os.path.exists(sam_megahit_path):
                continue
            else:
                # os.makedirs(sam_megahit_path,0o755) # megahit运行输出文件夹不能为已有的
                sam_r1_unmap = self.ref_remove_path + "/{}/{}.unmap.1.fastq.gz".format(sam_name,sam_name) 
                sam_r2_unmap = self.ref_remove_path + "/{}/{}.unmap.2.fastq.gz".format(sam_name,sam_name)
                sam_contig = sam_megahit_path + "/{}.contigs.fa".format(sam_name)
                sam_contig_500 = sam_megahit_path + "/{}.contigs_500.fa".format(sam_name)
                cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
                                                       self.megahit_cmd,
                                                       sam_r1_unmap,
                                                       sam_r2_unmap,
                                                       sam_megahit_path,
                                                       sam_name,
                                                       sam_contig,
                                                       sam_contig_500)
                os.system(command=cmd)
                print(cmd,file=self.log)
                print("Megahit finished: {}, check result dir: {}".format(sam_name,sam_megahit_path))
        print("All megahit finished, check dir: {}".format(self.megahit_path))
    
    def prodigal(self):
        """
        @功能: 基因预测
        """
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_prodigal_path = os.path.join(self.prodigal_path,sam_name)
            if os.path.exists(sam_prodigal_path):
                continue
            else:
                os.makedirs(sam_prodigal_path,0o755)
                sam_prot = sam_prodigal_path + "/{}_prot.faa".format(sam_name)
                sam_nucl = sam_prodigal_path + "/{}_nucl.fna".format(sam_name)
                sam_gff = sam_prodigal_path + "/{}_genes.gff".format(sam_name)
                sam_stat = sam_prodigal_path + "/{}.stat".format(sam_name)
                sam_contig_500 = self.megahit_path + "/{}/{}.contigs_500.fa".format(sam_name,sam_name)
                cmd = "{} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.prodigal_cmd,
                                                    sam_prot,
                                                    sam_nucl,
                                                    sam_gff,
                                                    sam_stat,
                                                    sam_contig_500)
                os.system(command=cmd)
                print(cmd,file=self.log)
                print("Prodigal finished: {}, check result dir: {}".format(sam_name,sam_prodigal_path))
        print("All prodigal finished, check dir: {}".format(self.prodigal_path))
    
    def cdhit(self):
        """
        @功能: cdhit去冗余,非冗余基因集构建bwa索引
        """
        # 合并各个样本的基因集
        prot_cat = self.cdhit_path + "/prot.faa"
        nucl_cat = self.cdhit_path + "/nucl.fna"
        prot_nonerude = self.cdhit_path + "/prot_nonerude.faa"
        nucl_nonerude = self.cdhit_path + "/nucl_nonerude.fna"
        nonerude_list = self.cdhit_path + "/prot_nonerude.list"
        geneset_bwa = self.cdhit_path + "/geneset_bwa"
        geneset_length = self.cdhit_path + "/geneset_length.txt"
        if os.path.exists(geneset_length):
            print("CD-HIT 结果已经存在,本次不重新运行! check dir: {}".format(self.cdhit_path))
            pass
        else:
            cmd1 = "cat {}/*/*_prot.faa > {}".format(self.prodigal_path,prot_cat)
            os.system(command=cmd1)
            print(cmd1,file=self.log)
            cmd2 = "cat {}/*/*_nucl.fna > {}".format(self.prodigal_path,nucl_cat)
            os.system(command=cmd2)
            print(cmd2,file=self.log)
            cmd = "{} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                    self.cdhit_cmd,
                                                    prot_cat,
                                                    prot_nonerude,
                                                    nonerude_list,
                                                    nucl_cat,
                                                    nucl_nonerude,
                                                    geneset_bwa,
                                                    geneset_length)
            os.system(command=cmd)
            print(cmd,file=self.log)
            print("CD-HIT finished, check dir: {}".format(self.cdhit_path))
    
    def emapper(self):
        """
        @功能: emapper功能注释, 生成KO/KEGG信息表
        """
        prot_nonerude = self.cdhit_path + "/prot_nonerude.faa"
        eggnog_profile = self.emapper_path + "/eggnog"
        eggnog_ann = self.emapper_path + "/eggnog.emapper.annotations"
        eggnog_KO = self.emapper_path + "/KEGG_KO.txt"
        eggnog_path = self.emapper_path + "/KEGG_PATHWAY.txt"
        if os.path.exists(eggnog_path):
            print("eggnog-mapper 结果已经存在,本次不重新运行! check dir: {}".format(self.emapper_path))
            pass
        else:
            cmd = "{} {} {} {} {} {} {}".format(self.bash_env,
                                                self.emapper_cmd,
                                                prot_nonerude,
                                                eggnog_profile,
                                                eggnog_ann,
                                                eggnog_KO,
                                                eggnog_path)
            os.system(command=cmd)
            print(cmd,file=self.log)
            print("eggnog-mapper finished, check dir: {}".format(self.emapper_path))
    
    def sam_gene_count(self):
        """
        @功能: 各样本丰度计算；以及合并所有样本
        """
        count_list = []
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_count_path = os.path.join(self.samcount_path,sam_name)
            if os.path.exists(sam_count_path):
                count_file = sam_count_path + "/{}.count".format(sam_name)
                count_list.append(count_file)
                continue
            else:
                os.makedirs(sam_count_path,0o755)
                geneset_bwa = self.cdhit_path + "/geneset_bwa"
                sam_r1_unmap = self.ref_remove_path + "/{}/{}.unmap.1.fastq.gz".format(sam_name,sam_name)
                sam_r2_unmap = self.ref_remove_path + "/{}/{}.unmap.2.fastq.gz".format(sam_name,sam_name)
                sam_count_bam = sam_count_path + "/{}_mapping_geneset.bam".format(sam_name)
                sam_count_txt = sam_count_path + "/{}.count".format(sam_name)
                cmd = "{} {} {} {} {} {} {} {}".format(self.bash_env,
                                                       self.count_cmd,
                                                       geneset_bwa,
                                                       sam_r1_unmap,
                                                       sam_r2_unmap,
                                                       sam_count_bam,
                                                       sam_name,
                                                       sam_count_txt)
                os.system(command=cmd)
                print(cmd,file=self.log)
                count_list.append(sam_count_txt)
                print("sample gene count finished, check dir: {}".format(self.samcount_path))
        # 合并样本
        # files = ['DBYX1FB.count', 'DBYX1FB.count.test1', 'DBYX1FB.count.test2']
        if os.path.exists(self.samcount_path + "/merged_file.txt"):
            pass
        else:
            df_list = [pd.read_csv(f, sep='\t') for f in count_list]
            df_merged = pd.concat(df_list).groupby('gene', as_index=False).sum()
            df_merged.fillna(0, inplace=True)
            merge_file = self.samcount_path + "/merged_file.txt"
            df_merged.to_csv(merge_file, sep='\t', index=False)
            print("样本基因丰度文件合并完毕, check file: {}".format(merge_file))
    
    def eggnog_kegg(self):
        """
        @功能: 利用eggnog-mapper方式得到的kegg信息获取KO/pathway结果表
        """
        eggnog_KO = self.emapper_path + "/KEGG_KO.txt"
        eggnog_path = self.emapper_path + "/KEGG_PATHWAY.txt"
        merge_file = self.samcount_path + "/merged_file.txt"
        out_ko = self.kegg_path + "/KO_samples.csv"
        out_pathway = self.kegg_path + "/pathway_samples.csv"
        if os.path.exists(out_pathway):
            pass
        else:
            cmd = "{} {} -kk {} -kp {} -mt {} -ok {} -op {}".format(self.py_env,
                                                                    self.kegg_cmd,
                                                                    eggnog_KO,
                                                                    eggnog_path,
                                                                    merge_file,
                                                                    out_ko,
                                                                    out_pathway)
            os.system(command=cmd)
            print(cmd,file=self.log)
            print("eggnog-kegg finished, check dir: {}".format(self.kegg_path))
        
    def main(self):
        self.fastp_trim()
        self.ref_remove()
        self.metaphlan()
        self.megahit()
        self.prodigal()
        self.cdhit()
        self.emapper()
        self.sam_gene_count()
        self.eggnog_kegg()
        self.log.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="PipeMetagenome",
                                     usage="python pipe_metagenome.py --fastq_list fq.list --output_dir result --ref /root/database/Canis_GCF_000002285.5/Canis_GCF_000002285_5",
                                     description="Pipeline of metagenome")
    # fastq_list必须保留表头格式固定为: #Sample R1  R2
    parser.add_argument("-l", "--fastq_list", dest="fqList", required=True, type=str, help="raw fq list")
    parser.add_argument("-o", "--output_dir", dest="outDir", required=True, type=str, help="result output")
    parser.add_argument("-r", "--ref", dest="ref", required=True, type=str, help="ref genome bowtie2 index")
    args = parser.parse_args()
    run = PipeMetagenome()
    run.main()
    # run.fastp_trim()
    # run.ref_remove()
    # run.metaphlan()
    # run.megahit()
    # run.prodigal()
    # run.cdhit()
    # run.emapper()
    # run.sam_gene_count()
    # run.eggnog_kegg()
