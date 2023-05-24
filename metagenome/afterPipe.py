#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   afterPipe.py
@Time    :   2023/05/22 18:23:57
@Author  :   Wanjin.Hu 
@Version :   1.0
@Contact :   wanjin.hu@diprobio.com
@Description :  pipe_metagenome运行完之后的操作
'''

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(prog="AfterMetagenomePipe",
                                 usage="",
                                 description="After pipeline of metagenome")
parser.add_argument("-i", "--input_dir", dest="inDir", required=True, type=str, help="input dir")
parser.add_argument("-o", "--output_dir", dest="outDir", required=True, type=str, help="output dir")
parser.add_argument("-g", "--group", dest="group", required=True, type=str, help="which group you use")
parser.add_argument("-l", "--lda", dest="lda", required=True, type=str, help="set lda score")

args = parser.parse_args()

class AfterPipe(object):
    """
    @功能: pipe_metagenome运行完之后的操作
    """

    def __init__(self):
        super(AfterPipe).__init__()
        if not os.path.exists(args.outDir):
            os.mkdir(args.outDir)
        # 程序相关
        self.py_env = "/usr/bin/python"
        self.py_lefse = "/root/pipe_script/general/lefse/lefse.py"
        # input_dir包含以下固定文件
        self.group = args.inDir + "/group.txt"
        self.sp_raw = args.inDir + "/00_merged_abundance_table.txt"
        self.phylum = args.inDir + "/01_metaphlan_phylum.txt"
        self.class_ = args.inDir + "/02_metaphlan_class.txt"
        self.order = args.inDir + "/03_metaphlan_order.txt"
        self.family = args.inDir + "/04_metaphlan_family.txt"
        self.genus = args.inDir + "/05_metaphlan_genus.txt"
        self.species = args.inDir + "/06_metaphlan_species.txt"
        self.pathway = args.inDir + "/pathway_samples.xls"
        self.ko = args.inDir + "/KO_samples.xls"
        
    def lefse_format(self,df,group):
        """
        @功能: 
        @参数df: factor组成表[pd读取的]: 行名×列名[factor × sample]
        @参数group: 分组表[pd读取的]
        """
        df = df.rename(columns={df.columns[0]: 'sample'})
        df_t = df.set_index('sample').T
        merged_df = pd.merge(group, df_t, left_index=True, right_index=True)
        merged_df.insert(1, column='sample', value=merged_df.index)
        return merged_df.T
        
    def pre_format(self,g_name):
        """
        @功能: 分析前的文件格式整理,包括:
        @参数g_name: 确定哪个分组文件[字符串]
        @fun_1: pathway和KO绝对丰度表转相对丰度表
        @fun_2: pathway的level1和2绝对丰度和相对丰度表整理
        @fun_3: pathway和KO表整理成lefse分析的格式表[用相对丰度表]
        @fun_4: 物种数据整理成lefse分析的格式表
        """
        ## fun_1
        ko_r = args.outDir + "/KO_r.txt"
        ko_a = args.outDir + "/KO_a.txt"
        df_ko = pd.read_csv(self.ko, sep='\t', index_col=2)
        df_ko = df_ko.drop(df_ko.columns[[0,1]], axis=1)
        df_ko.to_csv(ko_a,sep='\t')
        sums_ko = df_ko.sum(axis=0)
        df_ko_relative = df_ko.div(sums_ko)
        df_ko_relative.to_csv(ko_r,sep='\t')
        pathway_r = args.outDir + "/pathway_r.txt"
        pathway_a = args.outDir + "/pathway_a.txt"
        df_pathway = pd.read_csv(self.pathway, sep='\t', index_col=3)
        df_pathway = df_pathway.drop(df_pathway.columns[[0,1,2]], axis=1)
        df_pathway.to_csv(pathway_a,sep='\t')
        sums_pathway = df_pathway.sum(axis=0)
        df_pathway_relative = df_pathway.div(sums_pathway)
        df_pathway_relative.to_csv(pathway_r,sep='\t')
        ## fun_2
        lev_2_a = args.outDir + "/pathway_lev2_a.txt"
        lev_2_r = args.outDir + "/pathway_lev2_r.txt"
        df_1 = pd.read_csv(self.pathway, sep='\t')
        lev_2 = df_1.groupby('level2').sum()
        sums_lev_2 = lev_2.sum(axis=0)
        lev_2_relative = lev_2.div(sums_lev_2)
        n_rows = lev_2.shape[0]
        new_cols = [f'level2_{i}' for i in range(1, n_rows+1)]
        lev_2.insert(0, column='name', value=new_cols)
        lev_2_relative.insert(0, column='name', value=new_cols)
        lev_2.to_csv(lev_2_a,sep='\t')
        lev_2_relative.to_csv(lev_2_r,sep='\t')
        lev_1_a = args.outDir + "/pathway_lev1_a.txt"
        lev_1_r = args.outDir + "/pathway_lev1_r.txt"
        lev_1 = df_1.groupby('level1').sum()
        sums_lev_1 = lev_1.sum(axis=0)
        lev_1_relative = lev_1.div(sums_lev_1)
        n_rows_1 = lev_1.shape[0]
        new_cols_1 = [f'level1_{i}' for i in range(1, n_rows_1+1)]
        lev_1.insert(0, column='name', value=new_cols_1)
        lev_1_relative.insert(0, column='name', value=new_cols_1)
        lev_1.to_csv(lev_1_a,sep='\t')
        lev_1_relative.to_csv(lev_1_r,sep='\t')
        
        ## fun_3
        lefse_dir = args.outDir + "/lefse"
        if not os.path.exists(lefse_dir):
            os.mkdir(lefse_dir)
        lefse_sub_dir = "{}/{}".format(lefse_dir,g_name)
        if not os.path.exists(lefse_sub_dir):
            os.mkdir(lefse_sub_dir)
        group = pd.read_csv(self.group,sep="\t",index_col=0)
        group_sub = group[[g_name]]
        # pathway的level3
        df1 = pd.read_csv(pathway_r, sep='\t')
        merged_df = self.lefse_format(df=df1,group=group_sub)
        df1_out = lefse_sub_dir + "/lefse_pathway_lev3.txt"
        merged_df.to_csv(df1_out, sep='\t', header=False)

        # pathway的level2
        df2 = pd.read_csv(lev_2_r, sep='\t')
        df2.drop(columns=df2.columns[0], axis=1, inplace=True)
        merged_df = self.lefse_format(df=df2,group=group_sub)
        df2_out = lefse_sub_dir + "/lefse_pathway_lev2.txt"
        merged_df.to_csv(df2_out, sep='\t', header=False)
        
        # pathway的level1
        df3 = pd.read_csv(lev_1_r, sep='\t')
        df3.drop(columns=df3.columns[0], axis=1, inplace=True)
        merged_df = self.lefse_format(df=df3,group=group_sub)
        df3_out = lefse_sub_dir + "/lefse_pathway_lev1.txt"
        merged_df.to_csv(df3_out, sep='\t', header=False)
        
        # KO
        df4 = pd.read_csv(ko_r, sep='\t')
        merged_df = self.lefse_format(df=df4,group=group_sub)
        df4_out = lefse_sub_dir + "/lefse_KO.txt"
        merged_df.to_csv(df4_out, sep='\t', header=False)
        
        ## fun_4
        # phylum -> species
        tmpFile = args.inDir + "/tmp_lefse.txt"
        cmd1 = "grep -E -v \"^#|t__\" {} > {}".format(self.sp_raw,tmpFile)
        os.system(command=cmd1)
        df5 = pd.read_csv(tmpFile, sep='\t')
        merged_df = self.lefse_format(df=df5,group=group_sub)
        df5_out = lefse_sub_dir + "/lefse_sp.txt"
        merged_df.to_csv(df5_out, sep='\t', header=False)
        cmd2 = "rm -rf {}".format(tmpFile)
        os.system(command=cmd2)

    def lefse_run(self,g_name,lda):
        """
        @功能: 运行lefse:
        @参数g_name: 确定哪个分组文件[字符串]
        @参数lda: lda值
        """
        lefse_dir = "{}/lefse/{}".format(args.outDir,g_name)
        lefse_list = ["lefse_KO","lefse_pathway_lev1","lefse_pathway_lev2","lefse_pathway_lev3","lefse_sp"]
        for i in lefse_list:
            lefse_in = "{}/{}.txt".format(lefse_dir,i)
            lefse_out = "{}/{}".format(lefse_dir,i)
            cmd = "{} {} -i {} -o {} -l {}".format(self.py_env,
                                                   self.py_lefse,
                                                   lefse_in,
                                                   lefse_out,
                                                   lda)
            os.system(command=cmd)
        
    def main(self):
        self.pre_format(g_name=args.group)
        self.lefse_run(g_name=args.group,lda=args.lda)
        
if __name__ == "__main__":
    run = AfterPipe()
    run.main()
    