#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   lefse.py
@Time    :   2023/04/23 13:22:09
@Author  :   Wanjin.Hu 
@Version :   1.0
@Contact :   wanjin.hu@diprobio.com
@Description :  cladogram的结果用AI调整下画布大小
'''

import argparse
import os

parser = argparse.ArgumentParser(prog="Lefse",
                                 usage="",
                                 description="Lefse analyse")
parser.add_argument("-i", "--input", dest="input", required=True, type=str, help="lefse input file, such as: lefse_input.txt")
parser.add_argument("-o", "--output_dir", dest="outDir", required=False, type=str, default="lefse", help="result output dir")
parser.add_argument("-l", "--lda", dest="lda", required=False, type=str, default=2.0, help="LDA score")
args = parser.parse_args()

class Lefse(object):
    def __init__(self):
        super(Lefse,self).__init__()
        self.bash_env = "/usr/bin/bash"
        self.cmd = "/root/pipe_script/general/lefse/lefse.sh"
        if os.path.exists(args.outDir):
            print("{} 文件夹已存在,请删除后再进行分析!".format(args.outDir))
            quit()
        if not os.path.exists(args.outDir):
            os.mkdir(args.outDir)
    
    def lefse_run(self):
        """
        @功能: 
        """
        lefse_in = args.outDir + "/lefse_input.in"
        lefse_res = args.outDir + "/lefse_input.res"
        lefse_diff = args.outDir + "/lefse_diff.pdf"
        lefse_cladogram = args.outDir + "/lefse_cladogram.pdf"
        lefse_single_dir = args.outDir + "/lefse_single"
        cmd = "{} {} {} {} {} {} {} {} {}".format(self.bash_env,
                                                  self.cmd,
                                                  args.input,
                                                  lefse_in,
                                                  lefse_res,
                                                  lefse_diff,
                                                  lefse_cladogram,
                                                  lefse_single_dir,
                                                  args.lda)
        os.system(command=cmd)
        print("Lefse finished!")

    def main(self):
        self.lefse_run()

if __name__ == "__main__":
    run = Lefse()
    run.main()