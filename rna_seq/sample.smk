# ======================================================================================================================
# Project: BioPipe
# Script : sample.smk.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.27
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
import pandas as pd
import gzip
caseinfo=pd.read_csv("conf/sample.csv").set_index(["sample", "unit"], drop=False)
for index,info in caseinfo.iterrows():
    with gzip.open(info["fq1"]) as f:
        f.readline()
        this_len=len(f.readline().decode().rstrip("\n"))
    caseinfo.loc[index,"read_len"]=this_len

    # print(index)
    # print(info)

# print(caseinfo)