#!/usr/bin/python
import sys
import os


def process(blsmm_pheno,predixcan_pheno):
    IDs=[]
    f_pheno=open("./GTEX_AE/Wholeblood/WGCNA_12K_6modules/GTEx_WB_MSSNG.fam","r")
    line_pheno=f_pheno.readline().strip()
    while line_pheno:
        ID = line_pheno.split(" ")[0]
        IDs.append(ID)
        line_pheno=f_pheno.readline()
    f_pheno.close()
    
    count=0
    fw=open(predixcan_pheno,"w")
    fw.write("#FID\tIID\tP\n")
    f_bslmm = open(blsmm_pheno,"r")
    f_bslmm_line = f_bslmm.readline().strip()
    while f_bslmm_line:
        if f_bslmm_line !='NA':
            fw.write(IDs[count]+'\t'+IDs[count]+'\t'+f_bslmm_line+'\n')
        count+=1
        f_bslmm_line=f_bslmm.readline().strip()
    if(count != len(IDs)):
        print(str(count)+"\t"+str(len(IDs))+"\tNumber of IDs in predixcan does not match that in real pheno")
    fw.close()
    f_bslmm.close()


if __name__ == "__main__":
    if(len(sys.argv)==1):
        print("two inputs: [1] a file path for bslmm predicted pheno; [2] a file path for predixcan pheno") 
    else:
        process(sys.argv[1],sys.argv[2])
        