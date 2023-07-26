#!/bin/python
#20210109

import pandas as pd
import sys


def fam_files(X_path,pheno_file_name,target_folder):
    wk="/work/long_lab/qli/Autoencoder/GTEX_AE/Wholeblood/WGCNA_12K_6modules/"
    gtex_pheno_df = pd.read_csv(X_path,index_col=0,sep=" ") #set the 1st line as index
    name_index = list(gtex_pheno_df.index) #get index to list
    name_cols = gtex_pheno_df.columns.tolist() #change columns name to list 
    fam=pd.read_csv(target_folder+"/"+pheno_file_name+'.fam',sep=" ",header=None)
    fam.columns=["FID","IID","F","M","SEX","P"] #rename fam columns names
    fam_index = fam["IID"].tolist()
    position = []
    for item in name_index:
        if item in fam_index:
            position.append(fam_index.index(item)) 
    
    gtex_gene_index = name_cols.index(pheno_file_name)
    for i in range(len(position)):
        fam.loc[position[i],"P"] = gtex_pheno_df.iloc[i,gtex_gene_index]
    fam.to_csv(target_folder+"/"+pheno_file_name+'.fam',sep=" ",header=False,index=False)


def process(a,b,c):
    M = a #66_Div8_ID670_R0.85
    wk = "/work/long_lab/qli/Autoencoder/GTEX_AE/Wholeblood/WGCNA_12K_6modules/"
    target_folder=wk+M.split("_")[0]+"/"+c
    X_path = wk+"gene_expression_X_h_X_apo/M"+M+"_input_head_gcta.txt" 
    X_apo_path = wk+"gene_expression_X_h_X_apo/M"+M+"_output_head_gcta.txt"
    pheno_file_name=b
    
    if("apo" in target_folder):
        fam_files(X_apo_path,pheno_file_name,target_folder)
    else:
        fam_files(X_path,pheno_file_name,target_folder)


if __name__ == "__main__":
    if(len(sys.argv)<=1):
        print("three inputs: [1]:66_Div8_ID670_R0.85; [2]:name of the phenotype file; [3]:name of the final folder")  #66_Div8_ID670_R0.85 ENSG00000001036.13 34_genes_1M_X_bslmm/34_genes_1M_X_apo_bslmm
    else:
        process(sys.argv[1],sys.argv[2],sys.argv[3]) ##66_Div8_ID670_R0.85