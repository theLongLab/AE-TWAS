##20210109
##generate position range information for the genes in a module
##will be used to extract SNPs information from bed,bim,fam file

#!/bin/python
import pandas as pd
import sys
import re
import os
from os import path

index_list=["34","59","66","72","74","75"]
for index in index_list:
    gtex_dict={"34":"34_Div8_ID670_R0.95","59":"59_Div8_ID670_R0.85","66":"66_Div8_ID670_R0.85","72":"72_Div8_ID670_R0.94","74":"74_Div8_ID670_R0.88","75":"75_Div8_ID670_R0.88"}
    wk="./WGCNA_12K_6modules"
    file_name1 = wk+"/"+str(index)+".txt"
    gene_name_list=[]
    with open(file_name1,"r") as f:
        line=f.readline().strip().replace('\"','')
        while line:
            if line not in gene_name_list:
                gene_name_list.append(line)
            line=f.readline().strip().replace('\"','')
        
    print("Number of genes in the file list is ",len(gene_name_list))

    gtf=pd.read_table("./gencode.v26.annotation.3genes.gtf", sep='\t',header=None)       


    Folder1 = str(index)+"_genes_1M_X_bslmm"
    if not os.path.exists(wk+'/'+str(index)+"/"+Folder1):
        new = wk+'/'+str(index)+"/"+Folder1
        cmd = 'mkdir -p '+new
        os.system(cmd)
        
    Folder2 = str(index)+"_genes_1M_X_apo_bslmm"
    if not os.path.exists(wk+'/'+str(index)+"/"+Folder2):
        new = wk+'/'+str(index)+"/"+Folder2
        cmd = 'mkdir -p '+new
        os.system(cmd)

    sbatch_cmd = open(wk+'/'+str(index)+"/"+str(index)+"_genes_1M_bslmm_all.cmd","w")   

    num=len(gene_name_list)
    jobs_cmd=500
    slurm_gene=[]
    for i in range(num):
        slurm_gene.append(i)
        if (i>0 and (i% jobs_cmd)==0) or (i==(num-1)):
            batch_num = int(i/jobs_cmd)
            if ((i% jobs_cmd)!=0) and i==(num-1):
                batch_num = int(i/jobs_cmd) +1
            cmd_file = open(wk+'/'+str(index)+"/"+str(index)+"_genes_1M_bslmm"+str(batch_num)+".cmd","w")
            cmd_file.write('#!/bin/bash\n')
            cmd_file.write('#SBATCH --job-name='+str(index)+'_b'+str(batch_num)+'_bslmm\n')
            cmd_file.write('#SBATCH --chdir='+wk+'/'+str(index)+"/"+Folder1+'\n')
            cmd_file.write('#SBATCH --error='+str(index)+'_b'+str(batch_num)+'_bslmm.error\n')
            cmd_file.write('#SBATCH --out='+str(index)+'_b'+str(batch_num)+'_bslmm.out\n')
            cmd_file.write('#SBATCH --mem=100G\n')
            cmd_file.write('#SBATCH --nodes=1\n')
            cmd_file.write('#SBATCH --ntasks=1\n')
            cmd_file.write('#SBATCH --cpus-per-task=1\n')
            cmd_file.write('#SBATCH --time=1-0:0:0\n')
            cmd_file.write('#SBATCH --partition=bigmem\n')
            cmd_file.write('echo $(date)\n')
            for j in range(len (slurm_gene)):
                gene_index = slurm_gene[j]
                cmd_file.write('./plink --bfile ./WGCNA_12K_6modules/GTEx_WB_MSSNG --extract '+gene_name_list[gene_index]+'.txt --range --make-bed --out '+gene_name_list[gene_index]+'\n')
                cmd_file.write("ln -s "+wk+"/"+str(index)+"/"+Folder1+"/"+gene_name_list[gene_index]+".bed "+wk+"/"+str(index)+"/"+Folder2+"/"+gene_name_list[gene_index]+".bed\n")
                cmd_file.write("ln -s "+wk+"/"+str(index)+"/"+Folder1+"/"+gene_name_list[gene_index]+".bim "+wk+"/"+str(index)+"/"+Folder2+"/"+gene_name_list[gene_index]+".bim\n")
                cmd_file.write("cp "+wk+"/"+str(index)+"/"+Folder1+"/"+gene_name_list[gene_index]+".fam "+wk+"/"+str(index)+"/"+Folder2+"/"+gene_name_list[gene_index]+".fam\n")
                
                #change fam file to BSLMM fam file in X
                cmd_file.write('python ./Generate_BSLMM_fam.py '+gtex_dict.get(index)+' '+gene_name_list[gene_index]+' '+Folder1+'\n')
                #change fam file to BSLMM fam file in X_apo
                cmd_file.write('python ./Generate_BSLMM_fam.py '+gtex_dict.get(index)+' '+gene_name_list[gene_index]+' '+Folder2+'\n')
            
            cmd_file.write('echo Finishing $(date)')    
            cmd_file.close()
            slurm_gene = []
            sbatch_cmd.write("sbatch "+wk+'/'+str(index)+"/"+str(index)+"_genes_1M_bslmm"+str(batch_num)+".cmd\n")
            
    sbatch_cmd.close()
