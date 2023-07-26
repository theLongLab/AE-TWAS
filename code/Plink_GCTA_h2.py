##20201222 v1 
##Generate cmd files to calcuate h^2 for each gene using GCTA

#!/usr/bin/python
import sys
import random
import heapq
import os
import sys

def process(m, n, t):
    ##read gene names from txt file
    wk='./PathToDict/'
    index = str(m)
    Folder = index+"_genes_1M_seperate_plink_files"
    file_name1 = wk +str(index)+".txt"
    gene_name_list=[]
    with open(file_name1,"r") as f:
        line=f.readline().strip().replace('\"','')
        while line:
            if line not in gene_name_list:
                gene_name_list.append(line)
            line=f.readline().strip().replace('\"','')
    print("Number of genes in the file list is ",len(gene_name_list))
        
    ##generate cmd files to run heritability 
    num = int (n)
    jobs_cmd= int (t)
    if not os.path.exists(wk+index+"/genes_gcta_cmd_2"):
        cmd = 'mkdir -p '+ wk+index+"/genes_gcta_cmd_2"
        os.system(cmd)
    
    slurm_gene=[]
    for i in range(num):
        slurm_gene.append(i)
        if (i>0 and (i% jobs_cmd)==0) or (i==(num-1)):
            f = open(wk+index+"/genes_gcta_cmd_2/"+str(i)+ '.cmd', 'w')
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --job-name=gcta'+index+str(i)+'\n')
            f.write('#SBATCH --error=gcta'+index+str(i)+'.error\n')
            f.write('#SBATCH --out=gcta'+index+str(i)+'.out\n')
            f.write('#SBATCH --mem-per-cpu=10G\n')
            f.write('#SBATCH --nodes=1\n')
            f.write('#SBATCH --ntasks=1\n')
            f.write('#SBATCH --cpus-per-task=2\n')
            f.write('#SBATCH --time=0-23:23:59\n')
            f.write('#SBATCH --partition=theia,parallel\n')
            for j in range(len (slurm_gene)):
                gene_index = slurm_gene[j]
                pheno_index = gene_index+1
                ss1 = './gcta_1.93.2beta/gcta64 --bfile '+gene_name_list[gene_index]+' --make-grm --out '+gene_name_list[gene_index]
                ss2 = './gcta_1.93.2beta/gcta64 --reml --grm '+gene_name_list[gene_index]+' --pheno '+wk+'/'+str(index)+'/M'+index+'X'+' --mpheno '+str(pheno_index)+' --out X.'+gene_name_list[gene_index]
                ss3 = './gcta_1.93.2beta/gcta64 --reml --grm '+gene_name_list[gene_index]+' --pheno '+wk+'/'+str(index)+'/M'+index+'X_apo_2'+' --mpheno '+str(pheno_index)+' --out X_apo_2.'+gene_name_list[gene_index]
                f.write(ss1+ '\n')
                f.write(ss2+ '\n')
                f.write(ss3+ '\n')
            f.close()
            slurm_gene = []
    
if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("inputs: module number, gene number in the module, gene number in each cmd")
    else:
        process(sys.argv[1], sys.argv[2],sys.argv[3])
