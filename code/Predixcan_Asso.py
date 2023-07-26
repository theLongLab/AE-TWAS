#!/usr/bin/python
import sys
import os


def process(gene_file_name_prefix):
    wk="./WGCNA_12K_6modules"
    file_name1 = wk+"/"+str(gene_file_name_prefix)+".txt"
    index=gene_file_name_prefix
    gene_name_list=[]
    with open(file_name1,"r") as f:
        line=f.readline().strip().replace('\"','')
        while line:
            if line not in gene_name_list:
                gene_name_list.append(line)
            line=f.readline().strip().replace('\"','')
    
    Folder1 = str(index)+"_genes_1M_X_predixcan"
    if not os.path.exists(wk+'/'+str(index)+"/"+Folder1):
        new = wk+'/'+str(index)+"/"+Folder1
        cmd = 'mkdir -p '+new
        os.system(cmd)

    Folder2 = str(index)+"_genes_1M_X_apo_predixcan"
    if not os.path.exists(wk+'/'+str(index)+"/"+Folder2):
        new = wk+'/'+str(index)+"/"+Folder2
        cmd = 'mkdir -p '+new
        os.system(cmd)
        
    sbatch_cmd = open(wk+'/'+str(index)+"/"+str(index)+"_genes_1M_predixcan_all.cmd","w")    
    num=len(gene_name_list)
    jobs_cmd=500
    slurm_gene=[]
    for i in range(num):
        slurm_gene.append(i)
        if (i>0 and (i% jobs_cmd)==0) or (i==(num-1)):
            batch_num = int(i/jobs_cmd)
            if ((i% jobs_cmd)!=0) and i==(num-1):
                batch_num = int(i/jobs_cmd) +1
                
            cmd_file = open(wk+'/'+str(index)+"/"+str(index)+"_genes_1M_predixcan"+str(batch_num)+".cmd","w")
            cmd_file.write('#!/bin/bash\n')
            cmd_file.write('#SBATCH --job-name='+str(index)+'_p'+str(batch_num)+'\n')
            cmd_file.write('#SBATCH --error='+str(index)+'_p'+str(batch_num)+'.error\n')
            cmd_file.write('#SBATCH --out='+str(index)+'_p'+str(batch_num)+'.out\n')
            cmd_file.write('#SBATCH --mem=10G\n')
            cmd_file.write('#SBATCH --nodes=1\n')
            cmd_file.write('#SBATCH --ntasks=1\n')
            cmd_file.write('#SBATCH --cpus-per-task=1\n')
            cmd_file.write('#SBATCH --time=7-0:0:0\n')
            cmd_file.write('#SBATCH --partition=single,lattice,parallel,cpu2013,cpu2019,theia\n')
            for j in range(len (slurm_gene)):
                gene_index = slurm_gene[j] 
                cmd_file.write("python ./Bslmm_2_predixcan_pheno.py "+wk+'/'+str(index)+"/"+str(index)+'_genes_1M_X_bslmm/output/'+gene_name_list[gene_index]+".prdt.txt "+wk+'/'+str(index)+"/"+Folder1+"/"+gene_name_list[gene_index]+".predixcan.txt\n")
                cmd_file.write("python ./PrediXcan/Software/PrediXcanExample/PrediXcan.py --assoc --pheno "+wk+"/MSSNG.fam --pred_exp "+wk+'/'+str(index)+"/"+Folder1+"/"+gene_name_list[gene_index]+".predixcan.txt --logistic --output_prefix "+gene_name_list[gene_index]+'\n')

            cmd_file.close()
            cmd_file1 = open(wk+'/'+str(index)+"/"+str(index)+"_genes_1M_predixcan"+str(batch_num)+".apo.cmd","w")
            cmd_file1.write('#!/bin/bash\n')
            cmd_file1.write('#SBATCH --job-name='+str(index)+'_p'+str(batch_num)+'.apo\n')
            cmd_file1.write('#SBATCH --error='+str(index)+'_p'+str(batch_num)+'.apo.error\n')
            cmd_file1.write('#SBATCH --out='+str(index)+'_p'+str(batch_num)+'.apo.out\n')
            cmd_file1.write('#SBATCH --mem=10G\n')
            cmd_file1.write('#SBATCH --nodes=1\n')
            cmd_file1.write('#SBATCH --ntasks=1\n')
            cmd_file1.write('#SBATCH --cpus-per-task=1\n')
            cmd_file1.write('#SBATCH --time=7-0:0:0\n')
            cmd_file1.write('#SBATCH --partition=single,lattice,parallel,cpu2013,cpu2019,theia\n')
            for j in range(len (slurm_gene)):
                gene_index = slurm_gene[j] 
                cmd_file1.write("python ./Bslmm_2_predixcan_pheno.py "+wk+'/'+str(index)+"/"+str(index)+'_genes_1M_X_apo_bslmm/output/'+gene_name_list[gene_index]+".prdt.txt "+wk+'/'+str(index)+"/"+Folder2+"/"+gene_name_list[gene_index]+".predixcan.txt\n")
                cmd_file1.write("python ./PrediXcan/Software/PrediXcanExample/PrediXcan.py --assoc --pheno "+wk+"/MSSNG.fam --pred_exp "+wk+'/'+str(index)+"/"+Folder2+"/"+gene_name_list[gene_index]+".predixcan.txt --logistic --output_prefix "+gene_name_list[gene_index]+'\n')
            
            cmd_file1.close()
            
            slurm_gene = []
            sbatch_cmd.write("sbatch "+wk+'/'+str(index)+"/"+str(index)+"_genes_1M_predixcan"+str(batch_num)+".cmd\n")
            sbatch_cmd.write("sbatch "+wk+'/'+str(index)+"/"+str(index)+"_genes_1M_predixcan"+str(batch_num)+".apo.cmd\n")

    sbatch_cmd.close()



if __name__ == "__main__":
    if(len(sys.argv)==1):
        print("one input: [1] a file prefix contains all gene's name (34)")
    else:
        process(sys.argv[1])
