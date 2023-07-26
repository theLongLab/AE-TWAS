##### Functional verification of discovered genes
##### part 1
##### identify by protocol
## replace the disease name with CD, ASD, RA and T1D to get results of the rest diseases
setwd("data/real_data/")
disease = "SCZ"
AE_TWAS = read.csv(paste(disease,"/",disease,"all_res_apo_FDR0.05.csv",sep=""),header = T)
sum(AE_TWAS$p_corr <= 0.05) 
standard_TWAS = read.csv(paste(disease,"/",disease,"all_res_FDR0.05.csv",sep=""),header = T)
sum(standard_TWAS$p_corr <= 0.05)
AE_TWAS_sig_gene = AE_TWAS[which(AE_TWAS$p_corr <= 0.05),]
write.table(AE_TWAS_sig_gene,paste("Disgenet/",disease,"_AE_TWAS_sig_gene.csv",sep=""),sep=",",row.names=FALSE,quote = FALSE)
standard_TWAS_sig_gene = standard_TWAS[which(standard_TWAS$p_corr <= 0.05),]
write.table(standard_TWAS_sig_gene,paste("Disgenet/",disease,"_standard_TWAS_sig_gene.csv",sep=""),sep=",",row.names=FALSE,quote = FALSE)

##### identify by disgenet
CUI = "C0036341"
## replace the CUI number with that of CD, ASD, RA and T1D to get results of the rest diseases
disgenet = read.csv(paste("Disgenet/",disease,"_",CUI,"_gda.csv",sep=""),header = T)

AE_unique_total = setdiff(AE_TWAS_sig_gene$gene_symbol,standard_TWAS_sig_gene$gene_symbol)
length(AE_unique_total) 
standard_unique_total = setdiff(standard_TWAS_sig_gene$gene_symbol,AE_TWAS_sig_gene$gene_symbol)
length(standard_unique_total) 
AE_standard_intersect = intersect(AE_TWAS_sig_gene$gene_symbol,standard_TWAS_sig_gene$gene_symbol)
length(AE_standard_intersect) 

AE_success = intersect(AE_TWAS_sig_gene$gene_symbol,disgenet$Gene)
length(AE_success)
standard_success = intersect(standard_TWAS_sig_gene$gene_symbol,disgenet$Gene)
length(standard_success) 
AE_unique = setdiff(AE_success,standard_success)
length(AE_unique)  
standard_unique = setdiff(standard_success,AE_success)
length(standard_unique)

##### part 2: hub genes in DisGeNET
## select hub genes using same cut-off for each module, keep different gene number before and after AE
# 2049=1179+19+17+817+9+8 
# 2350=680+74+68+1462+35+31 
same25cutoff_before = read.csv("data/gene_connectivity/same25cutoff/same25cutoff_before.csv",header = T)
same25cutoff_after = read.csv("data/gene_connectivity/same25cutoff/same25cutoff_after.csv",header = T)
A = same25cutoff_before$gene # before total 2049
B = same25cutoff_after$gene # after total 2350
C = setdiff(A,B) # before unique 554
D = setdiff(B,A) # after unique 855

setwd("data/real_data/Disgenet/")
disease = "SCZ"  
## replace the disease name with CD, ASD, RA and T1D to get results of the rest diseases
CUI = "C0036341"   
## replace the CUI number with that of CD, ASD, RA and T1D to get results of the rest diseases
real_data_disgenet = read.csv(paste(disease,"_",CUI,"_gda.csv",sep=""),header = T)
E = real_data_disgenet$Gene # disgenet all
success_before = length(intersect(A,E))
success_before
success_after = length(intersect(B,E))
success_after
success_before_unique = length(intersect(C,E))
success_before_unique
success_after_unique = length(intersect(D,E))
success_after_unique

##### part 3: hub genes in TWAS
setwd("data/real_data/")
disease = "SCZ"
real_data_supp = read.csv(paste(disease,"_supp_table.csv",sep=""),header = T)
real_data_supp_sub = real_data_supp[,c(3,7,8)]

AE_TWAS = real_data_supp_sub[which(real_data_supp_sub$AE.TWAS != " NA"),]
FF = AE_TWAS$Gene # AE-TWAS all
standard_TWAS = real_data_supp_sub[which(real_data_supp_sub$Standard.TWAS != " NA"),]
G = standard_TWAS$Gene # standard all
AE_TWAS_unique = real_data_supp_sub[which((real_data_supp_sub$AE.TWAS != " NA") & (real_data_supp_sub$Standard.TWAS == " NA")),]
H = AE_TWAS_unique$Gene # AE-TWAS unique
standard_TWAS_unique = real_data_supp_sub[which((real_data_supp_sub$AE.TWAS == " NA") & (real_data_supp_sub$Standard.TWAS != " NA")),]
I = standard_TWAS_unique$Gene # standard unique

length(intersect(A,G))
length(intersect(B,FF))
length(intersect(C,I))
length(intersect(D,H))






