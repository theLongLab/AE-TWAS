##### Calculation of gene pairs’ correlation 
setwd("data/pair_gene_cor_result/")
module_name = "34"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
mean(abs(input_corr$V1))  
mean(abs(output_corr$V1))  
sum(abs(input_corr$V1) < abs(output_corr$V1))/dim(input_corr)[1] 
1-sum(abs(input_corr$V1) < abs(output_corr$V1))/dim(input_corr)[1]  
## replace the module_name with 59, 66, 72, 74 and 75 to get results of the rest modules

##### Calculation of gene’s connectivity
library(WGCNA)
library(dplyr)
### input
setwd("data/gene_connectivity/")
gene_expression_wb = read.csv("gene_expression_input.csv",header = T)
softPower = 16
adjacency_wb = adjacency(gene_expression_wb, power = softPower)
dim(adjacency_wb) # 8191 8191
sum(colnames(adjacency_wb) == colnames(gene_expression_wb)) # 8191
colors = c(rep("black",4715),rep("red",74),rep("blue",68),
           rep("green",3268),rep("orange",35),rep("white",31))
IC = intramodularConnectivity(adjacency_wb, colors, scaleByMax = FALSE)
IC$gene = row.names(IC)
summary(IC$kTotal)
module_name = c(rep("34",4715),rep("59",74),rep("66",68),
                rep("72",3268),rep("74",35),rep("75",31))
module_name2 = data.frame(colnames(gene_expression_wb),module_name)
colnames(module_name2) = c("gene","module")
IC2 = merge(x=IC,y=module_name2,by="gene")
IC2 = IC2[order(desc(IC2$kTotal)),]
write.table(IC2,"gene_connectivity_input.csv",sep=",",row.names=FALSE,quote = FALSE)

### output
gene_expression_wb = read.csv("gene_expression_output.csv",header = T)
softPower = 16
adjacency_wb = adjacency(gene_expression_wb, power = softPower)
dim(adjacency_wb) # 8191 8191
sum(colnames(adjacency_wb) == colnames(gene_expression_wb)) # 8191
colors = c(rep("black",4715),rep("red",74),rep("blue",68),
           rep("green",3268),rep("orange",35),rep("white",31))
IC = intramodularConnectivity(adjacency_wb, colors, scaleByMax = FALSE)
IC$gene = row.names(IC)
summary(IC$kTotal)
module_name = c(rep("34",4715),rep("59",74),rep("66",68),
                rep("72",3268),rep("74",35),rep("75",31))
module_name2 = data.frame(colnames(gene_expression_wb),module_name)
colnames(module_name2) = c("gene","module")
IC2 = merge(x=IC,y=module_name2,by="gene")
IC2 = IC2[order(desc(IC2$kTotal)),]
write.table(IC2,"gene_connectivity_output.csv",sep=",",row.names=FALSE,quote = FALSE)

### merge
IC_input = read.csv("gene_connectivity_input.csv",header = T)
colnames(IC_input) = c("gene","kTotal_input","kWithin_input","kOut_input","kDiff_input","module")
IC_output = read.csv("gene_connectivity_output.csv",header = T)
colnames(IC_output) = c("gene","kTotal_output","kWithin_output","kOut_output","kDiff_output","module")
summary(IC_input$kWithin)
summary(IC_output$kWithin)
IC_merge = merge(x=IC_input,IC_output,by=c("gene","module"))
sum(IC_merge$kTotal_input < IC_merge$kTotal_output) # 3159
colnames(IC_merge)[1] = "gene_id"
write.table(IC_merge,"gene_connectivity_merge.csv",sep=",",row.names=FALSE,quote = FALSE)





