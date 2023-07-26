##### Calculation of expression heritability
setwd("data/expression_heritability/")
X34 = read.csv("X_34.csv",header = T)
X59 = read.csv("X_59.csv",header = T)
X66 = read.csv("X_66.csv",header = T)
X72 = read.csv("X_72.csv",header = T)
X74 = read.csv("X_74.csv",header = T)
X75 = read.csv("X_75.csv",header = T)
X = rbind(X34,X59,X66,X72,X74,X75)

X34apo = read.csv("X_34_apo.csv",header = T)
X59apo = read.csv("X_59_apo.csv",header = T)
X66apo = read.csv("X_66_apo.csv",header = T)
X72apo = read.csv("X_72_apo.csv",header = T)
X74apo = read.csv("X_74_apo.csv",header = T)
X75apo = read.csv("X_75_apo.csv",header = T)
X_apo = rbind(X34apo,X59apo,X66apo,X72apo,X74apo,X75apo)

h2 = data.frame(X$Variance.value.,X_apo$Variance.value.) # 8191
colnames(h2) = c("X","X_apo")
write.table(h2,"h2.csv",sep=",",quote = FALSE,row.names = FALSE)
sum(h2$X == h2$X_apo) # 251 unchanged 
sum(h2$X < h2$X_apo) # 3965 increased
sum(h2$X > h2$X_apo) # 3975 decreased
# 251+3965+3975 = 8191
round(mean(h2$X),3) # 0.061
round(mean(h2$X_apo),3) # 0.057
h2_update = h2[which(h2$X != h2$X_apo),]

h2_update_interval1 = h2_update[which((h2_update$X >= 0) & (h2_update$X < 0.05)),] 
dim(h2_update_interval1) # 3801
sum(h2_update_interval1$X < h2_update_interval1$X_apo) # 2921
sum(h2_update_interval1$X > h2_update_interval1$X_apo) # 880
round(mean(h2_update_interval1$X),4) # 0.0189
round(mean(h2_update_interval1$X_apo),4) # 0.0543
round(2921/3801,4) # 0.7685
round(880/3801,4) # 0.2315

h2_update_interval1_1 = h2_update[which((h2_update$X >= 0) & (h2_update$X < 0.01)),] 
dim(h2_update_interval1_1) # 1531
sum(h2_update_interval1_1$X < h2_update_interval1_1$X_apo) # 1452
sum(h2_update_interval1_1$X > h2_update_interval1_1$X_apo) # 79
round(mean(h2_update_interval1_1$X),4) # 0.0014
round(mean(h2_update_interval1_1$X_apo),4) # 0.0552
round(1452/1531,4) # 0.9484
round(79/1531,4) # 0.0516

h2_update_interval1_2 = h2_update[which((h2_update$X >= 0.01) & (h2_update$X < 0.025)),] 
dim(h2_update_interval1_2) # 814
sum(h2_update_interval1_2$X < h2_update_interval1_2$X_apo) # 605
sum(h2_update_interval1_2$X > h2_update_interval1_2$X_apo) # 209
round(mean(h2_update_interval1_2$X),4) # 0.0176
round(mean(h2_update_interval1_2$X_apo),4) # 0.0528
round(605/814,4) # 0.7432
round(209/814,4) # 0.2568

h2_update_interval1_3 = h2_update[which((h2_update$X >= 0.025) & (h2_update$X < 0.05)),] 
dim(h2_update_interval1_3) # 1456
sum(h2_update_interval1_3$X < h2_update_interval1_3$X_apo) # 864
sum(h2_update_interval1_3$X > h2_update_interval1_3$X_apo) # 592
round(mean(h2_update_interval1_3$X),4) # 0.0379
round(mean(h2_update_interval1_3$X_apo),4) # 0.0543
round(864/1456,4) # 0.5934
round(592/1456,4) # 0.4066

h2_update_interval2 = h2_update[which((h2_update$X >= 0.05) & (h2_update$X < 0.1)),] 
dim(h2_update_interval2) # 2398
sum(h2_update_interval2$X < h2_update_interval2$X_apo) # 861
sum(h2_update_interval2$X > h2_update_interval2$X_apo) # 1537
round(mean(h2_update_interval2$X),4) # 0.0727
round(mean(h2_update_interval2$X_apo),4) # 0.0592
round(861/2398,4) # 0.359
round(1537/2398,4) # 0.641

h2_update_interval3 = h2_update[which((h2_update$X >= 0.1) & (h2_update$X < 0.15)),] 
dim(h2_update_interval3) # 1187
sum(h2_update_interval3$X < h2_update_interval3$X_apo) # 167
sum(h2_update_interval3$X > h2_update_interval3$X_apo) # 1020
round(mean(h2_update_interval3$X),4) # 0.1215
round(mean(h2_update_interval3$X_apo),4) # 0.0671
round(167/1187,4) # 0.1407
round(1020/1187,4) # 0.8593

h2_update_interval4 = h2_update[which((h2_update$X >= 0.15) & (h2_update$X < 0.2)),] 
dim(h2_update_interval4) # 395
sum(h2_update_interval4$X < h2_update_interval4$X_apo) # 12
sum(h2_update_interval4$X > h2_update_interval4$X_apo) # 383
round(mean(h2_update_interval4$X),4) # 0.1698
round(mean(h2_update_interval4$X_apo),4) # 0.0727
round(12/395,4) # 0.0304
round(383/395,4) # 0.9696

h2_update_interval5 = h2_update[which((h2_update$X >= 0.2) & (h2_update$X <= 1)),] 
dim(h2_update_interval5) # 159
sum(h2_update_interval5$X < h2_update_interval5$X_apo) # 4
sum(h2_update_interval5$X > h2_update_interval5$X_apo) # 155
round(mean(h2_update_interval5$X),4) # 0.2409
round(mean(h2_update_interval5$X_apo),4) # 0.0717
round(4/159,4) # 0.0252
round(155/159,4) # 0.9748

