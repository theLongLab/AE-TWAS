library(ggplot2)
library(ggpubr)
library("cowplot")
library(plyr)
library(reshape2)

##### figure 2
setwd("data/Figures/")
data1 = read.csv("figure2a.csv",header = T)
colnames(data1) = c("interval","increase","decrease")
a = c("[0, 0.001)","[0.001, 0.025)","[0.025, 0.05)","[0.05, 0.1)","[0.1, 0.15)","[0.15, 0.2)","[0.2, 1]")
data2 = data.frame(group = rep(c("Increased", "Decreased"), each=7), x = rep(a,2), y = c(data1$increase,-data1$decrease))
f2a = ggplot(data2, aes(x=x, y=y, fill=group,label = y)) + 
  geom_bar(stat="identity", position="identity",width=0.4) +
  geom_text(aes(label=y, vjust=-0.1),size=3) + 
  ylim(-1600,2000)+labs(y = "Gene Count", x = "Heritability Interval before AE Transformation")+
  scale_fill_manual(values = c("skyblue","darkblue")) +
  theme(text = element_text(size = 10),
        legend.position = c(0.86, 0.82),
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

data2b = read.csv("data/expression_heritability/h2.csv",header = T)
data2b_x = data.frame(data2b$X,rep("Before",dim(data2b)[1]))
colnames(data2b_x) = c("h2","time")
data2b_xpo = data.frame(data2b$X_apo,rep("After",dim(data2b)[1]))
colnames(data2b_xpo) = c("h2","time")
data2b_new = rbind(data2b_x, data2b_xpo)

f2b = ggplot(data2b_new, aes(x=h2, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "", x = "Gene Heritability")+ylim(0,2000)+
  #geom_vline(data=mu34, aes(xintercept=grp.mean, color=time),linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = c(0.86, 0.82),
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

ggarrange(f2a, f2b, labels = c("A", "B"), ncol = 2, nrow = 1)


##### figure 3
setwd("data/gene_pairs_correlation /")
module_name = "34"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
input_corr$time = "Before"
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
output_corr$time = "After"
data34 = rbind(input_corr,output_corr)
data34$V1 = abs(data34$V1)
mu34 = ddply(data34, "time", summarise, grp.mean=mean(V1))

module_name = "59"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
input_corr$time = "Before"
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
output_corr$time = "After"
data59 = rbind(input_corr,output_corr)
data59$V1 = abs(data59$V1)
mu59 = ddply(data59, "time", summarise, grp.mean=mean(V1))

module_name = "66"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
input_corr$time = "Before"
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
output_corr$time = "After"
data66 = rbind(input_corr,output_corr)
data66$V1 = abs(data66$V1)
mu66 = ddply(data66, "time", summarise, grp.mean=mean(V1))

module_name = "72"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
input_corr$time = "Before"
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
output_corr$time = "After"
data72 = rbind(input_corr,output_corr)
data72$V1 = abs(data72$V1)
mu72 = ddply(data72, "time", summarise, grp.mean=mean(V1))

module_name = "74"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
input_corr$time = "Before"
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
output_corr$time = "After"
data74 = rbind(input_corr,output_corr)
data74$V1 = abs(data74$V1)
mu74 = ddply(data74, "time", summarise, grp.mean=mean(V1))

module_name = "75"
input_corr = read.csv(paste(module_name,"_input_corr.csv",sep=""),header = F)
input_corr$time = "Before"
output_corr = read.csv(paste(module_name,"_output_corr.csv",sep=""),header = F)
output_corr$time = "After"
data75 = rbind(input_corr,output_corr)
data75$V1 = abs(data75$V1)
mu75 = ddply(data75, "time", summarise, grp.mean=mean(V1))

f3aa = ggplot(data34, aes(x=V1, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "Gene Pairs Count", x = "Module 34 (4715 genes)") + 
  geom_vline(data=mu34, aes(xintercept=grp.mean, color=time),
             linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = c(0.2, 0.8),
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

f3bb = ggplot(data59, aes(x=V1, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "", x = "Module 59 (74 genes)") + 
  geom_vline(data=mu59, aes(xintercept=grp.mean, color=time),
             linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

f3cc = ggplot(data66, aes(x=V1, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "", x = "Module 66 (68 genes)") + 
  geom_vline(data=mu66, aes(xintercept=grp.mean, color=time),
             linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

f3dd = ggplot(data72, aes(x=V1, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "Gene Pairs Count", x = "Module 72 (3268 genes)") + 
  geom_vline(data=mu72, aes(xintercept=grp.mean, color=time),
             linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

f3ee = ggplot(data74, aes(x=V1, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "", x = "Module 74 (35 genes)") + 
  geom_vline(data=mu74, aes(xintercept=grp.mean, color=time),
             linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.title=element_blank(),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

f3ff = ggplot(data75, aes(x=V1, fill = time, color=time)) +
  geom_histogram(binwidth=0.01, position="identity",alpha=0.7)+
  labs(y = "", x = "Module 75 (31 genes)") + 
  geom_vline(data=mu75, aes(xintercept=grp.mean, color=time),
             linetype="dashed") +
  scale_color_manual(values=c("darkcyan","bisque2")) + 
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme(text = element_text(size = 10),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey',linetype = 'dotted'))

ggarrange(f3aa,f3bb,f3cc,f3dd,f3ee,f3ff,labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)


##### figure 4
## T1D
total = c(86,80)
success = c(35,30)
total_unique = c(18,12)
unique_success = c(6,1)
zero_column = c("AE-TWAS","Standard TWAS")
first_column = success
second_column = round(success/total,4)*100
third_column = unique_success
fourth_column = round(unique_success/total_unique,4)*100
df_T1D = data.frame(zero_column,first_column,second_column,third_column,fourth_column)
colnames(df_T1D) = c("Method","Success","Success_Rate","Unique_Success","Unique_Success_Rate")
df_T1D$disease = "T1D"

## CD
total = c(16,16)
success = c(5,3)
total_unique = c(5,5)
unique_success = c(4,2)
zero_column = c("AE-TWAS","Standard TWAS")
first_column = success
second_column = round(success/total,4)*100
third_column = unique_success
fourth_column = round(unique_success/total_unique,4)*100
df_CD = data.frame(zero_column,first_column,second_column,third_column,fourth_column)
colnames(df_CD) = c("Method","Success","Success_Rate","Unique_Success","Unique_Success_Rate")
df_CD$disease = "CD"

## RA
total = c(40,38)
success = c(16,20)
total_unique = c(16,14)
unique_success = c(8,12)
zero_column = c("AE-TWAS","Standard TWAS")
first_column = success
second_column = round(success/total,4)*100
third_column = unique_success
fourth_column = round(unique_success/total_unique,4)*100
df_RA = data.frame(zero_column,first_column,second_column,third_column,fourth_column)
colnames(df_RA) = c("Method","Success","Success_Rate","Unique_Success","Unique_Success_Rate")
df_RA$disease = "RA"

## ASD
total = c(673,823)
success = c(23,26)
total_unique = c(255,405)
unique_success = c(15,18)
zero_column = c("AE-TWAS","Standard TWAS")
first_column = success
second_column = round(success/total,4)*100
third_column = unique_success
fourth_column = round(unique_success/total_unique,4)*100
df_ASD = data.frame(zero_column,first_column,second_column,third_column,fourth_column)
colnames(df_ASD) = c("Method","Success","Success_Rate","Unique_Success","Unique_Success_Rate")
df_ASD$disease = "ASD"

## SCZ
total = c(324,332)
success = c(46,40)
total_unique = c(147,155)
unique_success = c(23,17)
zero_column = c("AE-TWAS","Standard TWAS")
first_column = success
second_column = round(success/total,4)*100
third_column = unique_success
fourth_column = round(unique_success/total_unique,4)*100
df_SCZ = data.frame(zero_column,first_column,second_column,third_column,fourth_column)
colnames(df_SCZ) = c("Method","Success","Success_Rate","Unique_Success","Unique_Success_Rate")
df_SCZ$disease = "SCZ"

DisGeNet_combine = rbind(df_T1D,df_CD,df_RA,df_ASD,df_SCZ)

a = ggplot(data=DisGeNet_combine, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                                      y=Success, fill=as.factor(Method))) +
  geom_bar(position="dodge", stat="identity",width=0.5) +
  geom_text(aes(label=Success),position=position_dodge(width=0.5), 
            vjust=-0.5,size=3) +
  scale_fill_manual(values = c("tomato3", "tan")) +
  ylim(0, 50)+labs(x="",y="Validated") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

b = ggplot(data=DisGeNet_combine, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                                      y=Unique_Success, fill=as.factor(Method))) +
  geom_bar(position="dodge", stat="identity",width=0.5) +
  geom_text(aes(label=Unique_Success),position=position_dodge(width=0.5), 
            vjust=-0.5,size=3) +
  scale_fill_manual(values = c("tomato3", "tan")) +
  ylim(0, 50)+labs(x="",y="Unique Validated",fill ="Method")+
  theme_minimal(base_size = 15) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = c(0.8, 0.8),
        legend.title=element_blank(),)

c = ggplot(data=DisGeNet_combine, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                                      y=Success_Rate, fill=as.factor(Method))) +
  geom_bar(position="dodge", stat="identity",width=0.5) +
  geom_text(aes(label=Success_Rate),position=position_dodge(width=0.5), 
            vjust=-0.5,size=3) +
  scale_fill_manual(values = c("tomato3", "tan")) +
  ylim(0,100)+labs(x="Disease",y="Validated Rate(%)") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

d = ggplot(data=DisGeNet_combine, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                                      y=Unique_Success_Rate, fill=as.factor(Method))) +
  geom_bar(position="dodge", stat="identity",width=0.4) +
  geom_text(aes(label=Unique_Success_Rate),position=position_dodge(width=0.5), 
            vjust=-0.5,size=3) +
  scale_fill_manual(values = c("tomato3", "tan")) +
  ylim(0,100)+labs(x="Disease",y="Unique Validated Rate(%)") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

ggarrange(a,b,c,d,labels=c("A","B","C","D"),ncol=2,nrow=2)


##### figure 5
setwd("data/Figures/")
data5a = read.csv("figure5a.csv",header = T)
aggregate(value ~ module + variable, data = data5a, summary) # table11

f5a = ggplot(data=data5a,aes(x=as.factor(module),y=value,fill=variable)) + 
  geom_boxplot() +
  facet_wrap(~module,scale="free") +
  scale_fill_manual(values=c("darkcyan","bisque2")) + 
  theme_minimal(base_size = 10) +
  labs(x="Module",y="Total Connectivity") +
  theme(legend.position = "none",
        strip.text.x.top = element_blank())
f5a

data5b = read.csv("figure5b.csv",header = T)
aggregate(Value ~ Disease + Method, data = data5b, summary) # table12

f5b = ggplot(data=data5b,aes(x=factor(Disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")),y=Value,fill=Method)) +
  geom_boxplot()+
  scale_fill_manual(values=c("tomato3", "tan")) +
  ylim(0,450) + labs(x="",y="")+
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
f5b

data5c = read.csv("figure5c.csv",header = T)
f5c = ggplot(data=data5c, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                                        y=Unique_Success, fill=as.factor(Method))) +
  geom_bar(position="dodge",stat="identity",width=0.5)+
  geom_text(aes(label=Unique_Success),position=position_dodge(width=0.5),
            vjust=-0.5,size=3)+
  scale_fill_manual(values=c("darkcyan","bisque2"))+
  ylim(0,120)+labs(x="",y="Unique Validated")+
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
f5c

f5e = ggplot(data=data5c, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                              y=Unique_Success_Rate*100, fill=as.factor(Method))) +
  geom_bar(position="dodge",stat="identity",width=0.5)+
  geom_text(aes(label=Unique_Success_Rate*100),position=position_dodge(width=0.5),
            vjust=-0.5,size=3)+
  scale_fill_manual(values=c("darkcyan","bisque2"))+
  ylim(0,20)+labs(x="Disesse",y="Unique Validated Rate (%)")+
  theme_minimal(base_size = 10) +
  theme(legend.position = c(0.9, 0.9),
        legend.title=element_blank(),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.4,'cm'))
f5e


data5d = read.csv("figure5d.csv",header = T)
f5d = ggplot(data=data5d, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                              y=Unique_Success, fill=as.factor(Method))) +
  geom_bar(position="dodge",stat="identity",width=0.5)+
  geom_text(aes(label=Unique_Success),position=position_dodge(width=0.5),
            vjust=-0.5,size=3)+
  scale_fill_manual(values=c("tomato3", "tan"))+
  ylim(0,60)+labs(x="",y="")+
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
f5d


f5f = ggplot(data=data5d, aes(x=factor(disease,levels = c("ASD", "SCZ", "T1D", "CD","RA")), 
                              y=Unique_Success_Rate*100, fill=as.factor(Method))) +
  geom_bar(position="dodge",stat="identity",width=0.5)+
  geom_text(aes(label=Unique_Success_Rate*100),position=position_dodge(width=0.5),
            vjust=-0.5,size=3)+
  scale_fill_manual(values=c("tomato3", "tan"))+
  ylim(0,10)+labs(x="Disesse",y="")+
  theme_minimal(base_size = 10) +
  theme(legend.position = c(0.85, 0.9),
        legend.title=element_blank(),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.4,'cm'))
f5f

ggarrange(f5a,f5b,f5c,f5d,f5e,f5f,labels=c("A","B","C","D","E","F"),ncol=2,nrow=3)




