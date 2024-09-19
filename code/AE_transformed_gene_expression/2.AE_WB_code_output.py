#!/bin/python
##Load the trained model to get the code and output
##two inputs: 1st is the name of module (e.g. 13), 2nd is the path of trained module

##Step 1: Import packages
import os
import time
import sys
import pandas as pd
import numpy as np
from torchvision import transforms
import torch
from torch import nn
import torch.nn.functional as F
from torch.autograd import Variable
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pickle
import scipy.stats as stats
from sklearn import preprocessing

##Step 2: load data
M = sys.argv[1]
wk = "/work/long_lab/qli/WGCNA/GTEX/12K_TPM1_1-22/Trial3a_20201220/"
pheno_df = pd.read_csv(wk+"module_train_genes"+M+".csv",sep=",") ##wgcta split 
pheno_df2 = pd.read_csv(wk+"module_test_genes"+M+".csv",sep=",")
pheno_train = pheno_df.iloc[:,1:].to_numpy()
pheno_test = pheno_df2.iloc[:,1:].to_numpy()
scaler2 = preprocessing.MinMaxScaler()
pheno_train_norm = scaler2.fit_transform(pheno_train)
pheno_test_norm = scaler2.transform(pheno_test)

##Step 3: define dataloader
batch_size = np.shape(pheno_train_norm)[0]
test_batch_size = np.shape(pheno_test_norm)[0]

class InputDataset(Dataset):
    def __init__(self, np_array):
        #'Initialization'
        self.x = np_array
        self.n_sample = np.shape(np_array)[0]
    def __len__(self):
        #'Denotes the total number of samples'
        return self.n_sample
    def __getitem__(self, index):
        #'Generates one sample of data'
        return self.x[index]
        
pheno_train_set = InputDataset(pheno_train_norm)
pheno_train_set_loader = DataLoader(pheno_train_set,batch_size=batch_size,shuffle=True,num_workers=8)
pheno_test_set = InputDataset(pheno_test_norm)
pheno_test_set_loader = DataLoader(pheno_test_set,batch_size=test_batch_size,shuffle=False,num_workers=8)
        
        
##Step 4: define model
input_features=pheno_train_norm.shape[1]
output_features=input_features

h1 = int(input_features/8)

class Auto_Pheno(nn.Module):
    def __init__(self):
        super(Auto_Pheno,self).__init__() 
        # def the encoder function
        self.encoder = nn.Sequential(
            nn.Linear(input_features, h1),
            nn.Sigmoid()
        )
        #def the decoder function
        self.decoder = nn.Sequential(
            nn.Linear(h1,output_features),
            nn.Sigmoid()
        )       
    #def forward function
    def forward(self,x):
        y = self.encoder(x)
        x = self.decoder(y)
        return x,y
        
##Step 5: reload trained model
PATH=sys.argv[2]
model = Auto_Pheno()
# Load the model onto CPU
model.load_state_dict(torch.load(PATH, map_location=torch.device('cpu')))
# Ensure the model is on the CPU
model = model.to(torch.device('cpu'))

##Step 5.1: generate a final phenotype file
pheno_df_all = pd.concat([pheno_df,pheno_df2],axis=0)
pheno_df_all_np_before_norm = pheno_df_all.iloc[:,1:].to_numpy()
pheno_df_all_np = scaler2.fit_transform(pheno_df_all_np_before_norm)

##Step 5.2: setup a data loader for this final phenotype file
batch_size = np.shape(pheno_df_all_np)[0]
pheno_set = InputDataset(pheno_df_all_np)
pheno_set_loader = DataLoader(pheno_set,batch_size=batch_size,shuffle=False,num_workers=8)

##Step 6: apply whole model to the dataset
model.eval()
distance = nn.MSELoss() # for regression

output_pheno = np.zeros_like(pheno_df_all_np)
input_pheno = pheno_df_all_np
coder2 = np.zeros([pheno_df_all_np.shape[0],h1])
test_sum_loss = 0

for batch_count, pheno_data in enumerate(pheno_set_loader):
    pheno = Variable(pheno_data).float()  # Put test_pheno into CPU
    #=================forward==============
    output, coder = model.forward(pheno)
    loss = distance(output, pheno)
    test_sum_loss += loss.item()
    #===============precision============
    output2 = output.detach().numpy()  # Move tensor to numpy
    start_ind = batch_count * batch_size
    end_ind = batch_count * batch_size + output.shape[0]
    output_pheno[start_ind:end_ind] = output2
    input_pheno[start_ind:end_ind] = pheno_data.numpy()
    coder2[start_ind:end_ind] = coder.detach().numpy()
    
#================loss+average cor of (y,y_hat)======================
print('Loss: {:.4f}'.format(test_sum_loss))
res_r2_np_mean = r2_score(input_pheno, output_pheno)  # r2_score(y_true, y_pred)
print('The average R^2 between y and y_hat for final phenotypes is: {:.4f}'.format(res_r2_np_mean))

## Step 7: Output hidden code and X'
df1 = pheno_df_all.iloc[:, 0]
df1.columns = ['IID']
df1.index = pheno_df_all.iloc[:, 0]

output_pheno_df = pd.DataFrame(output_pheno)
output_pheno_df.columns = pheno_df_all.columns[1:]
output_pheno_df.index = pheno_df_all.iloc[:, 0]
out_pheno_df_res = pd.concat([df1, output_pheno_df], axis=1)

coder2_df = pd.DataFrame(coder2)
coder2_df.index = pheno_df_all.iloc[:, 0]
coder2_df_res = pd.concat([df1, coder2_df], axis=1)

R_2=round(res_r2_np_mean,2)
out_pheno_df_res.to_csv("/work/long_lab/qli/AE-TWAS/GTEX_AE/revision202408/gene_expression_X_h_X_apo/AE_addBatchNoise/"+M+"_Div8_ID670_R"+str(R_2)+"_output_gcta.txt",sep=" ",header=None,index=True)
out_pheno_df_res.to_csv("/work/long_lab/qli/AE-TWAS/GTEX_AE/revision202408/gene_expression_X_h_X_apo/AE_addBatchNoise/"+M+"_Div8_ID670_R"+str(R_2)+"_output_head_gcta.txt",sep=" ",header=True,index=True)

coder2_df_res.to_csv("/work/long_lab/qli/AE-TWAS/GTEX_AE/revision202408/gene_expression_X_h_X_apo/AE_addBatchNoise/"+M+"_Div8_ID670_R"+str(R_2)+"_code_gcta.txt",sep=" ",header=None,index=True)
coder2_df_res.to_csv("/work/long_lab/qli/AE-TWAS/GTEX_AE/revision202408/gene_expression_X_h_X_apo/AE_addBatchNoise/"+M+"_Div8_ID670_R"+str(R_2)+"_code_head_gcta.txt",sep=" ",header=True,index=True)

pheno_df_all.index=pheno_df_all.iloc[:,0]
pheno_df_all.to_csv("/work/long_lab/qli/AE-TWAS/GTEX_AE/revision202408/gene_expression_X_h_X_apo/AE_addBatchNoise/"+M+"_Div8_ID670_R"+str(R_2)+"_input_gcta.txt",sep=" ",header=None,index=True)
pheno_df_all.to_csv("/work/long_lab/qli/AE-TWAS/GTEX_AE/revision202408/gene_expression_X_h_X_apo/AE_addBatchNoise/"+M+"_Div8_ID670_R"+str(R_2)+"_input_head_gcta.txt",sep=" ",header=True,index=True)
