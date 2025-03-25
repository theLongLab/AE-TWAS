#!/bin/python

##Step 1: load packages
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
from datetime import datetime
date = datetime.now().date().strftime("%Y-%m-%d")

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
h1 = int(input_features/2)
h2 = int(input_features/4)
h3 = int(input_features/8)
class Auto_Pheno_shallow(nn.Module):
    def __init__(self):
        super(Auto_Pheno_shallow,self).__init__() 
        # def the encoder function
        self.encoder = nn.Sequential(
            nn.Linear(input_features, h1),
            nn.Sigmoid(),
            nn.Linear(h1, h2),
            nn.Sigmoid(),
            nn.Linear(h2, h3),
            nn.Sigmoid()
        )
        #def the decoder function
        self.decoder = nn.Sequential(
            nn.Linear(h3, h2),
            nn.Sigmoid(),
            nn.Linear(h2,h1),
            nn.Sigmoid(),
            nn.Linear(h1,output_features),
            nn.Sigmoid()
        )       
    #def forward function
    def forward(self,x):
        y = self.encoder(x)
        x = self.decoder(y)
        return x,y
        
##Step 5: train and test model
model = Auto_Pheno_shallow().cuda()
distance = nn.MSELoss() # for regression
optimizer = torch.optim.Adam(model.parameters(),lr=0.001) 
# amp.initialize(model, optimizer, opt_level="O1")
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=2000, gamma=0.5)
do_train = True
do_test = True

import time
start = time.time()
res_r2= []

num_epochs = 20000
train_test_log = open(wk+"/"+date+"_WB_"+M+".log",'w')

for epoch in range(num_epochs):
    if do_train:
        sum_loss = 0
        model.train()
        output_pheno_train = np.zeros_like(pheno_train_norm)
        input_pheno_train = np.zeros_like(pheno_train_norm)
        coder2 = np.zeros([pheno_train.shape[0],h3])
        for batch_count, pheno_data in enumerate(pheno_train_set_loader):
            train_pheno = Variable(pheno_data).float().cuda() # put train_pheno into GPU
            #=================forward==============
            output,coder = model.forward(train_pheno)
            loss = distance(output,train_pheno)
            sum_loss += loss.item()
            #===============precision============
            output2 = output.cpu().detach().numpy()
            start_ind = batch_count * batch_size
            end_ind = batch_count * batch_size + output.shape[0]
            output_pheno_train[start_ind:end_ind] = output2
            input_pheno_train[start_ind:end_ind] = pheno_data.cpu().numpy() 
            coder2[start_ind:end_ind] = coder.cpu().detach().numpy()
            #=================backward=============
            optimizer.zero_grad() # initial gradients
            loss.backward() #differential  
            optimizer.step() #update weights 
        scheduler.step()
        #================loss+average cor of (y,y_hat)======================
        if epoch % 100 ==0:
            train_test_log.write('LR: {}\n'.format(scheduler.get_lr()))
            train_test_log.write('epoch[{}/{}],loss:{:.4f}\n'.format(epoch+1,num_epochs,sum_loss))
            res_r2_np_mean = r2_score(input_pheno_train,output_pheno_train) #r2_score(y_true, y_pred)
            train_test_log.write('The average R^2 between y and y_hat for TRAIN phenotypes is: {:.4f}\n'.format(res_r2_np_mean))
            print("LR: ",scheduler.get_lr())
            print('epoch[{}/{}],loss:{:.4f}'.format(epoch+1,num_epochs,sum_loss))
#             res_r2_np_mean = r2_score(input_pheno_train,output_pheno_train) #r2_score(y_true, y_pred)
            print('The average R^2 between y and y_hat for TRAIN phenotypes is: {:.4f}'.format(res_r2_np_mean))
    if do_test:
        test_sum_loss = 0
        model.eval()
        output_pheno_test = np.zeros_like(pheno_test_norm)
        input_pheno_test = pheno_test_norm
        coder2 = np.zeros([pheno_test.shape[0],h3])
        for batch_count, pheno_data in enumerate(pheno_test_set_loader):
            test_pheno = Variable(pheno_data).float().cuda() # put test_pheno into GPU
            #=================forward==============
            output,coder = model.forward(test_pheno)
            loss = distance(output,test_pheno)
            test_sum_loss += loss.item()
            #===============precision============
            output2 = output.cpu().detach().numpy()
            start_ind = batch_count * batch_size
            end_ind = batch_count * batch_size + output.shape[0]
            output_pheno_test[start_ind:end_ind] = output2
            input_pheno_test[start_ind:end_ind] = pheno_data.cpu().numpy() 
            coder2[start_ind:end_ind] = coder.cpu().detach().numpy()
        #================loss+average cor of (y,y_hat)======================
        if epoch % 100 ==0:
            train_test_log.write('LR: {}\n'.format(scheduler.get_lr()))
            train_test_log.write('epoch[{}/{}],loss:{:.4f}\n'.format(epoch+1,num_epochs,test_sum_loss))
            res_r2_np_mean = r2_score(input_pheno_test,output_pheno_test) #r2_score(y_true, y_pred)
            train_test_log.write('The average R^2 between y and y_hat for TEST phenotypes is: {:.4f}\n'.format(res_r2_np_mean))         
            print("LR: ",scheduler.get_lr())
            print('epoch[{}/{}],loss:{:.4f}'.format(epoch+1,num_epochs,test_sum_loss))
#             res_r2_np_mean = r2_score(input_pheno_test,output_pheno_test) #r2_score(y_true, y_pred)
            print('The average R^2 between y and y_hat for TEST phenotypes is: {:.4f}'.format(res_r2_np_mean))
train_test_log.close()

##Step 6: save model
torch.save(model.state_dict(), wk+M)
