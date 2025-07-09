# -*- coding: utf-8 -*

"""

-Estimate parameters for candidate mechanical models in epithelial tissue
-Select the most predictive model from the candidate models

As presented in
Goshi Ogita, Takefumi Kondo, Keisuke Ikawa, Tadashi Uemura, Shuji Ishihara & Kaoru Sugimura (2022):
"Image based parameter inference for epithelial mechanics"

Parameter names are modified in this code for readability.

Paper  λ0 λ0 λ0 λ1 λ1 λ1 μ0  μ1  φ0    φ1    k
Code   a0 a1 a2 b0 b1 b2 a_r b_r phi_a phi_b k

 
"""
#%% import modules
import numpy as np
from pandas import DataFrame
import statsmodels.api as sm
import sys
import os

import lib.LoadData as LD
import lib.ParameterInf_lib as PE
import lib.Preprocessing as PP
import lib.Output as OP


#%%　select input file
filename = '.{0}Sample{0}sample{0}Vertex{0}VDat_sample.dat'.format(os.sep)
SampleName = filename.split(os.sep)[-3]

#%% make directory to save the estimate results 
[Output, OutputSample] = OP.MakeOutputDirectory(SampleName)

#%% setting for preprocessing
ExcludeOutlier = True
ExcludeShortEdge = True
AreaNormalization = True
lmin = 3 #If you apply this method to simulation data, you should tune this. (Equal or little bit more than the threshold length for cell rearrangement is fine?)
area_threshold_factor = 2

#%% construct mechanical model as a matrix
ParameterNames = ["a0","a1","a2","b0","b1","b2","k"]
def calc_L(l,theta,A):
        L = np.zeros([len(l) + len(A), len(ParameterNames)]) 
        L[:len(l),0]  = 1
        L[:len(l),1] = np.sin(2 * theta)
        L[:len(l),2] = np.cos(2 * theta)
        L[:len(l),3] = l
        L[:len(l),4] = l * np.sin(2 * theta)
        L[:len(l),5] = l * np.cos(2 * theta)
        L[len(l):,6] = -A
        return L 

#%% list of candidate models
Models = ["A", "B", "C", "D", "E"]
Model2SubjComb = {"A": [0,1,2,3,4,5],
                  "B": [0,1,2,5],
                  "C": [0,1,5],
                  "D": [2,5],
                  "E": [5]}

#%% load data
Data = LD.loaddata( filename )

#%% preprocessing
[x,y,edge,cell,Rnd,CELL_NUMBER,E_NUM,V_NUM,INV_NUM,R_NUM,stl,title,E_ex] = \
    PP.Preprocessing(Data, ExcludeOutlier, ExcludeShortEdge, AreaNormalization,lmin,)

#%% diplay tisse geometry
LD.DrawCells(x, y, edge, cell,OutputSample)

#%% calculate the coefficient matrix C
[C,C_NUM,X_NUM,[V_in,E_in,C_in],Rnd,INV_NUM] = \
    PE.CalcMatrixC(x,y,edge,cell,E_NUM,CELL_NUMBER,R_NUM,INV_NUM,Rnd,ERR_MAX=2.0e-12)

#%% calculate the cell shape features
l = np.array([edge[i].dist for i in E_in])
theta = np.array([edge[i].degree for i in E_in])          
A = np.array([cell[i].area for i in C_in])

#%% calculate cell shape feature matrix L 
L = calc_L(l,theta,A)

#%% devide matrix
CL = np.dot(C,L)    
X = CL[:,1:CL.shape[1]]
Y = -CL[:,0]
    
#%% make lists to save the estimation results of candidate models
ParamList = np.full((len(Models),X.shape[1]+1),np.nan)
ParamList[:,0] = 1
AIC = np.full(len(Models),np.nan)
TList = np.zeros((len(E_in),len(Models)))
PList = np.zeros((len(C_in),len(Models)))
#%% fitting
for mi,model in enumerate(Models): 
    # pick up columns of matrix X corresponding to a candidate model
    SubjComb = Model2SubjComb[model]
    X_mi = X[:,SubjComb]
    L_mi = L[:,np.append([0],np.array(SubjComb)+1)]
    
    # fitting
    res = sm.OLS(Y,X_mi).fit()
    ParamList[mi,np.array(SubjComb)+1] = res.params
    # obtain aic
    AIC[mi] = res.aic
    
    # calculate tension and pressure
    F = L_mi.dot(np.append([1],res.params))
    TList[:,mi] = F[:len(E_in)]
    P = F[len(E_in):]
    Pmean = np.mean(P)
    PList[:,mi] = P - Pmean

#%% model selection
BestModelID = np.argmin(AIC)
minAIC = AIC[BestModelID]
BestModel = Models[BestModelID]
params = DataFrame([ParamList[BestModelID]],columns=ParameterNames,index=[SampleName])
T = TList[:,BestModelID]
P = PList[:,BestModelID]

#%% calculate tension aniosotropy R_T and THETA
[RT,THETA] = OP.CalcRT(theta,T)

#%% save estimated parameters as csv file
OP.SaveParams(SampleName, BestModel, params, minAIC,SampleName,Output)
#%% plot correlation between tension and length/orientation of cell junctions
OP.PlotTL(l,theta,T,SampleName,OutputSample)
 #%% plot correlation between pressure and area of cells
OP.PlotPA(A,P,SampleName,OutputSample) 

## end of file.
