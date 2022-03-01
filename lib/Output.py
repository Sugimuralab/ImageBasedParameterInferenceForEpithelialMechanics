# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 16:56:59 2022

@author: Goshi
"""
import os


import numpy as np
import matplotlib.pyplot as plt

fontsize_label = 24
fontsize_axis = 15

def MakeOutputDirectory(SampleName,Output=""):
    """
    結果を出力するディレクトリを作成する
    Make the directry for output results
    """
    if Output == "":
        cwd = os.getcwd()
        Output = os.path.join(cwd,"output"+os.sep)
    OutputSample = os.path.join(Output,SampleName+os.sep)
    os.makedirs(OutputSample,exist_ok=True)
    return [Output, OutputSample]        


def R_alpha(x1, x2):
        """
        角度合成する
        Composit x1*sin(phi) + x2*cos(phi) -> r*cos(phi - alpha)
        
        """
        r = np.linalg.norm([x1,x2])
        alpha = np.arctan2(x1, x2)
        
        if alpha.values<0:
            alpha = alpha/2 + np.pi
        else:
            alpha = alpha/2
        
        return [r, alpha]

def CalcRT(theta,T):
    """
    張力異方性の強さと向きを計算する
    Calculate the strength and oriention of tension anisotropy
    """
    COS = np.cos(2 * theta)
    SIN = np.sin(2 * theta)
    SNmean = np.mean([T * COS, T * SIN], axis = 1)
    Smean = np.mean(T) 
    Nmean = np.mean([COS, SIN], axis = 1)
    ANISO = SNmean - Smean * Nmean
    THETA = np.arctan2(ANISO[1],ANISO[0])/2
    if THETA < 0:
        THETA = THETA + np.pi
    RT = np.linalg.norm(ANISO)/Smean
    return [RT, THETA]

def ConvertAnisotropicParameters(params):
    ar,phi_a = R_alpha(params.a1,params.a2)
    br,phi_b = R_alpha(params.b1/params.b0,params.b2/params.b0)
    return [ar,phi_a, br, phi_b]

#%%推定パラメータの出力
def SaveParams(Samplename,BestModel,params,minAIC,SampleName,Output):
    """
    csvファイルとして推定パラメータを保存する
    Save the estimated parameters to a csv file
    """
    [ar,phi_a,br,phi_b] = ConvertAnisotropicParameters(params)    
    params["ar"] = ar
    params["phi_a"] = phi_a
    params["br"] = br
    params["phi_b"] = phi_b
    #params["AIC"] = minAIC
    #params["Model"] = BestModel
    params = params.reindex(columns=["Model","AIC","a0","ar","phi_a","b0","br","phi_b","k"])
    params.to_csv(os.path.join(Output ,"EstParams.csv") ,header=True, mode="a")

def PlotTL(l, theta, T, SampleName, OutputSample):
    """
    細胞接着面の張力と長さ/向きの相関をプロットする
    Plot correlation between tension and length/orientation of cell junctions
    """
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    
    marker = list()
    for ti in theta:
        if ti <= np.pi/8:
            marker.append('.r')
        elif ti >= 7*np.pi/8:
            marker.append('.r')
        elif np.pi/8 < ti <= 3*np.pi/8:
            marker.append('.g')
        elif 3*np.pi/8 < ti <= 5*np.pi/8:
            marker.append('.b')
        else:
            marker.append('.m')
    
    for i in range(len(l)):
        plt.plot(l[i],T[i],marker[i])

    #Axis
    plt.ylim([0,1.2])
    plt.xlim([0,1.75])
    plt.xticks([.50,1.0,1.5],fontsize=fontsize_axis)
    plt.yticks(fontsize=fontsize_axis)

    plt.xlabel("Edge length / $\\sqrt{<A>}$",fontsize=fontsize_label)
    plt.ylabel("Predicted tension\n[a.u.]",fontsize=fontsize_label)

    plt.tight_layout()
    plt.savefig(OutputSample+SampleName+"_lT.png",dpi=400)
    plt.close()

def PlotPA(A,P,SampleName, OutputSample):
    """
    細胞の圧力と面積の相関をプロットする
    Plot correlation between pressure and area of cells
    """
    
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')

    plt.plot(A,P,".k")
    plt.xlabel("Cell area / $<A>$",fontsize=fontsize_label)
    plt.ylabel("Predicted $\Delta$Pressure\n [a.u.]",fontsize=fontsize_label)
    plt.xticks([.5,1.0,1.5],fontsize=fontsize_axis)
    plt.yticks(fontsize=fontsize_axis)
    plt.tight_layout()


    plt.savefig(OutputSample+SampleName+"_AP.pdf")
    plt.close()
