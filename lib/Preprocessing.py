# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 16:44:12 2022

@author: Goshi
"""

#%% Import modules
import numpy as np
import copy as cp

import lib.ParameterInf_lib as PE

def scale_converter(x,y,edge,cell,sc_mean = 1.0,sc = 0,flip=""):
    """
    空間スケールを変換する
    Convert the spatial scale
    """
    E_NUM = len(edge)
    CELL_NUMBER = len(cell)
    
    dist_error = 1.0e-10
    
    if(sc == 0):
        lmean = np.mean([edge[i].dist for i in range(E_NUM)])
        sc = sc_mean/lmean
    if(flip == "v"):
        xsc = -sc*x
        ysc = sc*y
    elif(flip == "h"):
        xsc = sc*x
        ysc = -sc*y
    else:
        xsc = sc*x
        ysc = sc*y
    
    edgesc = cp.deepcopy(edge)
    for i in range(E_NUM):
        edgesc[i].x1 = xsc[edgesc[i].junc1]
        edgesc[i].y1 = ysc[edgesc[i].junc1]
        edgesc[i].x2 = xsc[edgesc[i].junc2]
        edgesc[i].y2 = ysc[edgesc[i].junc2]
        edgesc[i].set_distance(x,y)
        if( edgesc[i].dist - sc*edge[i].dist > dist_error):
            print("dist_error = %f > %f" %(edge[i].dist - sc*edge[i].dist,dist_error))
        
    cellsc = cp.deepcopy(cell)
    for i in range(CELL_NUMBER):
        cellsc[i].area = sc*sc*cell[i].area
        
    return [xsc,ysc,edgesc,cellsc,sc]

def OutlierDetector(ListWOutlier,maximum_mag):
    """
    外れ値（maximum_mag * median < ）を検出する
    Detect the outlier (maximum_mag * median <)
    """
    MAD = np.median(ListWOutlier)
    Upper = (maximum_mag) * MAD < ListWOutlier
    print("Total: {}".format(len(ListWOutlier)))
    print("Upper Outleir: {}".format(sum(Upper)))
    return  np.where( Upper  )

def Preprocessing(Data,ExcludeOutlier, ExcludeShortEdge, AreaNormalization, lmin=3, area_threshold_factor=2):
    """
    前処理を行う：短い細胞接着面の除外・大きい細胞の除外・空間スケールの変換
    Preprocess the input data: exlude short junctions and large cells, rescale the spatial scale
    """
    [x,y,edge,cell,Rnd,CELL_NUMBER,E_NUM,V_NUM,INV_NUM,R_NUM,stl,title] = Data
    [V_in, E_in, C_in, RndV, RndE, RndC] = PE.GetInsiceVEC(edge, cell, V_NUM, CELL_NUMBER, E_NUM, INV_NUM, Rnd)


    if ExcludeShortEdge:
        short_edges = [i for i in E_in if edge[i].dist <lmin]
        E_ex = short_edges

    if ExcludeOutlier:
        OutlierC_in = OutlierDetector([cell[i].area for i in C_in], area_threshold_factor)
        OutlierVertices = list({ji for ci in C_in[OutlierC_in] for ji in cell[ci].junc} )
        OutlierV_in = np.where([(ji in OutlierVertices) for ji in V_in])[0]
        RndV = np.append(RndV,OutlierV_in)

    if AreaNormalization:
        [V_in, E_in, C_in, RndV, RndE, RndC] = PE.GetInsiceVEC(edge, cell, V_NUM, CELL_NUMBER, E_NUM, INV_NUM, [RndV,RndE,RndC])
        non_dim = 1/np.sqrt(np.median([cell[i].area for i in C_in]))
        [x,y,edge,cell,_] = scale_converter(x,y,edge,cell, sc=non_dim)
    return [x,y,edge,cell,Rnd,CELL_NUMBER,E_NUM,V_NUM,INV_NUM,R_NUM,stl,title,E_ex]
