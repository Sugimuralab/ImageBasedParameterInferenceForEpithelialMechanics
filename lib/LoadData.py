#!/usr/bin/env python
# coding:UTF-8

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import   Polygon
from matplotlib.collections import PatchCollection
import scipy.sparse as sp


'''
   Classes
     EDGE
     CELL
'''

class EDGE:
    def __init__(self):
        self.junc1 = -1
        self.junc2 = -1
        self.x1 = 0
        self.y1 = 0
        self.x2 = 0
        self.y2 = 0
        self.in_out = 'N'
        self.ncell = [-1, -1]
        self.dx = 0.0
        self.dy = 0.0
        self.dist = 0.0
        self.degree = 0.0
        self.E_peri = 0

    def set_distance(self,x,y):
        self.dx = self.x1 - self.x2
        self.dy = self.y1 - self.y2
        self.dist = np.sqrt( self.dx*self.dx + self.dy*self.dy )
        #remove abs
        self.degree = np.arctan2(self.dy,self.dx)

class CELL:
    def __init__(self):
        self.jnum = 0
        self.junc = np.empty(0,dtype = np.int32)
        self.in_out = 'N'
        self.edge = []
        self.NON_EDGE = np.nan
        self.area = 0

    def initialize_edge(self, NON_EDGE = 10000000000):
        self.NON_EDGE = NON_EDGE
        self.edge = np.ones(self.jnum, dtype = int)*NON_EDGE
        self.area = 0
'''
   Functions
     Set_cellNeighbors(tedge, tcell, E_NUM)  ## set edge.ncell and cell.edge
     loaddata(filename)
     GetMatrix_ForceEstimation(x,y,edge,cell,E_NUM,CELL_NUMBER,R_NUM,INV_NUM,Rnd,ERR_MAX):
     Show_L_curve(MM,B,G,E_NUM,CELL_NUMBER,mus):
'''


def Set_cellNeighbors(tedge, tcell, non_edge):
    """
    細胞の隣接関係を格納する
    Store the cell adjacencies
    """
    ejs = np.vstack( ([e.junc1 for e in tedge],[e.junc2 for e in tedge]) ).T
    cis = np.empty( (4,0), int)
    for (cc, c) in enumerate(tcell):
        c.initialize_edge( NON_EDGE = non_edge )
        cjuncs = np.vstack( (c.junc, np.roll(c.junc,-1)) )
        cjuncs = np.vstack( (cjuncs, cc*np.ones(c.jnum, int)) )
        cjuncs = np.vstack( (cjuncs, np.arange(c.jnum)) )
        cis = np.hstack( (cis,cjuncs) )

    for (e, ej) in enumerate(ejs):
        ne1 = np.intersect1d( np.where( ej[0] == cis[0]), np.where( ej[1] == cis[1]) )
        ne2 = np.intersect1d( np.where( ej[1] == cis[0]), np.where( ej[0] == cis[1]) )
        if (len(ne1)!=1 and len(ne2)!=1): assert('inconsistent')
        tedge[e].ncell[0] = int(cis[2,ne1])
        tedge[e].ncell[1] = int(cis[2,ne2])
        tcell[int(cis[2,ne1])].edge[int(cis[3,ne1])] =  e
        tcell[int(cis[2,ne2])].edge[int(cis[3,ne2])] = -e



def loaddata(filename,CHECK = False):
    """
    データを読み込み、変数に格納する
    Load data and store its information to variables 
    """
    print('... Load data from', filename)
    stl = 0
    R_NUM=0
    CELL_NUMBER = 0
    E_NUM = 0
    V_NUM = 0
    INV_NUM = 0
    RndV = np.empty(0, dtype = np.int32)
    RndE = np.empty(0, dtype = np.int32)
    RndC = np.empty(0, dtype = np.int32)

    title = ''

    for line in open(filename, 'r'):

        if '#' in line:
            if line[0:2] == '# ' and title == '':
                itemList = line[:-1].split()
                title = itemList[-1]
                print('title = \"',title,'\"')

            elif 'C_NUM' in line:
                itemList = line[:-1].split()
                CELL_NUMBER = int(itemList[2])
                print('CELL_NUMBER = ',  CELL_NUMBER)
                cell = [ CELL() for i in range(CELL_NUMBER)]  # Declare cell here

            elif 'IN_CNUM' in line:
                itemList = line[:-1].split()
                IN_CNUM = int(itemList[2])
                print('IN_CNUM = ', IN_CNUM )

            elif 'EX_CNUM' in line:
                itemList = line[:-1].split()
                R_NUM = int(itemList[2])
                print('R_NUM = ', R_NUM )

            elif '# V_NUM' in line:
                itemList = line[:-1].split()
                V_NUM = int(itemList[2])
                INV_NUM=V_NUM-R_NUM
                # N = 2*V_NUM   # !! Notice: Not 2*(V_NUM+R_NUM) !!
                x = np.zeros(V_NUM)  # Declare x
                y = np.zeros(V_NUM)  # Decalre y

            elif '# E_NUM' in line:
                itemList = line[:-1].split()
                E_NUM = int(itemList[2])
                print('E_NUM = ', E_NUM)
                edge = [ EDGE() for i in range(E_NUM)]  # Declare edge here

            elif '### Scale' in line:
                itemList = line[:-1].split()
                stl = int(itemList[2])
                print('scale' .stl)

        else:
            if 'V[' in line:
                itemList = line[:-1].split()
                vid = int(itemList[0].replace('V[','').replace(']',''))
                if vid >= len(x) : assert('')
                x[vid] = np.float64(itemList[1])
                y[vid] = np.float64(itemList[2])
                if 'Ext' in line:   # IDs of verteces at boundary
                    rid = int(itemList[0].replace('V[','').replace(']','') )
                    RndV = np.append(RndV, rid )

            elif 'E[' in line:
                #print(line)
                itemList = line[:-1].split()
                eid = int(itemList[0].replace('E[','').replace(']',''))
                if eid >= len(edge):   assert('edge Number is larger than expected')
                edge[eid].junc1 = int( itemList[1] )
                edge[eid].junc2 = int( itemList[2] )
                edge[eid].x1 = x[edge[eid].junc1]
                edge[eid].y1 = y[edge[eid].junc1]
                edge[eid].x2 = x[edge[eid].junc2]
                edge[eid].y2 = y[edge[eid].junc2]
                edge[eid].set_distance(x,y)

                if 'Ext' in line:  # IDs of edges at boundary
                    eid = int(itemList[0].replace('E[','').replace(']','') )
                    RndE = np.append(RndE, eid )
                    edge[eid].in_out = 'o'
                else:
                    edge[eid].in_out = 'i'


            elif 'C[' in line:
                itemList = line[:-1].split()
                cid = int(itemList[0].replace('C[','').replace(']',''))
                if cid >= len(cell):  assert('cell Number is larger than expected')
                cell[cid].jnum = int(itemList[1])
                for i in range(cell[cid].jnum):
                    cell[cid].junc = np.append( cell[cid].junc, np.int32(itemList[i+3]) )

                if 'Ext' in line:  # IDs of cells at boundary
                    cid = int(itemList[0].replace('C[','').replace(']','') )
                    RndC = np.append(RndC, cid )
                    cell[cid].in_out = 'o'
                else:
                    cell[cid].in_out = 'i'

    ### end of for line

    #    Rnd = np.append( RndJ, RndE,  RndC)

    Rnd = (RndV, RndE,RndC)
    #print(len(RndJ),RndJ)
    #print(len(RndE),RndE)
    #print(len(RndC),RndC)

    ### cell area : For checking validity of data
    ff=0
    for c in cell:
        area = 0.0
        darea = 0.0
        xx = x[c.junc]
        sx = np.roll(xx, -1)
        yy = y[c.junc]
        sy = np.roll(yy, -1)
        area = 0.5*sum(xx * sy -yy * sx)
        c.area = area
        if area<= 0.0:
            ff=1
            print('!! error: %c %d %f\n' & (c.in_out,i,area) )

    if ff==1:
        assert('')


    ###  set neibouring cells of edge: edges of cells %
    if CHECK:
        Set_cellNeighbors(edge, cell,E_NUM)

        ### check consistency of edge-cell relation 1
        ct = np.zeros(CELL_NUMBER, int);
        for e in edge:
            if e.ncell[0] == -1 or  e.ncell[1] == -1:
                assert('edge-cell error')
            ct[e.ncell[0]] += 1
            ct[e.ncell[1]] += 1

        for i,c in enumerate(cell):
            if c.in_out == 'o' and ct[i] != c.jnum-1:
                print('cell = %d count inconsistent  %d != %d ' % (i,ct[i],c.jnum-1) )
            elif c.in_out == 'i' and ct[i] !=  c.jnum:
                print('cell = %d count inconsistent  %d != %d ' % (i,ct[i],c.jnum) )


    ## Return data
    return [x,y,edge,cell,Rnd,CELL_NUMBER,E_NUM,V_NUM,INV_NUM,R_NUM,stl,title]

###############################################################


def DrawCells(x,y,edge,cell,OutputSample):
    """
    入力データの細胞配置を描画
    Draw the cell geometry of an input data
    """
    print('...    Show cells               ')
    cfig, ax = plt.subplots()
    for e in edge:
        tx = [x[e.junc1], x[e.junc2]]
        ty = [y[e.junc1], y[e.junc2]]
        plt.plot( tx, ty, '-k' )
    ax.set_aspect('equal', 'datalim')
    cfig.savefig(OutputSample+"CheckCells.png")
    


