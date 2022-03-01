import numpy as np
import sys
import scipy.sparse as sp


def GetInsiceVEC(edge,cell,V_NUM,CELL_NUMBER,E_NUM,INV_NUM,Rnd):
    """
    細胞頂点・細胞接着面・細胞のうち、推定に使うもの（V_in, E_in, C_in）と使わないもの（RndV, RndE, RndC）を指定
    Designate cell vertices, cell junctions, cells used (V_in, E_in, C_in) and
    and not used (RndV, RndE, RndC) for parameter estimation
    """
    RndV = np.array(Rnd[0],np.int)
    RndE = np.array(Rnd[1],np.int)
    RndC = np.array(Rnd[2],np.int)
    V_NUM = INV_NUM + len(RndV)        
    
    #Cin : Cells inside tissue
    C_in = set(range(CELL_NUMBER))-set(RndC)
    #V_in : Vertices which are not included with Border Cells(RndC) 
    RndV = set(RndV)
    for i in RndC:
        RndV = RndV|set(cell[i].junc)
    V_in = set(range(V_NUM)) - RndV
    #Exclude edges not to use
        
    #E_in : Edges which are not included with Border Cells(RndC)
    for i in range(E_NUM):
        if edge[i].junc1 in RndV and edge[i].junc2 in RndV:
            RndE = np.append(RndE,i)
    RndE = set(RndE)
    E_in = set(range(E_NUM)) - RndE
    #Convert set to np.array
    V_in = np.array(list(V_in),np.int)
    RndV = np.array(list(RndV),np.int)
    if (len(V_in)+len(RndV) != V_NUM):
        print("V_in : %d + RndJ : %d != V_NUM %d, inconsistent!" %(len(V_in),len(RndV),V_NUM))
    E_in = np.array(list(E_in),np.int)
    RndE = np.array(list(RndE),np.int)
    C_in = np.array(list(C_in),np.int)
    RndC = np.array(list(RndC),np.int)
    return [V_in, E_in, C_in, RndV, RndE, RndC]

def CalcMatrixC(x,y,edge,cell,E_NUM,CELL_NUMBER,R_NUM,INV_NUM,Rnd,ERR_MAX,SPARSE = False,E_ex = []):
    """
    係数行列Cを計算する
    Calculate the coefficient matrix C
    
    """
    print('...    Calculate Coefficeint Matrix        ')
    RndJ = np.array(Rnd[0],np.int)
    RndE = np.array(Rnd[1],np.int)
    RndC = np.array(Rnd[2],np.int)
    V_NUM = INV_NUM + len(RndJ)        
    
    #Cin : Cells inside tissue
    C_in = set(range(CELL_NUMBER))-set(RndC)
    #V_in : Vertices which are not included with Border Cells(RndC) 
    RndJ = set(RndJ)
    for i in RndC:
        RndJ = RndJ|set(cell[i].junc)
    V_in = set(range(V_NUM)) - RndJ
    #Exclude edges not to use
    if(len(E_ex) != 0):        
        for i in E_ex:            
            if(edge[i].junc1 in V_in):
                V_in.remove(edge[i].junc1)
                RndJ.add(edge[i].junc1)
            if(edge[i].junc2 in V_in):
                V_in.remove(edge[i].junc2)
                RndJ.add(edge[i].junc2)
        
    #E_in : Edges which are not included with Border Cells(RndC)
    for i in range(E_NUM):
        if edge[i].junc1 in RndJ and edge[i].junc2 in RndJ:
            RndE = np.append(RndE,i)
    RndE = set(RndE)
    E_in = set(range(E_NUM)) - RndE
    #Convert set to np.array
    V_in = np.array(list(V_in),np.int)
    RndJ = np.array(list(RndJ),np.int)
    if (len(V_in)+len(RndJ) != V_NUM):
        print("V_in : %d + RndJ : %d != V_NUM %d, inconsistent!" %(len(V_in),len(RndJ),V_NUM))
    E_in = np.array(list(E_in),np.int)
    RndE = np.array(list(RndE),np.int)
    C_in = np.array(list(C_in),np.int)
    RndC = np.array(list(RndC),np.int)
    #From Here : Same as GetMatrix_ForceEstimation
    INV_NUM = V_NUM - len(RndJ)
    C_NUM = 2*(INV_NUM);       # 条件数
    X_NUM = E_NUM + CELL_NUMBER; # 未知変数

    MX = np.zeros( (V_NUM,X_NUM) );  #  V_NUM = INV_NUM+R_NUM
    MY = np.zeros( (V_NUM,X_NUM) );

    for (i, te) in enumerate(edge):
        MX[te.junc1,i] = -te.dx/te.dist;
        MY[te.junc1,i] = -te.dy/te.dist;
        MX[te.junc2,i] =  te.dx/te.dist;
        MY[te.junc2,i] =  te.dy/te.dist;

    for (i,c) in enumerate(cell):
        pn = c.jnum-1
        ij = c.junc

        hx = 0.0
        hy = 0.0
        if len(ij) == pn:
            print('\'not consistent polygon class. Exit.\'');
            sys.exit()

        if  not (ij[0] in RndJ):
            MX[ij[0],E_NUM+i] = 0.5*( y[ij[1]]-y[ij[pn]] )
            MY[ij[0],E_NUM+i] = 0.5*( x[ij[pn]]-x[ij[1]] )
            hx = hx + MX[ij[0],E_NUM+i]
            hy = hy + MY[ij[0],E_NUM+i]

        for j in range(1,pn):
            MX[ij[j],E_NUM+i] = 0.5*( y[ij[j+1]]-y[ij[j-1]] )
            MY[ij[j],E_NUM+i] = 0.5*( x[ij[j-1]]-x[ij[j+1]] )
            hx = hx + MX[ij[j],E_NUM+i]
            hy = hy + MY[ij[j],E_NUM+i]

        if not (ij[pn] in RndJ):
            MX[ij[pn],E_NUM+i] = 0.5*( y[ij[0]]-y[ij[pn-1]] )
            MY[ij[pn],E_NUM+i] = 0.5*( x[ij[pn-1]]-x[ij[0]] )
            hx = hx + MX[ij[pn],E_NUM+i]
            hy = hy + MY[ij[pn],E_NUM+i]

    #To Here : Same as GetMatrix_ForceEstimation
    
    #MX = np.delete( MX, RndJ, 0 )
    #MY = np.delete( MY, RndJ, 0 )
    MX = MX[V_in,:]
    MY = MY[V_in,:]
    MMC = np.concatenate( (MX, MY), 0) #Matrix for checking
    ME = MMC[:,:E_NUM]
    MC = MMC[:,E_NUM:]
    #ME = np.delete( ME, RndE, 1)
    #MC = np.delete( MC, RndC, 1)
    ME = ME[:,E_in]
    MC = MC[:,C_in]
    MM = np.concatenate( (ME,MC), 1)
    
    #Re-define Condition & Unknowns
    C_NUM = 2*( V_NUM - len(RndJ)) 
    X_NUM = (E_NUM - len(RndE)) + (CELL_NUMBER - len(RndC))
    
    if MM.shape != (C_NUM, X_NUM):
        print('MM.shape', MM.shape)
        print('\'Not Valid Martix: incorrect size\'');
        sys.exit()

    
    ###    Check the validity of Matrix MM
    #RndE2 : edges only one of whose vertices are V_in
    RndE2 = np.empty(0,np.int)
    for i in E_in:
        if edge[i].junc1 in RndJ:
            RndE2 = np.append(RndE2,i)
        elif edge[i].junc2 in RndJ:
            RndE2 = np.append(RndE2,i)

    RndE2 = np.concatenate((RndE,RndE2))
    
    ze = np.arange( E_NUM )
    ze = np.delete( ze, RndE2, 0 )
    zc = np.arange( CELL_NUMBER )
    zc = np.delete( zc, RndC, 0 )


    CM = np.sum( MX, 0 )
    cc = np.hstack( (CM[ze],CM[zc+E_NUM]) )
    err_x = np.abs( np.sum(cc) )

    
    CM = np.sum(MY,0)
    cc =  np.hstack( (CM[ze],CM[zc+E_NUM]) )
    err_y = np.abs( np.sum(cc) )

    # not boundary vertces
    ix = np.delete( x, RndJ, 0)
    iy = np.delete( y, RndJ, 0)
    cx = np.hstack( (iy,-ix) )
    ccx = cx.dot(MMC)
    ccx = np.delete( ccx, E_NUM+RndC, 0 )
    ccx = np.delete( ccx, RndE2, 0 )
    err_ang = np.abs( np.sum(ccx) )

    #     print( '%e %e %e' % (err_y, err_x,err_ang) )
    # ;;;

    ce = np.hstack( ( np.zeros((E_NUM - len(RndE))),np.ones((CELL_NUMBER - len(RndC)))) )
    err_iso = sum(abs(MM.dot(ce)))

    print('## Either err_ang or err_iso should be zero ( >ERR_MAX= %.2e ) ' % (ERR_MAX) );
    print('## err_x= %e   err_y= %e   err_ang=  %e  err_iso= %e\n' % (err_x,err_y,err_ang,err_iso) );
    if max([err_x, err_y, err_ang, err_iso]) > ERR_MAX:
        print('larger than %e : errors %e %e %e %e\n' % (ERR_MAX,err_x,err_y,err_ang,err_iso) )
        sys.exit()
    
    if SPARSE:
        sMM = sp.csr_matrix(MM)
        return [sMM, C_NUM, X_NUM]
    else:
        return [MM,C_NUM,X_NUM,[V_in,E_in,C_in],[RndJ,RndE,RndC],INV_NUM]





