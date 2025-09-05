# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 19:19:07 2022

@author: Daniel
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from time import process_time


nx=20
ny=10

pi=np.pi

dt=0.05

MAXSIZE=ny+2

L2_p = L2_u = L2_v = 0.5

Res = 1e-6
P0 = 1.0
U0 = 1.0
V0 = 1.0
L = 8
h = 1.0
Ut = 0
w = 1.0
Beta = 0.92
Re = 50.0
A1 = 0.25
Uin=1
dx=L/nx
dy=h/ny



x=np.zeros(nx + 2) 
y=np.zeros(ny + 2)
U=np.zeros((nx + 2,ny + 2,3))
U1=np.zeros((nx + 2,ny + 2,3))
F=np.zeros((nx + 2,ny + 2,3))
G=np.zeros((nx + 2,ny + 2,3))
Flux=np.zeros((nx + 2,ny + 2,3))
Exact_Flux=np.zeros((nx + 2,ny + 2,3))
Error_Flux=np.zeros((nx + 2,ny + 2,3))
LHS=np.zeros((nx + 2,3,3,3))
RHS=np.zeros((nx + 2,3))
dU=np.zeros((nx + 2,ny + 2,3))
Ax=np.zeros((nx + 2,ny + 2,3,3))
Ay=np.zeros((nx + 2,ny + 2,3,3))
Bx=np.zeros((nx + 2,ny + 2,3,3))
By=np.zeros((nx + 2,ny + 2,3,3))
Cx=np.zeros((nx + 2,ny + 2,3,3))
Cy=np.zeros((nx + 2,ny + 2,3,3))



p_r=np.zeros((ny,nx))
u_r=np.zeros((ny,nx))
v_r=np.zeros((ny,nx))
X1=np.zeros((ny,nx))
Y1=np.zeros((ny,nx))


def Boundary_Condition():
    for i in range(1,nx+1): 
        U[i][0][0] = U[i][1][0]
        U[i][0][1] = -1 * U[i][1][1]
        U[i][0][2] = -1 * U[i][1][2]
        U[i][ny + 1][0] = U[i][ny][0]
        U[i][ny + 1][1] = 2 * Ut - 1 * U[i][ny][1]
        
        U[i][ny + 1][2] = -1 * U[i][ny][2]
        
    for j in range(1,ny+1): 
        U[0][j][0] = U[1][j][0];
        U[0][j][1] = 2 * Uin - 1 * U[1][j][1];
        U[0][j][2] = -1 * U[1][j][2];
        U[nx + 1][j][0] = -U[nx][j][0]; 
        U[nx + 1][j][1] = U[nx][j][1]; 
        U[nx + 1][j][2] = U[nx][j][2];
    U[1][0][0]=0;    


def Initialization():

    for i in range(0,nx+2): 
        x[i] = (i - 0.5) * dx
    for j in range(0,ny+2):
        y[j] = (j - 0.5) * dy
    for i in range(0,nx+2):
        for j in range(0,ny+2):
            U[i][j][0] = P0 * np.cos(pi * x[i]) * np.cos(pi * y[j])
            U[i][j][1] = U0 * np.sin(pi * x[i]) * np.sin(2 * pi * y[j])
            U[i][j][2] = V0 * np.sin(2 * pi * x[i]) * np.sin(pi * y[j])
            
            



def Flux_Integral():
    
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            F[i][j][0] =1 / (2 * dx * Beta) * (U[i + 1][j][1] - U[i -1][j][1])\
            - (A1) * (dy / dx)* (U[i - 1][j][0] - 2 *U[i][j][0]+ U[i +1][j][0])
            
            
            G[i][j][0] =1 / (2 * dy * Beta) * (U[i][j + 1][2] - U[i][j -1][2])\
            - A1 * (dx / dy)* (U[i][j - 1][0] - 2 *U[i][j][0]+ U[i][j +1][0])
            
            Flux[i][j][0] = -F[i][j][0] - G[i][j][0]
            
            
            F[i][j][1] =1 / (4 * dx)* (pow(U[i + 1][j][1], 2) - pow(U[i- 1][j][1], 2)\
            + 2 * U[i][j][1]* (U[i +1][j][1] - U[i - 1][j][1]))+ (U[i + 1][j][0] - U[i - 1][j][0])\
            / (2 * dx)- 1 / (Re * dx * dx)* (U[i + 1][j][1] - 2 *U[i][j][1]\
            + U[i -1][j][1])
                
                
            G[i][j][1] =1 / (4 * dy)* (U[i][j + 1][1] * U[i][j + 1][2] + U[i][j + 1][1] *\
            U[i][j][2]+ U[i][j][1] * U[i][j +1][2]- U[i][j][1] * U[i][j -1][2]\
            - U[i][j - 1][1] *U[i][j][2]- U[i][j - 1][1] *U[i][j - 1][2])\
            - 1 / (Re * dy * dy)* (U[i][j + 1][1] - 2 *U[i][j][1]+ U[i][j -1][1])
            
            
            Flux[i][j][1] = -F[i][j][1] - G[i][j][1]
            
            
            F[i][j][2] =1 / (4 * dx)* (U[i][j][1] * U[i + 1][j][2]+ U[i + 1][j][1] *\
            U[i][j][2]+ U[i + 1][j][1] * U[i+ 1][j][2]- U[i][j][1] * U[i -1][j][2]\
            - U[i - 1][j][1] *U[i][j][2]- U[i - 1][j][1] * U[i- 1][j][2])- 1 / (Re * dx * dx)\
            * (U[i + 1][j][2] - 2 *U[i][j][2]+ U[i -1][j][2])
            
            G[i][j][2] =1 / (4 * dy)* (pow(U[i][j + 1][2], 2) -pow(U[i][j - 1][2], 2)\
            + 2 * U[i][j][2]* (U[i][j+ 1][2] - U[i][j - 1][2]))+ 1 / (2 * dy) * (U[i][j + 1][0] -\
            U[i][j - 1][0])- 1 / (Re * dy * dy)* (U[i][j + 1][2] - 2 *U[i][j][2]\
            + U[i][j -1][2])                                      
            Flux[i][j][2] = -F[i][j][2] - G[i][j][2]
            
            


def Error_Computing():
    Error_p = Error_u = Error_v = 0;
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            Error_p = Error_p + pow(dU[i][j][0], 2);
            Error_u = Error_u + pow(dU[i][j][1], 2);
            Error_v = Error_v + pow(dU[i][j][2], 2);
            
    global L2_p,L2_u,L2_v        
    L2_p = np.sqrt(Error_p / (nx * ny));
    L2_u = np.sqrt(Error_u / (nx * ny));
    L2_v = np.sqrt(Error_v / (nx * ny));
    
def Implicit_Jacobian():
    for j in range(1,ny+1):
        for i in range(1,nx+1):
            Ax[i][j][0][0] = (dy / dx) * A1
            Ax[i][j][0][1] = -1 / (dx * (2 * Beta))
            Ax[i][j][1][0] = (-1 / dx) * 0.5
            Ax[i][j][1][1] = -1 / dx* (0.5 * (U[i - 1][j][1] + U[i][j][1]) + 1 / (Re* dx))
            
            
            Ax[i][j][2][1] = -(1 / dx )* 0.25 * (U[i - 1][j][2] +U[i][j][2])
            
            Ax[i][j][2][2] = -1 / dx* (0.25 * (U[i - 1][j][1] + U[i][j][1]) + 1 / (Re* dx))
            
            
            Ay[i][j][0][0] = (dx / dy) * A1
            Ay[i][j][0][2] = -1 / dy * 1 / (2 * Beta)
            Ay[i][j][1][1] = -1 / dy* (0.25 * (U[i][j][2] + U[i][j - 1][2]) + 1 / (Re* dy))
            
            
            Ay[i][j][1][2] = -1 / dy * 0.25 * (U[i][j][1] + U[i][j -1][1])
            
            Ay[i][j][2][0] = -1 / dy * 0.5
            Ay[i][j][2][2] = -1 / dy* (0.5 * (U[i][j - 1][2] + U[i][j][2]) + 1 / (Re* dy))
            
            
            Bx[i][j][0][0] = -2 * A1 * dy / dx
            Bx[i][j][0][1] = 0
            Bx[i][j][1][0] = 0
            Bx[i][j][1][1] = -1 / dx * (0.5 * U[i - 1][j][1] - 1 / (Re *dx))\
            + 1 / dx * (0.5 * U[i + 1][j][1] + 1 / (Re *dx))
            
            Bx[i][j][2][1] = -(1 / dx) * 0.25 * (U[i - 1][j][2] +U[i][j][2])\
            + (1 / dx) * 0.25 * (U[i][j][2] + U[i +1][j][2])
            
            Bx[i][j][2][2] = -1 / dx* (0.25 * (U[i - 1][j][1] + U[i][j][1]) - 1 / (Re* dx))\
            + 1 / dx* (0.25 * (U[i][j][1] + U[i +1][j][1])+ 1 / (Re * dx))
                       
            By[i][j][0][0] = -2 * A1 * dx / dy
            By[i][j][0][2] = 0
            By[i][j][1][1] = -1 / dy* (0.25 * (U[i][j][2] + U[i][j - 1][2]) - 1 / (Re* dy))\
            + 1 / dy* (0.25 * (U[i][j + 1][2] +U[i][j][2])+ 1 / (Re * dy))
            
            By[i][j][1][2] = -1 / dy * 0.25 * (U[i][j][1] + U[i][j -1][1])\
            + 1 / dy * 0.25 * (U[i][j + 1][1] + U[i][j][1])
            By[i][j][2][0] = 0
            By[i][j][2][2] = -1 / dy * (0.5 * U[i][j - 1][2] - 1 / (Re *dy))\
            + 1 / dy * (0.5 * U[i][j + 1][2] + 1 / (Re *dy))
            
            Cx[i][j][0][0] = (dy / dx) * A1
            Cx[i][j][0][1] = 1 / dx * 1 / (2 * Beta)
            Cx[i][j][1][0] = 1 / dx * 0.5
            Cx[i][j][1][1] = 1 / dx* (0.5 * (U[i + 1][j][1] + U[i][j][1]) - 1 / (Re* dx))
            
            Cx[i][j][2][1] = 1 / dx * 0.25 * (U[i][j][2] + U[i +1][j][2])
            
            Cx[i][j][2][2] = 1 / dx* (0.25 * (U[i][j][1] + U[i + 1][j][1]) - 1 / (Re* dx))
            Cy[i][j][0][0] = (dx / dy) * A1
            Cy[i][j][0][2] = 1 / dy * 1 / (2 * Beta)
            Cy[i][j][1][1] = 1 / dy* (0.25 * (U[i][j + 1][2] + U[i][j][2]) - 1 / (Re* dy))
            Cy[i][j][1][2] = 1 / dy * 0.25 * (U[i][j + 1][1] +U[i][j][1])            
            Cy[i][j][2][0] = 1 / dy * 0.5
            Cy[i][j][2][2] = 1 / dy* (0.5 * (U[i][j + 1][2] + U[i][j][2]) - 1 / (Re* dy))
            

def CopyVec(Source, Target):
    Target[0] = Source[0]
    Target[1] = Source[1]
    Target[2] = Source[2]

def Copy3x3(Source,Target):
    Target[0][0] = Source[0][0]
    Target[0][1] = Source[0][1]
    Target[0][2] = Source[0][2]
    Target[1][0] = Source[1][0]
    Target[1][1] = Source[1][1]
    Target[1][2] = Source[1][2]
    Target[2][0] = Source[2][0]
    Target[2][1] = Source[2][1]
    Target[2][2] = Source[2][2]

def Mult3x3(A,B, C):
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0]
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1]
    C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2]
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0]
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1]
    C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2]
    C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0]
    C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1]
    C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2]

def MultVec( A,  Vec, residult):
    residult[0] = A[0][0] * Vec[0] + A[0][1] * Vec[1] + A[0][2] * Vec[2]
    residult[1] = A[1][0] * Vec[0] + A[1][1] * Vec[1] + A[1][2] * Vec[2]
    residult[2] = A[2][0] * Vec[0] + A[2][1] * Vec[1] + A[2][2] * Vec[2]

def Add3x3( A,  B, Factor, C):

    C[0][0] = A[0][0] + Factor * B[0][0]
    C[0][1] = A[0][1] + Factor * B[0][1]
    C[0][2] = A[0][2] + Factor * B[0][2]
    C[1][0] = A[1][0] + Factor * B[1][0]
    C[1][1] = A[1][1] + Factor * B[1][1]
    C[1][2] = A[1][2] + Factor * B[1][2]
    C[2][0] = A[2][0] + Factor * B[2][0]
    C[2][1] = A[2][1] + Factor * B[2][1]
    C[2][2] = A[2][2] + Factor * B[2][2]
    
def AddVec( A,  B, Factor, C):

    C[0] = A[0] + Factor * B[0]
    C[1] = A[1] + Factor * B[1]
    C[2] = A[2] + Factor * B[2]

def Invert3x3( Block, Inverse):
    DetInv = 1.0/ (+Block[0][0] * Block[1][1] * Block[2][2]\
    + Block[0][1] * Block[1][2] * Block[2][0]\
    + Block[0][2] * Block[1][0] * Block[2][1]\
    - Block[0][2] * Block[1][1] * Block[2][0]\
    - Block[0][1] * Block[1][0] * Block[2][2]\
    - Block[0][0] * Block[1][2] * Block[2][1])
    Inverse[0][0] = +DetInv\
    * (Block[1][1] * Block[2][2] - Block[2][1] * Block[1][2]);
    Inverse[1][0] = -DetInv\
    * (Block[1][0] * Block[2][2] - Block[2][0] * Block[1][2]);
    Inverse[2][0] = +DetInv\
    * (Block[1][0] * Block[2][1] - Block[2][0] * Block[1][1]);
    Inverse[0][1] = -DetInv\
    * (Block[0][1] * Block[2][2] - Block[2][1] * Block[0][2]);
    Inverse[1][1] = +DetInv\
    * (Block[0][0] * Block[2][2] - Block[2][0] * Block[0][2]);
    Inverse[2][1] = -DetInv\
    * (Block[0][0] * Block[2][1] - Block[2][0] * Block[0][1]);
    Inverse[0][2] = +DetInv\
    * (Block[0][1] * Block[1][2] - Block[1][1] * Block[0][2]);
    Inverse[1][2] = -DetInv\
    * (Block[0][0] * Block[1][2] - Block[1][0] * Block[0][2]);
    Inverse[2][2] = +DetInv\
    * (Block[0][0] * Block[1][1] - Block[1][0] * Block[0][1]);
    

def SolveBlockTri(LHS, RHS,iNRows):
    
    
    Inv=np.zeros((3,3))
    for j in range(0,iNRows-1):
        Invert3x3(LHS[j][1], Inv)
        
        Temp=np.zeros((3,3))
        Mult3x3(Inv, LHS[j][2], Temp)
        Copy3x3(Temp, LHS[j][2])
        
        Temp=np.zeros(3)
        MultVec(Inv, RHS[j], Temp)
        CopyVec(Temp, RHS[j])
       
        A =LHS[j+1][0]
        B= LHS[j+1][1]
        C =LHS[j ][2]
        Temp=np.zeros((3,3))
        Temp2=np.zeros((3,3))
        TVec=np.zeros(3)
        TVec2=np.zeros(3)
        Mult3x3(A, C, Temp)
        Add3x3(B, Temp, -1., Temp2)
        Copy3x3(Temp2, B)
        MultVec(A, RHS[j], TVec)
        AddVec(RHS[j + 1], TVec, -1., TVec2)
        CopyVec(TVec2, RHS[j + 1])
        
    j = iNRows - 1;
    Invert3x3(LHS[j][1], Inv);
    
    Temp=np.zeros(3)
    MultVec(Inv, RHS[j], Temp)
    CopyVec(Temp, RHS[j])
    
    for j in range(iNRows - 2,-1,-1):
        C=LHS[j][2]
        RHS[j][0] -= (C[0][0]*RHS[j+1][0] +
        C[0][1]*RHS[j+1][1] +
        C[0][2]*RHS[j+1][2]);
        RHS[j][1] -= (C[1][0]*RHS[j+1][0] +
        C[1][1]*RHS[j+1][1] +
        C[1][2]*RHS[j+1][2]);
        RHS[j][2] -= (C[2][0]*RHS[j+1][0] +
        C[2][1]*RHS[j+1][1] +
        C[2][2]*RHS[j+1][2]);
       




            
def ADI_X():
    for m in range(1,ny+1):
        for i in range(1,nx+1):
            LHS[i][0][0][0] = dt * Ax[i][m][0][0]
            LHS[i][0][0][1] = dt * Ax[i][m][0][1]
            LHS[i][0][0][2] = dt * Ax[i][m][0][2]
            LHS[i][0][1][0] = dt * Ax[i][m][1][0]
            LHS[i][0][1][1] = dt * Ax[i][m][1][1]
            LHS[i][0][1][2] = dt * Ax[i][m][1][2]
            LHS[i][0][2][0] = dt * Ax[i][m][2][0]
            LHS[i][0][2][1] = dt * Ax[i][m][2][1]
            LHS[i][0][2][2] = dt * Ax[i][m][2][2]
            LHS[i][1][0][0] = 1 + dt * Bx[i][m][0][0]
            LHS[i][1][0][1] = dt * Bx[i][m][0][1]
            LHS[i][1][0][2] = dt * Bx[i][m][0][2]
            LHS[i][1][1][0] = dt * Bx[i][m][1][0]
            LHS[i][1][1][1] = 1 + dt * Bx[i][m][1][1]
            LHS[i][1][1][2] = dt * Bx[i][m][1][2]
            LHS[i][1][2][0] = dt * Bx[i][m][2][0]
            LHS[i][1][2][1] = dt * Bx[i][m][2][1]
            LHS[i][1][2][2] = 1 + dt * Bx[i][m][2][2]
            LHS[i][2][0][0] = dt * Cx[i][m][0][0]
            LHS[i][2][0][1] = dt * Cx[i][m][0][1]
            LHS[i][2][0][2] = dt * Cx[i][m][0][2]
            LHS[i][2][1][0] = dt * Cx[i][m][1][0]
            LHS[i][2][1][1] = dt * Cx[i][m][1][1]
            LHS[i][2][1][2] = dt * Cx[i][m][1][2]
            LHS[i][2][2][0] = dt * Cx[i][m][2][0]
            LHS[i][2][2][1] = dt * Cx[i][m][2][1]
            LHS[i][2][2][2] = dt * Cx[i][m][2][2]
            RHS[i][0] = dt * Flux[i][m][0]
            RHS[i][1] = dt * Flux[i][m][1]
            RHS[i][2] = dt * Flux[i][m][2]
        
        for l in range(0,3):
            for k in range(0,3): 
                LHS[0][0][l][k] = 0
                LHS[nx + 1][2][l][k] = 0
            
        
        LHS[0][1][0][0] = -1.0
        LHS[0][1][0][1] = 0.0
        LHS[0][1][0][2] = 0.0
        LHS[0][1][1][0] = 0.0
        LHS[0][1][1][1] = 1.0
        LHS[0][1][1][2] = 0.0
        LHS[0][1][2][0] = 0.0
        LHS[0][1][2][1] = 0.0
        LHS[0][1][2][2] = 1.0
        LHS[0][2][0][0] = 1.0
        LHS[0][2][0][1] = 0.0
        LHS[0][2][0][2] = 0.0
        LHS[0][2][1][0] = 0.0
        LHS[0][2][1][1] = 1.0
        LHS[0][2][1][2] = 0.0
        LHS[0][2][2][0] = 0.0
        LHS[0][2][2][1] = 0.0
        LHS[0][2][2][2] = 1.0
        RHS[0][0] = 0
        RHS[0][1] = 0
        RHS[0][2] = 0
        LHS[nx + 1][1][0][0] = 1.0
        LHS[nx + 1][1][0][1] = 0.0
        LHS[nx + 1][1][0][2] = 0.0
        LHS[nx + 1][1][1][0] = 0.0
        LHS[nx + 1][1][1][1] = 1.0
        LHS[nx + 1][1][1][2] = 0.0
        LHS[nx + 1][1][2][0] = 0.0
        LHS[nx + 1][1][2][1] = 0.0
        LHS[nx + 1][1][2][2] = 1.0
        LHS[nx + 1][0][0][0] = -1.
        LHS[nx + 1][0][0][1] = 0.
        LHS[nx + 1][0][0][2] = 0.
        LHS[nx + 1][0][1][0] = 0.
        LHS[nx + 1][0][1][1] = 1.
        LHS[nx + 1][0][1][2] = 0.
        LHS[nx + 1][0][2][0] = 0.
        LHS[nx + 1][0][2][1] = 0.
        LHS[nx + 1][0][2][2] = 1.
        RHS[nx + 1][0] = 0
        RHS[nx + 1][1] = 0
        RHS[nx + 1][2] = 0
        SolveBlockTri(LHS, RHS, nx + 2)
        for i in range(0,nx+2):
            dU[i][m][0] = RHS[i][0]
            dU[i][m][1] = RHS[i][1]
            dU[i][m][2] = RHS[i][2]
        

def ADI_Y():
    for p in range(1,nx+1):
        for j in range(1,ny+1):
            LHS[j][0][0][0] = dt * Ay[p][j][0][0];
            LHS[j][0][0][1] = dt * Ay[p][j][0][1];
            LHS[j][0][0][2] = dt * Ay[p][j][0][2];
            LHS[j][0][1][0] = dt * Ay[p][j][1][0];
            LHS[j][0][1][1] = dt * Ay[p][j][1][1];
            LHS[j][0][1][2] = dt * Ay[p][j][1][2];
            LHS[j][0][2][0] = dt * Ay[p][j][2][0];
            LHS[j][0][2][1] = dt * Ay[p][j][2][1];
            LHS[j][0][2][2] = dt * Ay[p][j][2][2];
            LHS[j][1][0][0] = 1 + dt * By[p][j][0][0];
            LHS[j][1][0][1] = dt * By[p][j][0][1];
            LHS[j][1][0][2] = dt * By[p][j][0][2];
            LHS[j][1][1][0] = dt * By[p][j][1][0];
            LHS[j][1][1][1] = 1 + dt * By[p][j][1][1];
            LHS[j][1][1][2] = dt * By[p][j][1][2];
            LHS[j][1][2][0] = dt * By[p][j][2][0];
            LHS[j][1][2][1] = dt * By[p][j][2][1];
            LHS[j][1][2][2] = 1 + dt * By[p][j][2][2];
            LHS[j][2][0][0] = dt * Cy[p][j][0][0];
            LHS[j][2][0][1] = dt * Cy[p][j][0][1];
            LHS[j][2][0][2] = dt * Cy[p][j][0][2];
            LHS[j][2][1][0] = dt * Cy[p][j][1][0];
            LHS[j][2][1][1] = dt * Cy[p][j][1][1];
            LHS[j][2][1][2] = dt * Cy[p][j][1][2];
            LHS[j][2][2][0] = dt * Cy[p][j][2][0];
            LHS[j][2][2][1] = dt * Cy[p][j][2][1];
            LHS[j][2][2][2] = dt * Cy[p][j][2][2];
            RHS[j][0] = dU[p][j][0];
            RHS[j][1] = dU[p][j][1];
            RHS[j][2] = dU[p][j][2];
        
        for l in range(0,3):
            for k in range(0,3):
                LHS[0][0][l][k] = 0;
                LHS[ny + 1][2][l][k] = 0;
        
        
        LHS[0][1][0][0] = -1.;
        LHS[0][1][0][1] = 0.;
        LHS[0][1][0][2] = 0.;
        LHS[0][1][1][0] = 0.;
        LHS[0][1][1][1] = 1.;
        LHS[0][1][1][2] = 0.;
        LHS[0][1][2][0] = 0.;
        LHS[0][1][2][1] = 0.;
        LHS[0][1][2][2] = 1.;
        LHS[0][2][0][0] = 1.;
        LHS[0][2][0][1] = 0.;
        LHS[0][2][0][2] = 0.;
        LHS[0][2][1][0] = 0.;
        LHS[0][2][1][1] = 1.;
        LHS[0][2][1][2] = 0.;
        LHS[0][2][2][0] = 0.;
        LHS[0][2][2][1] = 0.;
        LHS[0][2][2][2] = 1.;
        RHS[0][0] = 0.;
        RHS[0][1] = 0.;
        RHS[0][2] = 0.;
        LHS[ny + 1][1][0][0] = 1.;
        LHS[ny + 1][1][0][1] = 0.;
        LHS[ny + 1][1][0][2] = 0.;
        LHS[ny + 1][1][1][0] = 0.;
        LHS[ny + 1][1][1][1] = 1.;
        LHS[ny + 1][1][1][2] = 0.;
        LHS[ny + 1][1][2][0] = 0.;
        LHS[ny + 1][1][2][1] = 0.;
        LHS[ny + 1][1][2][2] = 1.;
        LHS[ny + 1][0][0][0] = -1.;
        LHS[ny + 1][0][0][1] = 0.;
        LHS[ny + 1][0][0][2] = 0.;
        LHS[ny + 1][0][1][0] = 0.;
        LHS[ny + 1][0][1][1] = 1.;
        LHS[ny + 1][0][1][2] = 0.;
        LHS[ny + 1][0][2][0] = 0.;
        LHS[ny + 1][0][2][1] = 0.;
        LHS[ny + 1][0][2][2] = 1.;
        RHS[ny + 1][0] = 0.;
        RHS[ny + 1][1] = 0.;
        RHS[ny + 1][2] = 0.;
        SolveBlockTri(LHS, RHS, ny + 2);
        for j in range(1,ny+1):
            dU[p][j][0] = RHS[j][0];
            dU[p][j][1] = RHS[j][1];
            dU[p][j][2] = RHS[j][2];

def Mean_Pressure():
    P_sum=0
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            P_sum = (U[i][j][0]) + P_sum
            
            
    P_avg = P_sum / (nx * ny);
    
def Implicit_Solver():
    global Iterationi,l2pii,l2uii,l2vii
    Iterationi=np.array([])
    l2pii=np.array([])
    l2uii=np.array([])
    l2vii=np.array([])
    Iter=0
    #t1=process_time()
    while (L2_p >= Res or L2_u >= Res or L2_v >= Res):
        Boundary_Condition();
        P_sum = 0;
        Flux_Integral();
        Implicit_Jacobian();
        ADI_X();
        ADI_Y();
        
        for k in range(0,3):
            for i in range(1,nx+1):
                for j in range(1,ny+1):
                        U[i][j][k] = U[i][j][k] + w * dU[i][j][k];
        
        Error_Computing();
        
        Mean_Pressure();
        Iter+=1;
        print(L2_p)
        Iterationi=np.append(Iterationi,Iter)
        l2pii=np.append(l2pii,L2_p)
        l2uii=np.append(l2uii,L2_u)
        l2vii=np.append(l2vii,L2_v)
    #t2=process_time()
    #print('CPU Time=',t2-t1)
    print(L2_p)
    print(L2_u)    
    print(L2_v)
    print(Iter)        
    
    
    


def Explicit_BC():
    for i in range(1,nx+1):  
        U[i][0][0] = U[i][1][0];
        U[i][0][1] = -1 * U[i][1][1];
        U[i][0][2] = -1 * U[i][1][2];
        U[i][ny + 1][0] = U[i][ny][0];
        U[i][ny + 1][1] = 2 * Ut - 1 * U[i][ny][1];
        U[i][ny + 1][2] = -1 * U[i][ny][2];
    
    for j in range(1,ny+1):
        U[0][j][0] = U[1][j][0];
        U[0][j][1] = 2 * Uin - 1 * U[1][j][1];
        U[0][j][2] = -1 * U[1][j][2];
        U[nx + 1][j][0] = -U[nx][j][0];
        U[nx + 1][j][1] = U[nx][j][1]; 
        U[nx + 1][j][2] = U[nx][j][2];
    
    U[1][0][0]=0;
    U1[1][0][0]=0;
    U[0][5][0]=0;


  

# def Explicit_BC():
#     for i in range(1,nx+1):  
#         U[i][0][0] = U[i][1][0];
#         U[i][0][1] = -1 * U[i][1][1];
#         U[i][0][2] = -1 * U[i][1][2];
#         U[i][ny + 1][0] = U[i][ny][0];
#         U[i][ny + 1][1] = 2 * Ut - 1 * U[i][ny][1];
#         U[i][ny + 1][2] = -1 * U[i][ny][2];
    
#     for j in range(1,ny+1):
#         U[0][j][0] = U[1][j][0];
#         U[0][j][1] = -1 * U[1][j][1];
#         U[0][j][2] = -1 * U[1][j][2];
#         U[nx + 1][j][0] = U[nx][j][0];
#         U[nx + 1][j][1] = -1 * U[nx][j][1];
#         U[nx + 1][j][2] = -1 * U[nx][j][2];
    
#     U[0][5][0]=0;
    
def Explicit_BC1():
    for i in range(1,nx+1): 
        U1[i][0][0] = U1[i][1][0];
        U1[i][0][1] = -1 * U1[i][1][1];
        U1[i][0][2] = -1 * U1[i][1][2];
        U1[i][ny + 1][0] = U1[i][ny][0];
        U1[i][ny + 1][1] = 2 * Ut - 1 * U1[i][ny][1];
        U1[i][ny + 1][2] = -1 * U1[i][ny][2];
    
    for j in range(1,ny+1):
        U1[0][j][0] = U1[1][j][0];
        U1[0][j][1] = 2 * Uin - 1 * U1[1][j][1];
        U1[0][j][2] = -1 * U1[1][j][2];
        U1[nx + 1][j][0] = -U1[nx][j][0];
        U1[nx + 1][j][1] = U1[nx][j][1];
        U1[nx + 1][j][2] = U1[nx][j][2];
    
    U[1][0][0]=0;
    U1[1][0][0]=0;
    U1[0][5][0]=0;
    
                                                  
# def Explicit_BC1():
#     for i in range(1,nx+1): 
#         U1[i][0][0] = U1[i][1][0];
#         U1[i][0][1] = -1 * U1[i][1][1];
#         U1[i][0][2] = -1 * U1[i][1][2];
#         U1[i][ny + 1][0] = U1[i][ny][0];
#         U1[i][ny + 1][1] = 2 * Ut - 1 * U1[i][ny][1];
#         U1[i][ny + 1][2] = -1 * U1[i][ny][2];
    
#     for j in range(1,ny+1):
#         U1[0][j][0] = U1[1][j][0];
#         U1[0][j][1] = -1 * U1[1][j][1];
#         U1[0][j][2] = -1 * U1[1][j][2];
#         U1[nx + 1][j][0] = U1[nx][j][0];
#         U1[nx + 1][j][1] = -1 * U1[nx][j][1];
#         U1[nx + 1][j][2] = -1 * U1[nx][j][2];
    
#     U1[0][5][0]=0;
    
    
    
    
    
    
                                                  
                                                  
def Explicit_RHS():
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            F[i][j][0] =1 / (2 * dx * Beta) * (U[i + 1][j][1] - U[i -1][j][1])\
            - (1.0 / 4) * (dy / dx)* (U[i - 1][j][0] - 2 *U[i][j][0]\
            + U[i +1][j][0]);
      
            G[i][j][0] =1 / (2 * dy * Beta) * (U[i][j + 1][2] - U[i][j -1][2])\
            - A1 * (dx / dy)* (U[i][j - 1][0] - 2 *U[i][j][0]+ U[i][j +1][0]);
            
            Flux[i][j][0] = -F[i][j][0] - G[i][j][0];
            
            F[i][j][1] = 1 / (4 * dx)* (pow(U[i + 1][j][1], 2) - pow(U[i- 1][j][1], 2)\
            + 2 * U[i][j][1]* (U[i +1][j][1] - U[i - 1][j][1]))+ (U[i + 1][j][0] - U[i - 1][j][0])\
            / (2 * dx)- 1 / (Re * dx * dx)* (U[i + 1][j][1] - 2 *U[i][j][1]+ U[i -1][j][1]);
            
 
            G[i][j][1] =1 / (4 * dy)* (U[i][j + 1][1] * U[i][j + 1][2]+ U[i][j + 1][1] *\
            U[i][j][2]+ U[i][j][1] * U[i][j +1][2]- U[i][j][1] * U[i][j -1][2]\
            - U[i][j - 1][1] *U[i][j][2]- U[i][j - 1][1] *U[i][j - 1][2])\
            - 1 / (Re * dy * dy)* (U[i][j + 1][1] - 2 *U[i][j][1]+ U[i][j -1][1]);
            
      
            Flux[i][j][1] = -F[i][j][1] - G[i][j][1];
            
            F[i][j][2] =1 / (4 * dx)* (U[i][j][1] * U[i + 1][j][2]+ U[i + 1][j][1] *\
            U[i][j][2]+ U[i + 1][j][1] * U[i+ 1][j][2]- U[i][j][1] * U[i -1][j][2]\
            - U[i - 1][j][1] *U[i][j][2]- U[i - 1][j][1] * U[i- 1][j][2])\
            - 1 / (Re * dx * dx)* (U[i + 1][j][2] - 2 *U[i][j][2]+ U[i -1][j][2]);
            
       
            G[i][j][2] =1 / (4 * dy)* (pow(U[i][j + 1][2], 2) -pow(U[i][j - 1][2], 2)\
            + 2 * U[i][j][2]* (U[i][j+ 1][2] - U[i][j - 1][2]))+ 1 / (2 * dy) * (U[i][j + 1][0] -\
            U[i][j - 1][0])- 1 / (Re * dy * dy)* (U[i][j + 1][2] - 2 *U[i][j][2]\
            + U[i][j -1][2]);
            
            Flux[i][j][2] = -F[i][j][2] - G[i][j][2];
            

                                                  
def Explicit_RHS1():
    for i in range(1,nx+1): 
        for j in range(1,ny+1):
            F[i][j][0] = 1 / (2 * dx * Beta)* (U1[i + 1][j][1] - U1[i - 1][j][1])\
            - (1.0 / 4) * (dy / dx)* (U1[i - 1][j][0] - 2 * U1[i][j][0]+ U1[i + 1][j][0]);
            
            G[i][j][0] = 1 / (2 * dy * Beta)* (U1[i][j + 1][2] - U1[i][j - 1][2])\
            - A1 * (dx / dy)* (U1[i][j - 1][0] - 2 * U1[i][j][0]+ U1[i][j + 1][0]);
            
            Flux[i][j][0] = -F[i][j][0] - G[i][j][0];
            
            F[i][j][1] = 1 / (4 * dx)* (pow(U1[i + 1][j][1], 2) - pow(U1[i - 1][j][1],2)\
            + 2 * U1[i][j][1]* (U1[i + 1][j][1] -U1[i - 1][j][1]))+ (U1[i + 1][j][0] - U1[i - 1][j][0]) / (2 * dx)\
            - 1 / (Re * dx * dx)* (U1[i + 1][j][1] - 2 * U1[i][j][1]+ U1[i - 1][j][1]);

            G[i][j][1] = 1 / (4 * dy)* (U1[i][j + 1][1] * U1[i][j + 1][2]\
            + U1[i][j + 1][1] * U1[i][j][2]+ U1[i][j][1] * U1[i][j + 1][2]\
            - U1[i][j][1] * U1[i][j - 1][2]- U1[i][j - 1][1] * U1[i][j][2]\
            - U1[i][j - 1][1] * U1[i][j - 1][2])- 1 / (Re * dy * dy)\
            * (U1[i][j + 1][1] - 2 * U1[i][j][1]+ U1[i][j - 1][1]);
            
   
            Flux[i][j][1] = -F[i][j][1] - G[i][j][1];
            
            
            F[i][j][2] = 1 / (4 * dx)* (U1[i][j][1] * U1[i + 1][j][2]\
            + U1[i + 1][j][1] * U1[i][j][2]+ U1[i + 1][j][1] * U1[i + 1][j][2]\
            - U1[i][j][1] * U1[i - 1][j][2]- U1[i - 1][j][1] * U1[i][j][2]\
            - U1[i - 1][j][1] * U1[i - 1][j][2])- 1 / (Re * dx * dx)\
            * (U1[i + 1][j][2] - 2 * U1[i][j][2]+ U1[i - 1][j][2]);
            

            G[i][j][2] = 1 / (4 * dy)* (pow(U1[i][j + 1][2], 2) - pow(U1[i][j - 1][2],2)\
            + 2 * U1[i][j][2]* (U1[i][j + 1][2] -U1[i][j - 1][2]))+ 1 / (2 * dy) * (U1[i][j + 1][0] - U1[i][j -1][0])\
            - 1 / (Re * dy * dy)* (U1[i][j + 1][2] - 2 * U1[i][j][2]+ U1[i][j - 1][2]);
            

            Flux[i][j][2] = -F[i][j][2] - G[i][j][2];
            
def Explicit_RK4(n):
    for i in range(1,nx+1):
        for j in range(1,ny+1): 
            U1[i][j][0] = U[i][j][0] + dt * Flux[i][j][0] / n;
            U1[i][j][1] = U[i][j][1] + dt * Flux[i][j][1] / n;
            U1[i][j][2] = U[i][j][2] + dt * Flux[i][j][2] / n;
            
    U1[0][5][0]=U[0][5][0]=0;    
    Explicit_BC1();
    Explicit_RHS1();
            

def Explicit_Solver():
    global Iterationi,l2pii,l2uii,l2vii
    Iterationi=np.array([])
    l2pii=np.array([])
    l2uii=np.array([])
    l2vii=np.array([])
    
    Initialization();
    L2_u = 1;
    L2_p = 1;
    L2_v = 1;
    Iter = 0;
    #t1=process_time()
    while (L2_p > Res or L2_u > Res or L2_v > Res):
        Iter+=1;
        Explicit_BC();
        Explicit_RHS();
        Explicit_RK4(4.0);
        Explicit_RK4(3.0);
        Explicit_RK4(2.0);
        
        L2_u = 0;
        L2_v = 0;
        L2_p = 0;
        for i in range(1,nx+1): 
            for j in range(1,ny+1):
                U1[i][j][0] = U[i][j][0] + dt * Flux[i][j][0];
                U1[i][j][1] = U[i][j][1] + dt * Flux[i][j][1];
                U1[i][j][2] = U[i][j][2] + dt * Flux[i][j][2];
                L2_p = L2_p + pow(dt * Flux[i][j][0], 2);
                L2_u = L2_u + pow(dt * Flux[i][j][1], 2);
                L2_v = L2_v + pow(dt * Flux[i][j][2], 2);
                
                U[i][j][0] = U1[i][j][0];
                U[i][j][1] = U1[i][j][1];
                U[i][j][2] = U1[i][j][2];
                
                U[i][j][0] = U[i][j][0] + w * (U[i][j][0] -U1[i][j][0]);
                
                U[i][j][1] = U[i][j][1] + w * (U[i][j][1] -U1[i][j][1]);
                
                U[i][j][2] = U[i][j][2] + w * (U[i][j][2] -U1[i][j][2]);
                
               
        
        L2_u = np.sqrt((L2_u) / (nx * ny));
        L2_v = np.sqrt((L2_v) / (nx * ny));
        L2_p = np.sqrt((L2_p) / (nx * ny));
        Iterationi=np.append(Iterationi,Iter)
        l2pii=np.append(l2pii,L2_p)
        l2uii=np.append(l2uii,L2_u)
        l2vii=np.append(l2vii,L2_v)
        print(L2_p)
        
        
        Mean_Pressure();
    
    #t2=process_time()
    #print('CPU Time=',t2-t1)
       
    print(L2_p)
    print(L2_u)    
    print(L2_v)
    print(Iter)
   
def Print_Result() :
    k=-1
    for i in range(1,nx+1): 
        for j in range(1,ny+1): 
            k=k+1
            x1 = dx * (i - 0.5);
            y1 = dy * (j - 0.5);
            X1[:,i-1]=x1
            Y1[j-1,:]=y1
            p_r[j-1,i-1]=U[i][j][0];
            u_r[j-1,i-1]=U[i][j][1];
            v_r[j-1,i-1]=U[i][j][2];

def velocity_duct():
    u=np.zeros(ny+2)
    
    dp = (U[nx - 1][int(ny / 2)][0] - U[nx - 2][int(ny / 2)][0]) / dx;
    e = 0;
    for j in range(1,ny+1):
        y1 = (j - 0.5) * dy;
        u[j] = (-Re / 2.0) * dp * y1 * (h - y1);
        e = e + (U[nx ][j][1] - u[j]) * (U[nx ][j][1] - u[j]);
    
    l2 = np.sqrt(e / ny);
    print("L2 compared with exact solution",l2)

Nx = nx;
Ny = ny;
dx = L / Nx;
dy = h / Ny;

t1=process_time()
Initialization();

#Explicit_Solver();
Implicit_Solver();
t2=process_time()
velocity_duct()
Print_Result()

print('CPU Time=',t2-t1)

mpl.rc('text',usetex=True)


ITH=np.concatenate((Iterationi.reshape(-1,1),l2pii.reshape(-1,1),l2uii.reshape(-1,1),l2vii.reshape(-1,1)),axis=1)
III1=ITH
plt.plot(III1[:,0],III1[:,1:4])
leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
plt.legend(leg)
plt.grid()
plt.xlabel(r'$Iterations$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
#plt.savefig('./figs/'+'It10',format='eps') 
plt.show()

'''
filename33='Iteration'+str(nx)+str(dt)+'D'+'exp1'
ITH=np.concatenate((Iterationi.reshape(-1,1),l2pii.reshape(-1,1),l2uii.reshape(-1,1),l2vii.reshape(-1,1)),axis=1)
df33=pd.DataFrame(ITH)
df33.to_csv(filename33+'.csv',index=False,header=False)
'''
'''

III2=pd.read_csv('Iteration200.05D'+'.csv',sep=',',header=None)
III2=III2.to_numpy()
plt.plot(III2[:,0],III2[:,1:4])
leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
plt.legend(leg)
plt.grid()
plt.xlabel(r'$Iterations$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
plt.savefig('./figs/'+'It20D0.05i'+'.eps',format='eps') 
plt.show()   
'''




'''
III3=pd.read_csv('Iteration800.02Dexp1'+'.csv',sep=',',header=None)
III3=III3.to_numpy()
plt.plot(III3[:,0],III3[:,1:4])
leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
plt.legend(leg)
plt.grid()
plt.xlabel(r'$Iterations$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
plt.savefig('./figs/'+'It80D0.02exp1'+'.eps',format='eps') 
plt.show()   

'''

'''
III2=pd.read_csv('Iteration200.05D'+'.csv',sep=',',header=None)
III2=III2.to_numpy()
plt.plot(III2[:,0],III2[:,1:4])
leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
plt.legend(leg)
plt.grid()
plt.xlabel(r'$Iterations$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
plt.savefig('./figs/'+'ItD0.05i'+'.eps',format='eps') 
plt.show()   


III3=pd.read_csv('Iteration200.05Dexp'+'.csv',sep=',',header=None)
III3=III3.to_numpy()
plt.plot(III3[:,0],III3[:,1:4])
leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
plt.legend(leg)
plt.grid()
plt.xlabel(r'$Iterations$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
plt.savefig('./figs/'+'ItD0.05exp'+'.eps',format='eps') 
plt.show()   
'''









Xt=7

ii=7/dx+0.5
kk=np.floor(ii)
vv=kk+1

if X1[0,int(ii)]==Xt:
    UR=u_r[:,int(ii)]
else:    
    UR=(u_r[:,int(kk)]+u_r[:,int(vv)])/2
    
YY=Y1[:,0]

p1='VelocityD'+str(Xt)+str(nx)
plt.plot(UR,YY)
plt.xlabel('$u$')
plt.ylabel('$y$')
plt.axis('equal')
plt.grid()
#plt.savefig('./figs/'+p1+'.eps',format='eps')
plt.show()

   
'''
UR=UR.reshape(-1,1)
Y2=YY.reshape(-1,1)
UU=np.concatenate((UR,Y2),axis=1)


filename2='UU'+str(nx)+'D'
df=pd.DataFrame(UU)  
df.to_csv(filename2+'.csv',index=False,header=False)



AA=pd.read_csv('UU20D'+'.csv',sep=',',header=None)
AA=AA.to_numpy()

BB=pd.read_csv('UU40D'+'.csv',sep=',',header=None)
BB=BB.to_numpy()


CC=pd.read_csv('UU80D'+'.csv',sep=',',header=None)
CC=CC.to_numpy()


yyy=np.linspace(0,1,1000)
uuu=6*yyy*(1-yyy)

plt.plot(AA[:,0],AA[:,1])
plt.plot(BB[:,0],BB[:,1])
plt.plot(CC[:,0],CC[:,1])
plt.plot(uuu,yyy)
plt.grid()
leg=[r'$20\times 10$',r'$40\times 20$',r'$80\times 40$',r'Fully-Developed Exact Solution']
plt.legend(leg)
plt.xlabel(r'$u$')
plt.ylabel(r'$y$')
#plt.savefig('./figs/'+'symD'+'.eps',format='eps') 
plt.show()
'''