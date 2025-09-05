# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 19:19:07 2022

@author: Daniel
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


L = 1.0
h = 3

nx=10
ny=10*h

pi=np.pi

dt=0.05

MAXSIZE=ny+2

L2_p = L2_u = L2_v = 0.5

Res = 1e-6
P0 = 1.0
U0 = 1.0
V0 = 1.0


Ut = 1
w = 1.1
Beta = 1
Re = 100
A1 = 1
dx=L/nx
dy=h/ny



utum=1


x=np.zeros(nx + 2) 
y=np.zeros(ny + 2)
U=np.zeros((nx + 2,ny + 2,3))
U1=np.zeros((nx + 2,ny + 2,3))
F=np.zeros((nx + 2,ny + 2,3))
G=np.zeros((nx + 2,ny + 2,3))
Flux=np.zeros((nx + 2,ny + 2,3))
Exact_Flux=np.zeros((nx + 2,ny + 2,3))
Error_Flux=np.zeros((nx + 2,ny + 2,3))
LHS=np.zeros((ny + 2,3,3,3))
RHS=np.zeros((ny + 2,3))
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
        U[0][j][0] = U[1][j][0]
        U[0][j][1] = -1 * U[1][j][1]
        U[0][j][2] = -1 * U[1][j][2]
        U[nx + 1][j][0] = U[nx][j][0]
        U[nx + 1][j][1] = -1 * U[nx][j][1]
        U[nx + 1][j][2] = -1 * U[nx][j][2]
        
    #U[nx][2][0]=0 # Solving Pressure Drift Problem


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
            
            
    #U[nx][2][0]=0 #Solving Pressure Drift Problem



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
        U[0][j][1] = -1 * U[1][j][1];
        U[0][j][2] = -1 * U[1][j][2];
        U[nx + 1][j][0] = U[nx][j][0];
        U[nx + 1][j][1] = -1 * U[nx][j][1];
        U[nx + 1][j][2] = -1 * U[nx][j][2];
    
    U[0][5][0]=0;
    
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
        U1[0][j][1] = -1 * U1[1][j][1];
        U1[0][j][2] = -1 * U1[1][j][2];
        U1[nx + 1][j][0] = U1[nx][j][0];
        U1[nx + 1][j][1] = -1 * U1[nx][j][1];
        U1[nx + 1][j][2] = -1 * U1[nx][j][2];
    
    U1[0][5][0]=0;
    
                                                  

    
                                                  
                                                  
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
    
    Initialization();
    L2_u = 1;
    L2_p = 1;
    L2_v = 1;
    Iter = 0;
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
        
        
        
        Mean_Pressure();
        
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
            #p_r[i-1,j-1]=[x1,y1,U[i][j][0]];
            p_r[j-1,i-1]=U[i][j][0];
            u_r[j-1,i-1]=U[i][j][1];
            v_r[j-1,i-1]=U[i][j][2];
            

    
    
                
            
   


Nx = nx;
Ny = ny;
dx = L / Nx;
dy = h / Ny;


Initialization();

#Explicit_Solver();
Implicit_Solver();
Print_Result()



mpl.rc('text',usetex=True)

p1='pressure'+str(Ut)+'_'+str(nx)+'y'+str(h)+'.eps'
if Ut !=0:
    plt.contour(X1,Y1,p_r,colors='black',linestyles='solid')
ax1=plt.contourf(X1,Y1,p_r)
plt.colorbar(label='$P$')  
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis('equal')
#plt.savefig('./figs/'+p1,format='eps')
p1='pressure'+str(Ut)+'_'+str(nx)
plt.show()




u1='uvelocity'+str(Ut)+'_'+str(nx)+'y'+str(h)+'.eps'
if Ut !=0:
    plt.contour(X1,Y1,u_r,colors='black',linestyles='solid')
plt.contourf(X1,Y1,u_r)  
plt.colorbar(label='$u$')  
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis('equal') 
#plt.savefig('./figs/'+u1,format='eps')
plt.show()



v1='vvelocity'+str(Ut)+'_'+str(nx)+'y'+str(h)+'.eps'
if Ut !=0:
    plt.contour(X1,Y1,v_r,colors='black',linestyles='solid')
plt.contourf(X1,Y1,v_r)  
plt.colorbar(label='$v$')  
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis('equal') 
#plt.savefig('./figs/'+v1,format='eps')
plt.show()



uv='abs'+str(Ut)+'_'+str(nx)+'y'+str(h)+'.eps'
if Ut !=0:
    plt.streamplot(X1,Y1,u_r,v_r,color='black') 
plt.contourf(X1,Y1,np.sqrt(u_r**2+v_r**2))
plt.colorbar(label=r'$\sqrt{u^2+v^2}$')  
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis('equal')
#plt.savefig('./figs/'+uv,format='eps') 
plt.show()

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

# uv='abs'+str(Ut)+'_'+str(nx)+'y'+str(h)+'.eps'
# if Ut !=0:
#     plt.contour(X1,Y1,np.sqrt(u_r**2+v_r**2),100000,linewidths=0.8,colors='black',linestyles='solid')
# plt.contourf(X1,Y1,np.sqrt(u_r**2+v_r**2))  
# plt.colorbar(label=r'$\sqrt{u^2+v^2}$')  
# plt.xlabel('$x$')
# plt.ylabel('$y$')
# plt.axis('equal')
# #plt.savefig('./figs/'+uv,format='eps') 
# plt.show()

'''
filename='u_r_'+str(Ut)+'y'+str(h)+str(nx)
df1=pd.DataFrame(u_r)  
df1.to_csv(filename+'.csv',index=False,header=False)

filename='X'+'y'+str(h)+str(nx)
df2=pd.DataFrame(X1)  
df2.to_csv(filename+'.csv',index=False,header=False)

filename='Y'+'y'+str(h)+str(nx)
df2=pd.DataFrame(Y1)  
df2.to_csv(filename+'.csv',index=False,header=False)

AA1=pd.read_csv('u_r_120'+'.csv',sep=',',header=None)
AA1=AA1.to_numpy()

BB1=pd.read_csv('u_r_-120'+'.csv',sep=',',header=None)
BB1=BB1.to_numpy()

XX1=pd.read_csv('X20'+'.csv',sep=',',header=None)
XX1=XX1.to_numpy()


YY1=pd.read_csv('Y20'+'.csv',sep=',',header=None)
YY1=YY1.to_numpy()


plt.contour(XX1,YY1,AA1+BB1,colors='black',linestyles='solid')
plt.contourf(XX1,YY1,AA1+BB1)  
plt.colorbar(label=r'$u|_{U=1}+u|_{U=-1}$')  
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis('equal')
#plt.savefig('./figs/'+'uu1'+'.eps',format='eps') 
plt.show()
'''

'''
UU=(u_r[:,int(nx/2+1)]+u_r[:,int(nx/2-1)])/2    
# plt.plot(UU,Y1[:,int(nx/2+1)])
# plt.xlabel(r'$u$')
# plt.ylabel(r'$y$')
# plt.grid()
# plt.show()
UU=UU.reshape(-1,1)
Y2=Y1[:,int(nx/2+1)]
Y2=Y2.reshape(-1,1)
UU=np.concatenate((UU,Y2),axis=1)


filename2='UU'+str(Ut)+str(nx)+str(h)
df=pd.DataFrame(UU)  
df.to_csv(filename2+'.csv',index=False,header=False)



AA=pd.read_csv('UU1103'+'.csv',sep=',',header=None)
AA=AA.to_numpy()

BB=pd.read_csv('UU1203'+'.csv',sep=',',header=None)
BB=BB.to_numpy()


CC=pd.read_csv('UU1403'+'.csv',sep=',',header=None)
CC=CC.to_numpy()

# DD=pd.read_csv('UU1803'+'.csv',sep=',',header=None)
# DD=DD.to_numpy()


plt.plot(AA[:,0],AA[:,1])
plt.plot(BB[:,0],BB[:,1])
plt.plot(CC[:,0],CC[:,1])
#plt.plot(DD[:,0],DD[:,1])
plt.grid()
leg=[r'$10\times 30$',r'$20\times 60$',r'$40\times 120$']
plt.legend(leg)
plt.xlabel(r'$u$')
plt.ylabel(r'$y$')
#plt.savefig('./figs/'+'symh'+'.eps',format='eps') 
plt.show()
'''




'''
filename33='Iteration'+str(Ut)+str(nx)+'y'+str(h)
ITH=np.concatenate((Iterationi.reshape(-1,1),l2pii.reshape(-1,1),l2uii.reshape(-1,1),l2vii.reshape(-1,1)),axis=1)
df33=pd.DataFrame(ITH)
df33.to_csv(filename33+'.csv',index=False,header=False)



# III1=pd.read_csv('Iteration10'+'.csv',sep=',',header=None)
# III1=III1.to_numpy()
# plt.plot(III1[:,0],III1[:,1:4])
# leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
# plt.legend(leg)
# plt.grid()
# plt.xlabel(r'$Iterations$')
# plt.ylabel(r'$L_2$')
# plt.savefig('./figs/'+'It10',format='eps') 
# plt.show()



III2=pd.read_csv('Iteration120y3'+'.csv',sep=',',header=None)
III2=III2.to_numpy()
plt.plot(III2[:,0],III2[:,1:4])
leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
plt.legend(leg)
plt.grid()
plt.xlabel(r'$Iterations$')
plt.ylabel(r'$L_2$')
plt.yscale('log')
#plt.savefig('./figs/'+'It20y3'+'.eps',format='eps') 
plt.show()




# III0=pd.read_csv('Iteration020'+'.csv',sep=',',header=None)
# III0=III0.to_numpy()
# plt.plot(III0[:,0],III0[:,1:4])
# leg=[r'$L_{2}p$',r'$L_{2}u$',r'$L_{2}v$']
# plt.legend(leg)
# plt.grid()
# plt.xlabel(r'$Iterations$')
# plt.ylabel(r'$L_2$')
# plt.yscale('log')
# #plt.savefig('./figs/'+'It020'+'.eps',format='eps') 
# plt.show()


# CCC2=pd.read_csv('u_velocity'+'.csv',sep=',',header=None)
# CCC2=CCC2.to_numpy()


# CCC3=pd.read_csv('v_velocity'+'.csv',sep=',',header=None)
# CCC3=CCC3.to_numpy()


# llx=1
# lly=3
# nnx=320
# nny=3*nnx
# ddx=llx/nnx
# ddy=lly/nny
# XXX2=np.linspace(ddx/2,llx-ddx/2,nnx)
# YYY2=np.linspace(ddy/2,lly-ddy/2,nny)

# ZZZ2=np.zeros((nnx,nny))
# ZZZ3=np.zeros((nny,nnx))
# for i in range(0,nnx):
#     for j in range(0,nny):
#         ZZZ2[i,j]=CCC2[i*nny+j,2]
#         ZZZ3[j,i]=CCC3[i*nny+j,2]

# [X1,Y1]=np.meshgrid(XXX2,YYY2)
# ZZZ2=np.transpose(ZZZ2)
# plt.streamplot(X1,Y1,ZZZ2,ZZZ3,color='black',density=1) 
# #plt.contourf(X1,Y1,np.sqrt(u_r**2+v_r**2))
# #plt.streamplot(X1,Y1,ZZZ2,color='black')
# plt.contourf(X1,Y1,ZZZ2)  
# plt.colorbar(label=r'$\sqrt{u^2+v^2}$')  
# plt.xlabel('$x$')
# plt.ylabel('$y$')
# plt.axis('equal')
'''





