# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 18:12:53 2022

@author: Daniel
"""

import numpy as np


global Exact_FI_p,Exact_FI_u,Exact_FI_v
global nx,ny
global x,y,U,U1,F,G,Flux,Exact_Flux,Error_Flux
global Beta,Re,P0,U0,V0,L,h,A1,dx,dy
nx=20
ny=20

P0=1
U0=1
V0=1
L=1
h=1
Beta=1
Re=10
A1=0
dx=L/nx
dy=h/ny
pi=np.pi

x=np.zeros(nx+2)
y=np.zeros(ny+2)
U=np.zeros((nx+2,ny+2,3))
U1=np.zeros((nx+2,ny+2,3))
F=np.zeros((nx+2,ny+2,3))
G=np.zeros((nx+2,ny+2,3))
Flux=np.zeros((nx+2,ny+2,3))
Exact_Flux=np.zeros((nx+2,ny+2,3))
Error_Flux=np.zeros((nx+2,ny+2,3))

pi=np.pi
def mesh(xm,ym,Rei,Beti):
    global nx,ny
    global x,y,U,U1,F,G,Flux,Exact_Flux,Error_Flux
    global Beta,Re,P0,U0,V0,L,h,A1,dx,dy
    nx=xm
    ny=ym
 
    P0=1
    U0=1
    V0=1
    L=1
    h=1
    Beta=Beti
    Re=Rei
    A1=0
    dx=L/nx
    dy=h/ny
    
    
    x=np.zeros(nx+2)
    y=np.zeros(ny+2)
    U=np.zeros((nx+2,ny+2,3))
    U1=np.zeros((nx+2,ny+2,3))
    F=np.zeros((nx+2,ny+2,3))
    G=np.zeros((nx+2,ny+2,3))
    Flux=np.zeros((nx+2,ny+2,3))
    Exact_Flux=np.zeros((nx+2,ny+2,3))
    Error_Flux=np.zeros((nx+2,ny+2,3))


def Exact_Flux_Integral(i,j):
    global Exact_FI_p,Exact_FI_u,Exact_FI_v
    x=(i-0.5)*dx
    y=(j-0.5)*dy
    Cx=np.cos(pi*x)
    Sx=np.sin(pi*x)
    Cy=np.cos(pi*y)
    Sy=np.sin(pi*y)
    C2x=np.cos(2*pi*x)
    S2x=np.sin(2*pi*x)
    C2y=np.cos(2*pi*y)
    S2y=np.sin(2*pi*y)
    Exact_FI_p=-pi*(U0*Cx*S2y+V0*S2x*Cy)/Beta
    Exact_FI_u=P0*pi*Sx*Cy-U0*U0*pi*S2y*S2y*S2x-U0*V0*pi*Sx*S2x*(Cy*S2y+2*C2y*Sy)-5*U0*pi*pi*Sx*S2y/Re
    Exact_FI_v = P0 * pi * Cx * Sy - V0 * V0 * pi * S2x * S2x * S2y- U0 * V0 * pi * Sy * S2y * (Cx * S2x + 2 * C2x * Sx)- 5 * V0 * pi * pi * Sy * S2x / Re

def Flux_Integral():

    global E_p,E_v,E_u,U,G,U1,Flux,Exact_Flux,Error_Flux
    E_p = E_v = E_u = 0;
    for i in range (1,nx+1,1):
        for j in range(1,ny+1,1):
            Exact_Flux_Integral(i, j)
            F[i][j][0] =1 / (2 * dx * Beta) * (U[i + 1][j][1] - U[i -1][j][1])
            
            
            G[i][j][0]= 1 / (2 * dy * Beta) * (U[i][j + 1][2] - U[i][j -1][2])
        
            Flux[i][j][0] = (-F[i][j][0] - G[i][j][0])
            
            E_p = E_p+ pow((Flux[i][j][0] - Exact_FI_p),2)
            
            F[i][j][1] = (1 / (4 * dx))* (pow(U[i + 1][j][1], 2) -\
            pow(U[i- 1][j][1], 2)+ 2 * U[i][j][1]* (U[i +1][j][1] - U[i - 1][j][1]))\
            + (U[i + 1][j][0] - U[i - 1][j][0])/ (2 * dx)- (1 / (Re * dx * dx))\
            * (U[i + 1][j][1] - 2 *U[i][j][1]+ U[i -1][j][1])
            
            
            G[i][j][1] =(1 / (4 * dy))* (U[i][j + 1][1] * U[i][j + 1][2]\
            + U[i][j + 1][1] *U[i][j][2]+ U[i][j][1] * U[i][j + 1][2]\
            - U[i][j][1] * U[i][j -1][2]- U[i][j - 1][1] *U[i][j][2]\
            - U[i][j - 1][1] *U[i][j - 1][2])- (1 / (Re * dy * dy))\
            * (U[i][j + 1][1] - 2 *U[i][j][1]+ U[i][j -1][1])
            
            Flux[i][j][1] = -F[i][j][1] - G[i][j][1]
            E_u = E_u+ (Flux[i][j][1] - Exact_FI_u)* (Flux[i][j][1] - Exact_FI_u)
            
            F[i][j][2] =1 / (4 * dx)* (U[i][j][1] * U[i + 1][j][2]+ U[i + 1][j][1] *\
            U[i][j][2] + U[i + 1][j][1] * U[i+ 1][j][2] - U[i][j][1] * U[i - 1][j][2]\
            - U[i - 1][j][1] * U[i][j][2] - U[i - 1][j][1] * U[i- 1][j][2])\
            - 1 / (Re * dx * dx) * (U[i + 1][j][2] - 2 * U[i][j][2]+ U[i -1][j][2])
                              
            G[i][j][2] =1 / (4 * dy)* (pow(U[i][j + 1][2], 2) -pow(U[i][j - 1][2], 2)\
            + 2 * U[i][j][2]* (U[i][j+ 1][2] - U[i][j - 1][2]))+ 1 / (2 * dy) * (U[i][j + 1][0] -\
            U[i][j - 1][0])- 1 / (Re * dy * dy)* (U[i][j + 1][2] - 2 * U[i][j][2]\
            + U[i][j - 1][2])                                                                      
            
           
            
        
            Flux[i][j][2] = -F[i][j][2] - G[i][j][2]
            E_v = E_v+ (Flux[i][j][2] - Exact_FI_v)* (Flux[i][j][2] - Exact_FI_v)
            
            
            Exact_Flux[i][j][0] = Exact_FI_p
            Exact_Flux[i][j][1] = Exact_FI_u
            Exact_Flux[i][j][2] = Exact_FI_v
            Error_Flux[i][j][0] = abs(abs(Flux[i][j][0]) - abs(Exact_Flux[i][j][0]))
            
            Error_Flux[i][j][1] = abs(abs(Flux[i][j][1]) - abs(Exact_Flux[i][j][1]))
            
            Error_Flux[i][j][2] = abs(abs(Flux[i][j][2]) - abs(Exact_Flux[i][j][2]))
            
           


def Flux_Error():
    global L2_p,L2_u,L2_v
    L2_p =np.sqrt(E_p / (nx * ny));
    L2_u = np.sqrt(E_u / (nx * ny));
    L2_v = np.sqrt(E_v / (nx * ny));


def Initialization():
    global x,y,U

    for i in range(0,nx+2,1):
        x[i]=(i - 0.5) * dx;
        for j in range(0,ny+2,1):
            y[j] = (j - 0.5) * dy;
    for i in range(0,nx+2,1):
        for j in range(0,ny+2,1):
            U[i][j][0] = P0 * np.cos(pi * x[i]) * np.cos(pi * y[j]);
            U[i][j][1] = U0 * np.sin(pi * x[i]) * np.sin(2 * pi * y[j]);
            U[i][j][2] = V0 * np.sin(2 * pi * x[i]) * np.sin(pi * y[j]);

import csv
import pandas as pd
print('Correctness of Residual:')   
print('')     
Rer=[10,10**-5,10**5]
Meshr=[20,40]  
Head=['Mesh size','$L_{2p}$','$L_{2u}$','$L_{2v}$'] 
lm=[]
lp=[]
lu=[]
lv=[] 
lr=[]
lb=[]
Filename='ValidL2.csv'
for j in Rer:
    
    for i in Meshr:
        lb.append(Beta)
        lr.append(j)
        
        lm.append(str(i)+r' '+r'$ \times $'+r' '+str(i))
        mesh(i,i,j,1) 
        Initialization()
        Flux_Integral()
        Flux_Error()
        lp.append(L2_p)
        lu.append(L2_u)
        lv.append(L2_v)
        R=str(j)
        M=str(i)
        print('Reynolds='+R+'  ,Mesh='+M+'x'+M)
        print('L2P=',L2_p)
        print('L2v=',L2_v)
        print('L2u=',L2_u)
        print("-----------------------------")
    Data=pd.DataFrame({'Mesh size':lm,'Reynolds':lr,r'$ B $':lb,r'$ L_{2p} $':lp,r'$ L_{2u} $':lu,'$ L_{2v} $':lv})
    Data.to_csv(Filename,index=False)

lp1=[]
lu1=[]
lv1=[]


for i in range (0,len(lp),2):
    a=lp[i]/lp[i+1]
    lp1.append(a)
    b=lu[i]/lu[i+1]
    lu1.append(b)
    c=lv[i]/lv[i+1]
    lv1.append(c)
''' 
Filename='4.csv'
Data=pd.DataFrame({'Reynolds':Rer,r'$ \dfrac{L_{2p}|_{20\times 20}}{L_{2p}|_{40\times 40}} $':lp1,r'$ \dfrac{L_{2u}|_{20\times 20}}{L_{2u}|_{40\times 40}} $':lu1,r'$ \dfrac{L_{2v}|_{20\times 20}}{L_{2v}|_{40\times 40}} $':lv1})
Data.to_csv(Filename,index=False)    
'''    
    