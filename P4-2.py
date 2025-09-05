# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 19:48:30 2022

@author: Daniel
"""
import numpy as np
import pandas as pd
s=20
m=s
n=s

L=1
H=1
dx=L/m
dy=H/n
Pi=np.pi


X_perturb=10
Y_perturb=10

P0=1
u0=1
v0=1
Beta=1
Re=10

x=np.zeros(m+2)
y=np.zeros(n+2)

F=np.zeros((m+2,n+2,3))
G=np.zeros((m+2,n+2,3))
Flux_n=np.zeros((m+2,n+2,3))
Flux_n1=np.zeros((m+2,n+2,3))
Flux=np.zeros((m+2,n+2,3))
U=np.zeros((m+3,n+3,3))
U1=np.zeros((m+3,n+3,3))
P=np.zeros((m+2,n+2))
u=np.zeros((m+2,n+2))
v=np.zeros((m+2,n+2))

LHS_5=np.zeros(((n+2)*(m+2),(n+2)*(m+2),3,3))
RHS_T=np.zeros(((n+2)*(m+2),3))
Result=np.zeros(3)
RHS_1=np.zeros(((n+2)*(m+2),3))
dU=np.zeros(((n+2)*(m+2),3))
Final_RHS=np.zeros((m+2,n+2,3))
Error=np.zeros((m + 2,n + 2,4))

Ax=np.zeros((m + 2,n + 2,3,3))
Ay=np.zeros((m + 2,n + 2,3,3))
Bx=np.zeros((m + 2,n +2,3,3))

Cx=np.zeros((m + 2,n + 2,3,3))
Cy=np.zeros((m + 2,n + 2,3,3))
By=np.zeros((m + 2,n + 2,3,3))


LHS_1=np.zeros((3,3))
dU_1=np.zeros(3)
Result_1=np.zeros(3)



def P_u_v():

    for i in range(0,m+2,1):  
        for j in range(0,n+2,1):        
            P[i][j] = P0*np.cos(Pi*(i - 0.5)*dx)*np.cos(Pi*(j - 0.5)*dy)
            u[i][j] = u0*np.sin(Pi*(i - 0.5)*dx)*np.sin(2 * Pi*(j -0.5)*dy)
            v[i][j] = v0*np.sin(2 * Pi*(i - 0.5)*dx)*np.sin(Pi*(j -0.5)*dy)
            
            
        
def Cofi():

    for i in range(1,m+1):
    
        for j in range(1,n+1):
            
        
            Ax[i][j][0][0] = 0.0
            Ax[i][j][0][1] = -1.0 / (2.0*Beta*dx)
            Ax[i][j][0][2] = 0.0
            Ax[i][j][1][0] = -0.5/dx
            Ax[i][j][1][1] = -((u[i][j] + u[i - 1][j]) / (2.0*dx)) - (1.0/(Re*dx*dx))
            
            Ax[i][j][1][2] = 0.0
            Ax[i][j][2][0] = 0.0
            Ax[i][j][2][1] = -(v[i][j] + v[i - 1][j]) / (4.0*dx)
            Ax[i][j][2][2] = -((u[i][j] + u[i - 1][j]) / (4.0*dx)) - (1.0/(Re*dx*dx));
            
            Cx[i][j][0][0] = 0.0
            Cx[i][j][0][1] = ((1.0 / (2.0*Beta*dx)))
            Cx[i][j][0][2] = 0.0
            Cx[i][j][1][0] = 0.5 / dx
            Cx[i][j][1][1] = (((u[i][j] + u[i + 1][j]) / 2.0) - (1.0 /(Re*dx))) / dx
            
            Cx[i][j][1][2] = 0.0
            Cx[i][j][2][0] = 0.0
            Cx[i][j][2][1] = (v[i][j] + v[i + 1][j]) / (4.0*dx)
            Cx[i][j][2][2] = (((u[i][j] + u[i + 1][j]) / 4.0) - (1.0 /(Re*dx))) / dx
            
            Ay[i][j][0][0] = 0
            Ay[i][j][0][1] = 0;
            Ay[i][j][0][2] = -(1.0 / (2.0*Beta*dy))
            Ay[i][j][1][0] = 0
            Ay[i][j][1][1] = -((v[i][j] + v[i][j - 1]) / (4.0*dy)) - (1.0/ (Re*dy*dy))
            
            Ay[i][j][1][2] = -(u[i][j] + u[i][j - 1]) / (4.0*dy)
            Ay[i][j][2][0] = -0.5 / dy
            Ay[i][j][2][1] = 0
            Ay[i][j][2][2] = -((v[i][j] + v[i][j - 1]) / (2.0*dy)) - (1.0/ (Re*dy*dy))
            
            Cy[i][j][0][0] = 0
            Cy[i][j][0][1] = 0
            Cy[i][j][0][2] = (1.0 / (2.0*Beta*dy))
            Cy[i][j][1][0] = 0
            Cy[i][j][1][1] = (((v[i][j] + v[i][j + 1]) / (4.0*dy)) - (1.0/ (Re*dy*dy))) 
            
            Cy[i][j][1][2] = (u[i][j] + u[i][j + 1]) / (4.0*dy)
            Cy[i][j][2][0] = 1 / (2*dy)
            Cy[i][j][2][1] = 0;
            Cy[i][j][2][2] = (((v[i][j] + v[i][j + 1]) /( 2.0*dy)) - (1.0/ (Re*dy*dy)))
            
            Bx[i][j][0][0] = 0
            Bx[i][j][0][1] = 0
            Bx[i][j][0][2] = 0
            Bx[i][j][1][0] = 0
            Bx[i][j][1][1] = (((u[i + 1][j] - u[i - 1][j]) / 2.0) + (2.0 /(Re*dx))) / dx
            
            Bx[i][j][1][2] = 0
            Bx[i][j][2][0] = 0
            Bx[i][j][2][1] = (v[i + 1][j] - v[i - 1][j]) / (4.0*dx) 
            Bx[i][j][2][2] = (((u[i + 1][j] - u[i - 1][j]) / 4.0) + (2.0 /(Re*dx))) / dx
            
            By[i][j][0][0] = 0
            By[i][j][0][1] = 0
            By[i][j][0][2] = 0
            By[i][j][1][0] = 0
            By[i][j][1][1] = (((v[i][j + 1] - v[i][j - 1]) / (4.0*dy)) +(2.0 / (Re*dy*dy)))
            
            By[i][j][1][2] = (u[i][j + 1] - u[i][j - 1]) / (4.0*dy)
            By[i][j][2][0] = 0
            By[i][j][2][1] = 0
            By[i][j][2][2] = (((v[i][j + 1] - v[i][j - 1]) / (2.0*dy)) +(2.0 / (Re*dy*dy)))
            
def Flux_Integral_n():
    
    for i in range(0,m+1):
        for j in range(0,n+1):
            
            F[i][j][0] =1 / (2 * dx * Beta) * (U[i + 1][j][1] - U[i -1][j][1])
            
            
            G[i][j][0] =1/(2*dy*Beta)*(U[i][j+1][2]-U[i][j-1][2])

            
            Flux_n[i][j][0] = (-F[i][j][0] - G[i][j][0])
            
            F[i][j][1] =(1 / (4 * dx))* (U[i + 1][j][1]*U[i + 1][j][1] -\
            U[i - 1][j][1]*U[i - 1][j][1]+ 2 * U[i][j][1]* (U[i +1][j][1] - U[i - 1][j][1]))\
            + (U[i + 1][j][0] - U[i - 1][j][0])/ (2 * dx)- (1 / (Re * dx * dx))\
            * (U[i + 1][j][1] - 2 *U[i][j][1]+ U[i -1][j][1])
            
            
            G[i][j][1] =(1 / (4 * dy))* (U[i][j + 1][1] * U[i][j + 1][2]+ U[i][j + 1][1] *\
            U[i][j][2]+ U[i][j][1] * U[i][j +1][2]- U[i][j][1] * U[i][j -1][2]\
            - U[i][j - 1][1] *U[i][j][2]- U[i][j - 1][1] *U[i][j - 1][2])\
            - (1 / (Re * dy * dy))* (U[i][j + 1][1] - 2 *U[i][j][1]+ U[i][j -1][1])
            
            
            Flux_n[i][j][1] = -F[i][j][1] - G[i][j][1]
            
            F[i][j][2] =1 / (4 * dx)* (U[i][j][1] * U[i + 1][j][2]+ U[i + 1][j][1] *\
            U[i][j][2]+ U[i + 1][j][1] * U[i+ 1][j][2]- U[i][j][1] * U[i -1][j][2]\
            - U[i - 1][j][1] *U[i][j][2]- U[i - 1][j][1] * U[i- 1][j][2])\
            - 1 / (Re * dx * dx)* (U[i + 1][j][2] - 2 *U[i][j][2]+ U[i -1][j][2])
            
            G[i][j][2] =1 / (4 * dy)* (U[i][j + 1][2]*U[i][j + 1][2] -U[i][j - 1][2]*U[i][j - 1][2]\
            + 2 * U[i][j][2]* (U[i][j+ 1][2] - U[i][j - 1][2]))+ 1 / (2 * dy) * (U[i][j + 1][0] -\
            U[i][j - 1][0])- 1 / (Re * dy * dy)* (U[i][j + 1][2] - 2 *U[i][j][2]+ U[i][j -1][2])
            
            
            Flux_n[i][j][2] = -F[i][j][2] - G[i][j][2]
            
            
def Flux_Integral_n1():

    for i in range(0,m+1):
        for j in range(0,n+1):
            F[i][j][0] =1 / (2 * dx * Beta) * (U1[i + 1][j][1] - U1[i -1][j][1])
            
            
            G[i][j][0] =1 / (2 * dy * Beta) * (U1[i][j + 1][2] - U1[i][j- 1][2])
            
            
            Flux_n1[i][j][0] = (-F[i][j][0] - G[i][j][0])
            
            F[i][j][1] =(1 / (4 * dx))* (U1[i + 1][j][1]*U1[i + 1][j][1] -\
            U1[i - 1][j][1] *U1[i - 1][j][1]+ 2 * U1[i][j][1]* (U1[i +1][j][1] - U1[i - 1][j][1]))\
            + (U1[i + 1][j][0] - U1[i -1][j][0]) / (2 * dx)- (1 / (Re * dx * dx))\
            * (U1[i + 1][j][1] - 2* U1[i][j][1]+ U1[i -1][j][1])
            
            G[i][j][1] =(1 / (4 * dy))* (U1[i][j + 1][1] * U1[i][j + 1][2]\
            + U1[i][j + 1][1] *U1[i][j][2]+ U1[i][j][1] * U1[i][j+ 1][2]\
            - U1[i][j][1] * U1[i][j- 1][2]- U1[i][j - 1][1] *U1[i][j][2]\
            - U1[i][j - 1][1] *U1[i][j - 1][2])- (1 / (Re * dy * dy))\
            * (U1[i][j + 1][1] - 2* U1[i][j][1]+ U1[i][j- 1][1])
            
            
            Flux_n1[i][j][1] = -F[i][j][1] - G[i][j][1]
            
            
            F[i][j][2] =1 / (4 * dx)* (U1[i][j][1] * U1[i + 1][j][2]+ U1[i + 1][j][1] *\
            U1[i][j][2]+ U1[i + 1][j][1] *U1[i + 1][j][2]- U1[i][j][1] * U1[i -1][j][2]\
            - U1[i - 1][j][1] *U1[i][j][2]- U1[i - 1][j][1] *U1[i - 1][j][2])- 1 / (Re * dx * dx)\
            * (U1[i + 1][j][2] - 2* U1[i][j][2]+ U1[i -1][j][2])    
            
            G[i][j][2] =1 / (4 * dy)* (U1[i][j + 1][2]*U1[i][j + 1][2] -U1[i][j - 1][2]*U1[i][j - 1][2]\
            + 2 * U1[i][j][2]* (U1[i][j+ 1][2] - U1[i][j - 1][2]))+ 1 / (2 * dy) * (U1[i][j + 1][0] -\
            U1[i][j - 1][0])- 1 / (Re * dy * dy)* (U1[i][j + 1][2] - 2* U1[i][j][2]\
            + U1[i][j- 1][2])
                                                   
            Flux_n1[i][j][2] = -F[i][j][2] - G[i][j][2]                                     
            
            
            
            
def Flux_sum():
    Flux_Integral_n()
    Flux_Integral_n1()
    for i in range(1,m+1):
        for j in range(1,n+1):
            Flux[i][j][0]=Flux_n1[i][j][0] - Flux_n[i][j][0]
            Flux[i][j][1]=Flux_n1[i][j][1] - Flux_n[i][j][1]
            Flux[i][j][2]=Flux_n1[i][j][2] - Flux_n[i][j][2]
   


def Initialization():

    for i in range(0,m+2):
        x[i] = (i - 0.5) * dx
    for j in range(0,n+2):
        y[j] = (j - 0.5) * dy
    for i in range(0,m+2):
        for j in range(0,n+2):
            U[i][j][0] = P0 * np.cos(Pi * x[i]) * np.cos(Pi * y[j])
            U[i][j][1] = u0 * np.sin(Pi * x[i]) * np.sin(2 * Pi * y[j])
            U[i][j][2] = v0 * np.sin(2 * Pi * x[i]) * np.sin(Pi * y[j])
            U1[i][j][0]=U[i][j][0]
            U1[i][j][1]=U[i][j][1]
            U1[i][j][2]=U[i][j][2]
            if i==X_perturb and j==Y_perturb:
                U1[i][j][0]=U1[i][j][0]+10**-6
                U1[i][j][1]=U1[i][j][1]+10**-6
                U1[i][j][2]=U1[i][j][2]+10**-6
                
                
def Multivec(A,Vec):
    Result[0] = A[0][0]*Vec[0] + A[0][1]*Vec[1] + A[0][2]*Vec[2]
    Result[1] = A[1][0]*Vec[0] + A[1][1]*Vec[1] + A[1][2]*Vec[2]
    Result[2] = A[2][0]*Vec[0] + A[2][1]*Vec[1] + A[2][2]*Vec[2]
                
         
def Left_HS():
    dU[Y_perturb*(m+2)+X_perturb+1][0]=1e-6;
    dU[Y_perturb*(m+2)+X_perturb+1][1]=1e-6; 
    dU[Y_perturb*(m+2)+X_perturb+1][2]=1e-6;
    for j in range(1,n+1):
        for i in range(1,m+1):
            
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][0][0]=Bx[i][j][0][0]+By[i][j][0][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][0][1]=Bx[i][j][0][1]+By[i][j][0][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][0][2]=Bx[i][j][0][2]+By[i][j][0][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][1][0]=Bx[i][j][1][0]+By[i][j][1][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][1][1]=Bx[i][j][1][1]+By[i][j][1][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][1][2]=Bx[i][j][1][2]+By[i][j][1][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][2][0]=Bx[i][j][2][0]+By[i][j][2][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][2][1]=Bx[i][j][2][1]+By[i][j][2][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+1][2][2]=Bx[i][j][2][2]+By[i][j][2][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][0][0]=Ax[i][j][0][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][0][1]=Ax[i][j][0][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][0][2]=Ax[i][j][0][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][1][0]=Ax[i][j][1][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][1][1]=Ax[i][j][1][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][1][2]=Ax[i][j][1][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][2][0]=Ax[i][j][2][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][2][1]=Ax[i][j][2][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i][2][2]=Ax[i][j][2][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][0][0]=Cx[i][j][0][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][0][1]=Cx[i][j][0][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][0][2]=Cx[i][j][0][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][1][0]=Cx[i][j][1][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][1][1]=Cx[i][j][1][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][1][2]=Cx[i][j][1][2]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][2][0]=Cx[i][j][2][0]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][2][1]=Cx[i][j][2][1]
            LHS_5[j*(m+2)+i+1][j*(m+2)+i+2][2][2]=Cx[i][j][2][2]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][0][0]=Cy[i][j][0][0]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][0][1]=Cy[i][j][0][1]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][0][2]=Cy[i][j][0][2]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][1][0]=Cy[i][j][1][0]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][1][1]=Cy[i][j][1][1]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][1][2]=Cy[i][j][1][2]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][2][0]=Cy[i][j][2][0]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][2][1]=Cy[i][j][2][1]
            LHS_5[j*(m+2)+i+1][(j+1)*(m+2)+i+1][2][2]=Cy[i][j][2][2]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][0][0]=Ay[i][j][0][0]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][0][1]=Ay[i][j][0][1]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][0][2]=Ay[i][j][0][2]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][1][0]=Ay[i][j][1][0]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][1][1]=Ay[i][j][1][1]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][1][2]=Ay[i][j][1][2]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][2][0]=Ay[i][j][2][0]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][2][1]=Ay[i][j][2][1]
            LHS_5[j*(m+2)+i+1][(j-1)*(m+2)+i+1][2][2]=Ay[i][j][2][2]
            
            
    for i in range(0,(n+2)*(m+2)):
        for j in range(0,(n+2)*(m+2)):
            for y in range(0,3):
                for t in range(0,3):
                    LHS_1[y][t]=LHS_5[i][j][y][t]
                
                dU_1[y]=dU[j][y]
                
            Multivec(LHS_1,dU_1);
            RHS_T[i][0]+=Result[0];
            RHS_T[i][1]+=Result[1];
            RHS_T[i][2]+=Result[2];
        
    
    for i1 in range(0,m+2):
        for j1 in range(0,n+2):
            Final_RHS[i1-1][j1][0]=RHS_T[(j1)*(m+2)+(i1)][0]
            Final_RHS[i1-1][j1][1]=RHS_T[(j1)*(m+2)+(i1)][1]
            Final_RHS[i1-1][j1][2]=RHS_T[(j1)*(m+2)+(i1)][2]
               
            
            
def error():
    global L_inf_p,L_inf_u,L_inf_v,L_2_p,L_2_u,L_2_v,L2p,L2u,L2v
    L_inf_p=0
    L_inf_u=0
    L_inf_v=0
    L_2_p=0
    L_2_u=0
    L_2_v=0;

    for i in range(1,m):
        
        for j in range(1,n): 
        
            Error[i][j][1] = abs(abs(Final_RHS[i][j][0] ) -abs(Flux[i][j][0]))
            
            Error[i][j][2] = abs(abs(Final_RHS[i][j][1] ) -abs(Flux[i][j][1]))
            
            Error[i][j][3] =abs( abs(Final_RHS[i][j][2]) - abs(Flux[i][j][2]))
            
            if Error[i][j][1] > L_inf_p :
                L_inf_p=Error[i][j][1]
                
            if Error[i][j][2] > L_inf_u:
                L_inf_u=Error[i][j][2]
                
            if Error[i][j][3] > L_inf_v :
            	L_inf_v=Error[i][j][3]
            
            L_2_p+=Error[i][j][1]*Error[i][j][1];
            L_2_u+=Error[i][j][2]*Error[i][j][2];
            L_2_v+=Error[i][j][3]*Error[i][j][3];
            
        
        
    
    L2p=pow(L_2_p/(n*m),0.5)
    L2u=pow(L_2_u/(n*m),0.5)
    L2v=pow(L_2_v/(n*m),0.5)
               
            
            
            
            
            
            
            

            
            
Initialization()
Flux_sum()
P_u_v()
Cofi()
Left_HS()
error();

print(Error[:,:,1])
print(Error[:,:,2])
print(Error[:,:,3])
print(Error[:,:,1])



'''
# df=pd.DataFrame(Flux[:,:,1])  
# df.to_csv('Flux.csv',index=False)

df=pd.DataFrame(Error[:,:,1])  
df.to_csv('EP.csv',index=False,header=None)

df=pd.DataFrame(Error[:,:,2])  
df.to_csv('EU.csv',index=False,header=None)

df=pd.DataFrame(Error[:,:,3])  
df.to_csv('EV.csv',index=False,header=None)
'''
           
            
            
            
            
           
            
              
