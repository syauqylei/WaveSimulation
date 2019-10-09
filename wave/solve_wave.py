#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 16:48:46 2019

@author: rinn
"""
from numba import autojit, njit ,jit

@njit
def solve_wave_FD16(U,Up,h,c,ncpml,b,psi,phi):
    """ Calculate wave solution using 16th order finite difference
    """
    ny , nx = U.shape
    for i in range(8,ny-8):
        for j in range(8,nx-8):
            Up[i,j] = 2.0*U[i,j] - Up[i,j] + c[i-8,j-8]* \
                        ((-735*U[i-8,j]+15360*U[i-7,j]-156800*U[i-6,j]+1053696*U[i-5,j]-5350800*U[i-4,j]+22830080*U[i-3,j]-94174080*U[i-2,j]+538137600*U[i-1,j]-924708642*U[i+0,j]+538137600*U[i+1,j]-94174080*U[i+2,j]+22830080*U[i+3,j]-5350800*U[i+4,j]+1053696*U[i+5,j]-156800*U[i+6,j]+15360*U[i+7,j]-735*U[i+8,j])+ \
                         (-735*U[i,j-8]+15360*U[i,j-7]-156800*U[i,j-6]+1053696*U[i,j-5]-5350800*U[i,j-4]+22830080*U[i,j-3]-94174080*U[i,j-2]+538137600*U[i,j-1]-924708642*U[i,j+0]+538137600*U[i,j+1]-94174080*U[i,j+2]+22830080*U[i,j+3]-5350800*U[i,j+4]+1053696*U[i,j+5]-156800*U[i,j+6]+15360*U[i,j+7]-735*U[i,j+8]))/ \
                         (302702400*1.0*h**2)
    
    #CPML boundary in X-domain
    for i in range(8+ncpml,ny-ncpml-8):
        for j in range(8,ncpml+1):
            phi[i,j]=b[j-8]*phi[i,j]+(b[j-8]-1.0)*(U[i,j+1]-U[i,j])/h
            phi[i,-j-1]=b[j-8]*phi[i,-j-1]+(b[j-8]-1.0)*(U[i,-j-1]-U[i,-j-2])/h
        for j in range(8,ncpml):
            psi[i,j]=b[j-8]*psi[i,j]+(b[j-8]-1.0)*\
                                        ((U[i,j-1]-2*U[i,j]+U[i,j+1])/h/h \
                                         +(phi[i,j+1]-phi[i,j])/h)
            psi[i,-j-1]=b[j-8]*psi[i,-j-1]+(b[j-8]-1.0)*\
                                        ((U[i,-j-2]-2*U[i,-j-1]+U[i,-j])/h/h \
                                         +(phi[i,-j-1]-phi[i,-j-2])/h)
            Up[i,j] += c[i-8,j-8]*((phi[i,j+1]-phi[i,j])/h+psi[i,j])
            Up[i,-j-1] += c[i-8,-j-9]*((phi[i,-j-1]-phi[i,-j-2])/h+psi[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(8,ncpml+1):
        for j in range(8,nx-8):
            phi[i,j]=b[i-8]*phi[i,j]+(b[i-8]-1.0)*(U[i+1,j]-U[i,j])/h
            phi[-i-1,j]=b[i-8]*phi[-i-1,j]+(b[i-8]-1.0)*(U[-i-1,j]-U[-i-2,j])/h
    for i in range(8,ncpml):
        for j in range(8,nx-8):
            psi[i,j]=b[i-8]*psi[i,j]+(b[i-8]-1.0)*\
                                        ((U[i-1,j]-2*U[i,j]+U[i+1,j])/h/h \
                                         +(phi[i+1,j]-phi[i,j])/h)
            psi[-i-1,j]=b[i-8]*psi[-i-1,j]+(b[i-8]-1.0)*\
                                        ((U[-i-2,j]-2*U[-i-1,j]+U[-i,j])/h/h \
                                         +(phi[-i-1,j]-phi[-i-2,j])/h)
            
            Up[i,j] += c[i-8,j-8]*((phi[i+1,j]-phi[i,j])/h+psi[i,j])
            Up[-i-1,j] += c[-i-9,j-8]*((phi[-i-1,j]-phi[-i-2,j])/h+psi[-i-1,j])



@njit
def solve_wave_FD12(U,Up,h,c,ncpml,b,psi,phi):
    """ Calculate wave solution using 12th order finite difference
    """
    ny , nx = U.shape
    for i in range(6,ny-6):
        for j in range(6,nx-6):
            Up[i,j] = 2.0*U[i,j] - Up[i,j] + c[i-8,j-6]* \
                        ((-50*U[i-6,j]+864*U[i-5,j]-7425*U[i-4,j]+44000*U[i-3,j]-222750*U[i-2,j]+1425600*U[i-1,j]-2480478*U[i+0,j]+1425600*U[i+1,j]-222750*U[i+2,j]+44000*U[i+3,j]-7425*U[i+4,j]+864*U[i+5,j]-50*U[i+6,j])+ \
                         (-50*U[i,j-6]+864*U[i,j-5]-7425*U[i,j-4]+44000*U[i,j-3]-222750*U[i,j-2]+1425600*U[i,j-1]-2480478*U[i,j+0]+1425600*U[i,j+1]-222750*U[i,j+2]+44000*U[i,j+3]-7425*U[i,j+4]+864*U[i,j+5]-50*U[i,j+6]))/ \
                         (831600*1.0*h**2)
    
    #CPML boundary in X-domain
    for i in range(6+ncpml,ny-ncpml-6):
        for j in range(6,ncpml+1):
            phi[i,j]=b[j-1]*phi[i,j]+(b[j-1]-1.0)*(U[i,j+1]-U[i,j])/h
            phi[i,-j-1]=b[j-1]*phi[i,-j-1]+(b[j-1]-1.0)*(U[i,-j-1]-U[i,-j-2])/h
        for j in range(6,ncpml):
            psi[i,j]=b[j-1]*psi[i,j]+(b[j-1]-1.0)*\
                                        ((U[i,j-1]-2*U[i,j]+U[i,j+1])/h/h \
                                         +(phi[i,j+1]-phi[i,j])/h)
            psi[i,-j-1]=b[j-1]*psi[i,-j-1]+(b[j-1]-1.0)*\
                                        ((U[i,-j-2]-2*U[i,-j-1]+U[i,-j])/h/h \
                                         +(phi[i,-j-1]-phi[i,-j-2])/h)
            Up[i,j] += c[i-6,j-6]*((phi[i,j+1]-phi[i,j])/h+psi[i,j])
            Up[i,-j-1] += c[i-6,-j-7]*((phi[i,-j-1]-phi[i,-j-2])/h+psi[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(6,ncpml+1):
        for j in range(6,nx-6):
            phi[i,j]=b[i-1]*phi[i,j]+(b[i-1]-1.0)*(U[i+1,j]-U[i,j])/h
            phi[-i-1,j]=b[i-1]*phi[-i-1,j]+(b[i-1]-1.0)*(U[-i-1,j]-U[-i-2,j])/h
    for i in range(6,ncpml):
        for j in range(6,nx-6):
            psi[i,j]=b[i-1]*psi[i,j]+(b[i-1]-1.0)*\
                                        ((U[i-1,j]-2*U[i,j]+U[i+1,j])/h/h \
                                         +(phi[i+1,j]-phi[i,j])/h)
            psi[-i-1,j]=b[i-1]*psi[-i-1,j]+(b[i-1]-1.0)*\
                                        ((U[-i-2,j]-2*U[-i-1,j]+U[-i,j])/h/h \
                                         +(phi[-i-1,j]-phi[-i-2,j])/h)
            Up[i,j] += c[i-6,j-6]*((phi[i+1,j]-phi[i,j])/h+psi[i,j])
            Up[-i-1,j] += c[-i-7,j-6]*((phi[-i-1,j]-phi[-i-2,j])/h+psi[-i-1,j])


@njit
def solve_wave_FD8(U,Up,h,c,ncpml,b,psi,phi):
    """ Calculate wave solution using 8th finite difference
    """
    ny , nx = U.shape
    for i in range(4,ny-4):
        for j in range(4,nx-4):
            Up[i,j] = 2.0*U[i,j] - Up[i,j] + c[i-4,j-4]* \
                        ((-9*U[i,j-4]+128*U[i,j-3]-1008*U[i,j-2]+8064*U[i,j-1]-14350*U[i,j+0]+8064*U[i,j+1]-1008*U[i,j+2]+128*U[i,j+3]-9*U[i,j+4])+ \
                         (-9*U[i-4,j]+128*U[i-3,j]-1008*U[i-2,j]+8064*U[i-1,j]-14350*U[i+0,j]+8064*U[i+1,j]-1008*U[i+2,j]+128*U[i+3,j]-9*U[i+4,j]))/ \
                         (5040*1.0*h**2)
    
    #CPML boundary in X-domain
    for i in range(4+ncpml,ny-ncpml-4):
        for j in range(4,ncpml+1):
            phi[i,j]=b[j-1]*phi[i,j]+(b[j-1]-1.0)*(U[i,j+1]-U[i,j])/h
            phi[i,-j-1]=b[j-1]*phi[i,-j-1]+(b[j-1]-1.0)*(U[i,-j-1]-U[i,-j-2])/h
        for j in range(4,ncpml):
            psi[i,j]=b[j-1]*psi[i,j]+(b[j-1]-1.0)*\
                                        ((U[i,j-1]-2*U[i,j]+U[i,j+1])/h/h \
                                         +(phi[i,j+1]-phi[i,j])/h)
            psi[i,-j-1]=b[j-1]*psi[i,-j-1]+(b[j-1]-1.0)*\
                                        ((U[i,-j-2]-2*U[i,-j-1]+U[i,-j])/h/h \
                                         +(phi[i,-j-1]-phi[i,-j-2])/h)
            Up[i,j] += c[i-4,j-4]*((phi[i,j+1]-phi[i,j])/h+psi[i,j])
            Up[i,-j-1] += c[i-4,-j-5]*((phi[i,-j-1]-phi[i,-j-2])/h+psi[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(4,ncpml+1):
        for j in range(4,nx-4):
            phi[i,j]=b[i-1]*phi[i,j]+(b[i-1]-1.0)*(U[i+1,j]-U[i,j])/h
            phi[-i-1,j]=b[i-1]*phi[-i-1,j]+(b[i-1]-1.0)*(U[-i-1,j]-U[-i-2,j])/h
    for i in range(4,ncpml):
        for j in range(4,nx-4):
            psi[i,j]=b[i-1]*psi[i,j]+(b[i-1]-1.0)*\
                                        ((U[i-1,j]-2*U[i,j]+U[i+1,j])/h/h \
                                         +(phi[i+1,j]-phi[i,j])/h)
            psi[-i-1,j]=b[i-1]*psi[-i-1,j]+(b[i-1]-1.0)*\
                                        ((U[-i-2,j]-2*U[-i-1,j]+U[-i,j])/h/h \
                                         +(phi[-i-1,j]-phi[-i-2,j])/h)
            Up[i,j] += c[i-4,j-4]*((phi[i+1,j]-phi[i,j])/h+psi[i,j])
            Up[-i-1,j] += c[-i-5,j-4]*((phi[-i-1,j]-phi[-i-2,j])/h+psi[-i-1,j])
@njit
def solve_wave_FD4(U,Up,h,c,ncpml,b,psi,phi):
    """ Calculate wave solution using 4th order finite difference
    """
    ny , nx = U.shape
    for i in range(2,ny-2):
        for j in range(2,nx-2):
            Up[i,j] = 2.0*U[i,j] - Up[i,j] + c[i-2,j-2]* \
                        ((-1*U[i-2,j]+16*U[i-1,j]-30*U[i,j]+16*U[i+1,j]-1*U[i+2,j]) + \
                         (-1*U[i,j-2]+16*U[i,j-1]-30*U[i,j]+16*U[i,j+1]-1*U[i,j+2]))/ \
                         (12*1.0*h**2)
    #CPML boundary in X-domain
    for i in range(2+ncpml,ny-ncpml-2):
        for j in range(2,ncpml+1):
            phi[i,j]=b[j-1]*phi[i,j]+(b[j-1]-1.0)*(U[i,j+1]-U[i,j])/h
            phi[i,-j-1]=b[j-1]*phi[i,-j-1]+(b[j-1]-1.0)*(U[i,-j-1]-U[i,-j-2])/h
        for j in range(2,ncpml):
            psi[i,j]=b[j-1]*psi[i,j]+(b[j-1]-1.0)*\
                                        ((U[i,j-1]-2*U[i,j]+U[i,j+1])/h/h \
                                         +(phi[i,j+1]-phi[i,j])/h)
            psi[i,-j-1]=b[j-1]*psi[i,-j-1]+(b[j-1]-1.0)*\
                                        ((U[i,-j-2]-2*U[i,-j-1]+U[i,-j])/h/h \
                                         +(phi[i,-j-1]-phi[i,-j-2])/h)
            Up[i,j] += c[i-2,j-2]*((phi[i,j+1]-phi[i,j])/h+psi[i,j])
            Up[i,-j-1] += c[i-2,-j-3]*((phi[i,-j-1]-phi[i,-j-2])/h+psi[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(2,ncpml+1):
        for j in range(2,nx-2):
            phi[i,j]=b[i-1]*phi[i,j]+(b[i-1]-1.0)*(U[i+1,j]-U[i,j])/h
            phi[-i-1,j]=b[i-1]*phi[-i-1,j]+(b[i-1]-1.0)*(U[-i-1,j]-U[-i-2,j])/h
    for i in range(2,ncpml):
        for j in range(2,nx-2):
            psi[i,j]=b[i-1]*psi[i,j]+(b[i-1]-1.0)*\
                                        ((U[i-1,j]-2*U[i,j]+U[i+1,j])/h/h \
                                         +(phi[i+1,j]-phi[i,j])/h)
            psi[-i-1,j]=b[i-1]*psi[-i-1,j]+(b[i-1]-1.0)*\
                                        ((U[-i-2,j]-2*U[-i-1,j]+U[-i,j])/h/h \
                                         +(phi[-i-1,j]-phi[-i-2,j])/h)
            Up[i,j] += c[i-2,j-2]*((phi[i+1,j]-phi[i,j])/h+psi[i,j])
            Up[-i-1,j] += c[-i-3,j-2]*((phi[-i-1,j]-phi[-i-2,j])/h+psi[-i-1,j])
@njit
def solve_wave_FD2(U,Up,h,c,ncpml,b,psi,phi):
    """ Calculate wave solution using 2nd order finite difference
    """
    ny , nx = U.shape
    for i in range(1,ny-1):
        for j in range(1,nx-1):
            Up[i,j] = 2.0*U[i,j] - Up[i,j] + c[i-1,j-1]* \
                        (U[i-1,j]-2*U[i,j]+U[i+1,j]+ \
                         (U[i,j-1]-2*U[i,j]+U[i,j+1]))/h/h
    #CPML boundary in X-domain
    for i in range(1+ncpml,ny-ncpml-1):
        for j in range(1,ncpml+1):
            phi[i,j]=b[j-1]*phi[i,j]+(b[j-1]-1.0)*(U[i,j+1]-U[i,j])/h
            phi[i,-j-1]=b[j-1]*phi[i,-j-1]+(b[j-1]-1.0)*(U[i,-j-1]-U[i,-j-2])/h
        for j in range(1,ncpml):
            psi[i,j]=b[j-1]*psi[i,j]+(b[j-1]-1.0)*\
                                        ((U[i,j-1]-2*U[i,j]+U[i,j+1])/h/h \
                                         +(phi[i,j+1]-phi[i,j])/h)
            psi[i,-j-1]=b[j-1]*psi[i,-j-1]+(b[j-1]-1.0)*\
                                        ((U[i,-j-2]-2*U[i,-j-1]+U[i,-j])/h/h \
                                         +(phi[i,-j-1]-phi[i,-j-2])/h)
            Up[i,j] += c[i-1,j-1]*((phi[i,j+1]-phi[i,j])/h+psi[i,j])
            Up[i,-j-1] += c[i-1,-j-2]*((phi[i,-j-1]-phi[i,-j-2])/h+psi[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(1,ncpml+1):
        for j in range(1,nx-1):
            phi[i,j]=b[i-1]*phi[i,j]+(b[i-1]-1.0)*(U[i+1,j]-U[i,j])/h
            phi[-i-1,j]=b[i-1]*phi[-i-1,j]+(b[i-1]-1.0)*(U[-i-1,j]-U[-i-2,j])/h
    for i in range(1,ncpml):
        for j in range(1,nx-1):
            psi[i,j]=b[i-1]*psi[i,j]+(b[i-1]-1.0)*\
                                        ((U[i-1,j]-2*U[i,j]+U[i+1,j])/h/h \
                                         +(phi[i+1,j]-phi[i,j])/h)
            psi[-i-1,j]=b[i-1]*psi[-i-1,j]+(b[i-1]-1.0)*\
                                        ((U[-i-2,j]-2*U[-i-1,j]+U[-i,j])/h/h \
                                         +(phi[-i-1,j]-phi[-i-2,j])/h)
            Up[i,j] += c[i-1,j-1]*((phi[i+1,j]-phi[i,j])/h+psi[i,j])
            Up[-i-1,j] += c[-i-2,j-1]*((phi[-i-1,j]-phi[-i-2,j])/h+psi[-i-1,j])

@njit
def solve_wave_stereo(U,Up,Ux,Upx,Uy,Upy,h,c_1,c_2,c_3,ncpml,b,psi,psix,psiy,phi,phix,phiy):
    """ Calculate wave solution using stereo-modeling
        Source : "A central difference method with low numerical dispersion for solving
                    the scalar wave equation "(Dinghui Yang, Ping Tong and Xiaoying Deng,2012)
    """
    ny , nx = U.shape
    for i in range(1,ny-1):
        for j in range(1,nx-1):
            
            d4x = -12.0/(h**4)*(U[i,j-1]\
                                -2*U[i,j]\
                                +U[i,j+1])\
                        +6.0/(h**3)*(Ux[i,j+1]\
                                     -Ux[i,j-1])
    
            d4y = -12.0/(h**4)*(U[i-1,j]\
                                -2*U[i,j]\
                                +U[i+1,j])\
                        +6.0/(h**3)*(Uy[i+1,j]\
                                     -Uy[i-1,j])
            
            d2x2y = 1.0/(h**4)*(2.0*(U[i,j+1]\
                                +U[i,j-1]\
                                +U[i-1,j]\
                                +U[i+1,j]\
                                -2*U[i,j])\
                                -U[i+1,j+1]\
                                -U[i-1,j-1]\
                                -U[i+1,j-1]\
                                -U[i-1,j+1])\
                                    +1.0/(2.0*h**3)*(Ux[i+1,j+1]\
                                                     +Ux[i-1,j+1]\
                                                     -Ux[i-1,j-1]\
                                                     -Ux[i+1,j-1]\
                                                     -2.0*Ux[i,j+1]\
                                                     +2.0*Ux[i,j-1])\
                                    +1.0/(2.0*h**3)*(Uy[i+1,j+1]\
                                                     +Uy[i+1,j-1]\
                                                     -Uy[i-1,j-1]\
                                                     -Uy[i-1,j+1]\
                                                     -2.0*Uy[i+1,j]\
                                                     +2.0*Uy[i-1,j])
            
            Up[i,j] = 2.0*U[i,j] - Up[i,j] + c_1[i-1,j-1]*(U[i-1,j]+U[i+1,j]+U[i,j-1]+U[i,j+1]-4*U[i,j])/(h**2)\
                                                    - c_2[i-1,j-1]*(d4x+d4y) + c_3[i-1,j-1]*d2x2y
                            
            d5x = -90.0/(h**5)*(U[i,j+1]\
                                -U[i,j-1])\
                    +30.0/(h**4)*(Ux[i,j+1]\
                                  +4.0*Ux[i,j]\
                                  +Ux[i,j-1])

            dx4y = -6.0/(h**5)*(U[i+1,j+1]\
                                -U[i-1,j-1]\
                                -U[i+1,j-1]\
                                +U[i-1,j+1]\
                                +2.0*U[i,j-1]\
                                -2.0*U[i,j+1])\
                    +3.0/(h**4)*(Uy[i+1,j+1]\
                                 +Uy[i-1,j-1]\
                                 -Uy[i+1,j-1]\
                                 -Uy[i-1,j+1])

            d3x2y = -3.0/(2.0*h**5)*(U[i+1,j+1]\
                                     -U[i-1,j-1]\
                                     +U[i-1,j+1]\
                                     -U[i+1,j-1]\
                                     +2.0*U[i,j-1]\
                                     -2.0*U[i,j+1])\
                    +3.0/(2.0*h**4)*(Ux[i+1,j+1]\
                                     +Ux[i-1,j-1]\
                                     +Ux[i+1,j-1]\
                                     +Ux[i-1,j+1]\
                                     -2.0*Ux[i,j+1]\
                                     -2.0*Ux[i,j-1])
            
            Upx[i,j] = 2.0*Ux[i,j] - Upx[i,j] + c_1[i-1,j-1]*(Ux[i-1,j]+Ux[i+1,j]+Ux[i,j-1]+Ux[i,j+1]-4*Ux[i,j])/(h**2)\
                        - c_2[i-1,j-1]*(d5x+dx4y)+c_3[i-1,j-1]*d3x2y
                        
            d4xy = -6.0/(h**5)*(U[i+1,j+1]\
                                -U[i-1,j-1]\
                                +U[i+1,j-1]\
                                -U[i-1,j+1]\
                                +2.0*U[i-1,j]\
                                -2.0*U[i+1,j])\
                    +3.0/(h**4)*(Ux[i+1,j+1]\
                                 +Ux[i-1,j-1]\
                                 -Ux[i+1,j-1]\
                                 -Ux[i-1,j+1])
            
            d5y = -90.0/(h**5)*(U[i+1,j]\
                                -U[i-1,j])\
                    +30.0/(h**4)*(Uy[i+1,j]\
                                  +4.0*Uy[i,j]\
                                  +Uy[i-1,j])
            
            d2x3y = -3.0/(2.0*h**5)*(U[i+1,j+1]\
                                     -U[i-1,j-1]\
                                     +U[i+1,j-1]\
                                     -U[i-1,j+1]\
                                     +2.0*U[i-1,j]\
                                     -2.0*U[i+1,j])\
                         +3.0/(2.0*h**4)*(Uy[i+1,j+1]\
                                          +Uy[i-1,j-1]\
                                          +Uy[i+1,j-1]\
                                          +Uy[i-1,j+1]\
                                          -2.0*Uy[i+1,j]\
                                          -2.0*Uy[i-1,j])
            
            Upy[i,j] = 2.0*Uy[i,j] - Upy[i,j] + c_1[i-1,j-1]*(Uy[i-1,j]+Uy[i+1,j]+Uy[i,j-1]+Uy[i,j+1]-4*Uy[i,j])/(h**2)  \
                         - c_2[i-1,j-1]*(d4xy+d5y)+ c_3[i-1][j-1]*(d2x3y)
    '''
    #CPML boundary in X-domain
    for i in range(1+ncpml,ny-ncpml-1):
        for j in range(1,ncpml+1):
            phi[i,j]=b[j-1]*phi[i,j]+(b[j-1]-1.0)*(U[i,j+1]-U[i,j])/h
            phi[i,-j-1]=b[j-1]*phi[i,-j-1]+(b[j-1]-1.0)*(U[i,-j-1]-U[i,-j-2])/h
        for j in range(1,ncpml):
            psi[i,j]=b[j-1]*psi[i,j]+(b[j-1]-1.0)*\
                                        ((U[i,j-1]-2*U[i,j]+U[i,j+1])/h/h \
                                         +(phi[i,j+1]-phi[i,j])/h)
            psi[i,-j-1]=b[j-1]*psi[i,-j-1]+(b[j-1]-1.0)*\
                                        ((U[i,-j-2]-2*U[i,-j-1]+U[i,-j])/h/h \
                                         +(phi[i,-j-1]-phi[i,-j-2])/h)
            Up[i,j] += c_1[i-1,j-1]*((phi[i,j+1]-phi[i,j])/h+psi[i,j])
            Up[i,-j-1] += c_1[i-1,-j-2]*((phi[i,-j-1]-phi[i,-j-2])/h+psi[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(1,ncpml+1):
        for j in range(1,nx-1):
            phi[i,j]=b[i-1]*phi[i,j]+(b[i-1]-1.0)*(U[i+1,j]-U[i,j])/h
            phi[-i-1,j]=b[i-1]*phi[-i-1,j]+(b[i-1]-1.0)*(U[-i-1,j]-U[-i-2,j])/h
    for i in range(1,ncpml):
        for j in range(1,nx-1):
            psi[i,j]=b[i-1]*psi[i,j]+(b[i-1]-1.0)*\
                                        ((U[i-1,j]-2*U[i,j]+U[i+1,j])/h/h \
                                         +(phi[i+1,j]-phi[i,j])/h)
            psi[-i-1,j]=b[i-1]*psi[-i-1,j]+(b[i-1]-1.0)*\
                                        ((U[-i-2,j]-2*U[-i-1,j]+U[-i,j])/h/h \
                                         +(phi[-i-1,j]-phi[-i-2,j])/h)
            Up[i,j] += c_1[i-1,j-1]*((phi[i+1,j]-phi[i,j])/h+psi[i,j])
            Up[-i-1,j] += c_1[-i-2,j-1]*((phi[-i-1,j]-phi[-i-2,j])/h+psi[-i-1,j])
    
    #CPML boundary in X-domain
    for i in range(1+ncpml,ny-ncpml-1):
        for j in range(1,ncpml+1):
            phix[i,j]=b[j-1]*phix[i,j]+(b[j-1]-1.0)*(Ux[i,j+1]-Ux[i,j])/h
            phix[i,-j-1]=b[j-1]*phix[i,-j-1]+(b[j-1]-1.0)*(Ux[i,-j-1]-Ux[i,-j-2])/h
        for j in range(1,ncpml):
            psix[i,j]=b[j-1]*psix[i,j]+(b[j-1]-1.0)*\
                                        ((Ux[i,j-1]-2*Ux[i,j]+Ux[i,j+1])/h/h \
                                         +(phix[i,j+1]-phix[i,j])/h)
            psix[i,-j-1]=b[j-1]*psix[i,-j-1]+(b[j-1]-1.0)*\
                                        ((Ux[i,-j-2]-2*Ux[i,-j-1]+Ux[i,-j])/h/h \
                                         +(phix[i,-j-1]-phix[i,-j-2])/h)
            Upx[i,j] += c_1[i-1,j-1]*((phix[i,j+1]-phix[i,j])/h+psix[i,j])
            Upx[i,-j-1] += c_1[i-1,-j-2]*((phix[i,-j-1]-phix[i,-j-2])/h+psix[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(1,ncpml+1):
        for j in range(1+ncpml,nx-ncpml-1):
            phix[i,j]=b[i-1]*phix[i,j]+(b[i-1]-1.0)*(Ux[i+1,j]-Ux[i,j])/h
            phix[-i-1,j]=b[i-1]*phix[-i-1,j]+(b[i-1]-1.0)*(Ux[-i-1,j]-Ux[-i-2,j])/h
    for i in range(1,ncpml):
        for j in range(1+ncpml,nx-ncpml-1):
            psix[i,j]=b[i-1]*psix[i,j]+(b[i-1]-1.0)*\
                                        ((Ux[i-1,j]-2*Ux[i,j]+Ux[i+1,j])/h/h \
                                         +(phix[i+1,j]-phix[i,j])/h)
            psix[-i-1,j]=b[i-1]*psix[-i-1,j]+(b[i-1]-1.0)*\
                                        ((Ux[-i-2,j]-2*Ux[-i-1,j]+Ux[-i,j])/h/h \
                                         +(phix[-i-1,j]-phix[-i-2,j])/h)
            Upx[i,j] += c_1[i-1,j-1]*((phix[i+1,j]-phix[i,j])/h+psix[i,j])
            Upx[-i-1,j] += c_1[-i-2,j-1]*((phix[-i-1,j]-phix[-i-2,j])/h+psix[-i-1,j])
    
    #CPML boundary in X-domain
    for i in range(1+ncpml,ny-ncpml-1):
        for j in range(1,ncpml+1):
            phiy[i,j]=b[j-1]*phiy[i,j]+(b[j-1]-1.0)*(Uy[i,j+1]-Uy[i,j])/h
            phiy[i,-j-1]=b[j-1]*phiy[i,-j-1]+(b[j-1]-1.0)*(Uy[i,-j-1]-Uy[i,-j-2])/h
        for j in range(1,ncpml):
            psiy[i,j]=b[j-1]*psix[i,j]+(b[j-1]-1.0)*\
                                        ((Uy[i,j-1]-2*Uy[i,j]+Uy[i,j+1])/h/h \
                                         +(phiy[i,j+1]-phiy[i,j])/h)
            psiy[i,-j-1]=b[j-1]*psix[i,-j-1]+(b[j-1]-1.0)*\
                                        ((Uy[i,-j-2]-2*Uy[i,-j-1]+Uy[i,-j])/h/h \
                                         +(phiy[i,-j-1]-phiy[i,-j-2])/h)
            Upy[i,j] += c_1[i-1,j-1]*((phiy[i,j+1]-phiy[i,j])/h+psiy[i,j])
            Upy[i,-j-1] += c_1[i-1,-j-2]*((phiy[i,-j-1]-phiy[i,-j-2])/h+psiy[i,-j-1])
    
    #CPML boundary in Y-domain
    for i in range(1,ncpml+1):
        for j in range(1+ncpml,nx-ncpml-1):
            phiy[i,j]=b[i]*phiy[i,j]+(b[i-1]-1.0)*(Uy[i+1,j]-Uy[i,j])/h
            phiy[-i-1,j]=b[i-1]*phiy[-i-1,j]+(b[i-1]-1.0)*(Uy[-i-1,j]-Uy[-i-2,j])/h
    for i in range(1,ncpml):
        for j in range(1+ncpml,nx-ncpml-1):
            psiy[i,j]=b[i-1]*psiy[i,j]+(b[i-1]-1.0)*\
                                        ((Uy[i-1,j]-2*Uy[i,j]+Uy[i+1,j])/h/h \
                                         +(phiy[i+1,j]-phiy[i,j])/h)
            psiy[-i-1,j]=b[i-1]*psiy[-i-1,j]+(b[i-1]-1.0)*\
                                        ((Uy[-i-2,j]-2*Uy[-i-1,j]+Uy[-i,j])/h/h \
                                         +(phiy[-i-1,j]-phiy[-i-2,j])/h)
            Upy[i,j] += c_1[i-1,j-1]*((phiy[i+1,j]-phiy[i,j])/h+psiy[i,j])
            Upy[-i-1,j] += c_1[-i-2,j-1]*((phiy[-i-1,j]-phiy[-i-2,j])/h+psiy[-i-1,j])
    '''