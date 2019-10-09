#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 15:49:11 2019

@author: rinn
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from scipy import interpolate
from .solve_wave import *



class wave:
    
    def __init__(self,h,T,Vel):
        
        self.list_of_methods = {'stereo':1,'FD2': 1,'FD4' : 2,'FD8': 4, 'FD12' : 6,'FD16':8}
        #the chosen method to simulate
        self.method = "FD2"
        
        self.h = h
        self.T = T
        self.cfl = 0.50
        
        #adding boundary cells
        self.Vel = Vel
        self.ny , self.nx = self.Vel.shape
        
        #Calculate dt 
        MaxVel = max(map(max,Vel))
        self.dt = self.cfl*self.h/MaxVel
        self.n = math.ceil(self.T/self.dt)
        #CPML Boundary parameters
        self.n_cmpl = 32
        self.R=0.0001
        L=self.n_cmpl*self.h
        self.dmp_par = [math.exp(-(-3.0*MaxVel/2/L*math.log(self.R)*\
                        (L-i*self.h)*(L-i*self.h)/L/L)*self.dt) for i in range(self.n_cmpl+1)]
        
        # default ricker wavelet with 30 Hz frequency  
        # Source : https://wiki.seg.org/wiki/Dictionary:Ricker_wavelet
        self.src = np.array([(1.0-2.0*(math.pi*30.0*i*self.dt)**2) \
                             *math.exp(-(math.pi*30.0*i*self.dt)**2) \
                             for i in range(self.n)])    
        self.ext_cells = self.list_of_methods.get(self.method)
        #set the default location to the center of the array
        self.loc = self.ext_cells+self.n_cmpl+self.ny//2 , \
                    self.ext_cells+self.n_cmpl+self.nx//2
        self.n_iter = 0
        self.i_show = self.n

    def set_loc(self,i,j):
        """ Set the location of the source
        """
        self.ext_cells = self.list_of_methods.get(self.method)
        self.loc = self.n_cmpl+self.ext_cells+i , self.n_cmpl+self.ext_cells+j
    
    
    def set_i_show(self,i):
        """ Set n_show for diplaying the wavefield in every iteration 
            of a multiple of i
        """
        self.i_show = i
    
    def set_src_func(src_func):
        """ Set the source function that is used in the simulation
            src_func[0,:] is value of f(x)
            src_func[1,:] is value of x
        """
        t=np.array([ i*self.dt for i in range(n)])
        
        """ interpolate the source fucntion so it has the same calculated time step
            from the stability condiiton
        """
        f = interpolate.interp1d(src_func[0,:],src_func[1,:])
        self.src = f(t)
    
    def set_method(self,meth):
        """ Set a method to run the simulation. The input argument is a string.
            ex : 'FD2' or 'FD4'
        """
        self.method = meth
    
   
    def disp(self,i_show = 0):
        """ Display the wavefield map every a multiple of n-iteration
        """
        if i_show != 0  and self.n_iter % i_show == 0:
            plt.imshow(self.U,cmap="seismic",vmax=max(self.src),vmin=-max(self.src))
            plt.show()
    
    def extend_velocity(self):
        """ Extend velocity due to CMPL Boundary condition
        """
        
        self.ext_cells = self.list_of_methods.get(self.method)
        x = np.linspace(0,self.h*self.nx,self.nx)
        y = np.linspace(0,self.h*self.ny,self.ny)
        
        xnew = np.linspace(-self.h*self.n_cmpl,self.h*self.nx+self.h*self.n_cmpl,self.nx+2*self.n_cmpl)
        ynew = np.linspace(-self.h*self.n_cmpl,self.h*self.ny+self.h*self.n_cmpl,self.ny+2*self.n_cmpl)
        
        f = interpolate.RectBivariateSpline(y,x,self.Vel,kx=5,ky=5)
        
        self.ext_vel = f(ynew,xnew)
    
    def extend_wave_array(self):
        """ Extend wavefields array to approriate shape 
            for solve_wave module operation
        """
         
        self.extend_velocity()
        
        self.record = np.zeros([self.n,self.nx])
        self.U = np.zeros_like(self.ext_vel)
        self.Up = np.zeros_like(self.ext_vel)
        n_ext = self.list_of_methods.get(self.method)
        ny , nx = self.U.shape
        
        ext_arr_y = np.zeros([ny,2*n_ext])
        ext_arr_x = np.zeros([2*n_ext,nx+2*n_ext])

        if self.method != 'stereo' :
            self.U = np.concatenate((ext_arr_y,self.U),axis=1)
            self.U = np.concatenate((ext_arr_x,self.U),axis=0)

            self.phi = np.zeros_like(self.U)
            self.psi = np.zeros_like(self.U)

            self.Up = np.concatenate((ext_arr_y,self.Up),axis=1)
            self.Up = np.concatenate((ext_arr_x,self.Up),axis=0)
        
        if self.method == 'stereo' :
            self.Ux = np.zeros_like(self.ext_vel)
            self.Uy = np.zeros_like(self.ext_vel)
            self.Upx = np.zeros_like(self.ext_vel)
            self.Upy = np.zeros_like(self.ext_vel)
            
            self.U = np.concatenate((ext_arr_y,self.U),axis=1)
            self.U = np.concatenate((ext_arr_x,self.U),axis=0)
            self.Up = np.concatenate((ext_arr_y,self.Up),axis=1)
            self.Up = np.concatenate((ext_arr_x,self.Up),axis=0)
            
            self.Ux = np.concatenate((ext_arr_y,self.Ux),axis=1)
            self.Ux = np.concatenate((ext_arr_x,self.Ux),axis=0)
            self.Upx = np.concatenate((ext_arr_y,self.Upx),axis=1)
            self.Upx = np.concatenate((ext_arr_x,self.Upx),axis=0)
            
            self.Uy = np.concatenate((ext_arr_y,self.Uy),axis=1)
            self.Uy = np.concatenate((ext_arr_x,self.Uy),axis=0)
            self.Upy = np.concatenate((ext_arr_y,self.Upy),axis=1)
            self.Upy = np.concatenate((ext_arr_x,self.Upy),axis=0)
            
            self.phi = np.zeros_like(self.U)
            self.psi = np.zeros_like(self.U)
            self.phix = np.zeros_like(self.U)
            self.psix = np.zeros_like(self.U)
            self.phiy = np.zeros_like(self.U)
            self.psiy = np.zeros_like(self.U)
            
    def run_fd(self):
        """ run the wave simulation
        """
        self.extend_wave_array()
        c1 = self.ext_vel*self.ext_vel*self.dt*self.dt
        start = timer()
        for i in range(self.n):
            self.n_iter += 1
            
            self.U[self.loc] += self.src[i]
            
            if self.method == 'FD2':
                solve_wave_FD2(self.U,self.Up,self.h,c1,self.n_cmpl,self.dmp_par,self.psi,self.phi)
            if self.method == 'FD4':
                solve_wave_FD4(self.U,self.Up,self.h,c1,self.n_cmpl,self.dmp_par,self.psi,self.phi)
            if self.method == 'FD8':
                solve_wave_FD8(self.U,self.Up,self.h,c1,self.n_cmpl,self.dmp_par,self.psi,self.phi)
            if self.method == 'FD12':
                solve_wave_FD12(self.U,self.Up,self.h,c1,self.n_cmpl,self.dmp_par,self.psi,self.phi)
            if self.method == 'FD16':
                solve_wave_FD16(self.U,self.Up,self.h,c1,self.n_cmpl,self.dmp_par,self.psi,self.phi)

            #update
            ptr=self.Up
            self.Up=self.U
            self.U=ptr
            self.record[i] = self.U[self.n_cmpl+self.ext_cells,self.n_cmpl+self.ext_cells:-self.n_cmpl-self.ext_cells]
            self.disp(self.i_show)
        
        end = timer()
        
        print("Runtime ",end-start ,"s")

    def run_stereo(self):
        """ run the wave simulation with stereo-modeling
        """
        self.extend_wave_array()
        c1 = (self.ext_vel*self.dt)**2
        c2 = ((self.ext_vel*self.dt*self.h)**2-(self.ext_vel*self.dt)**4)/12.0
        c3 = (self.ext_vel*self.dt)**4/6.0


        start = timer()
        for i in range(self.n):
            
            self.U[self.loc] += self.src[i]
            solve_wave_stereo(self.U,self.Up,self.Ux,self.Upx,self.Uy,self.Upy,self.h,c1,c2,c3
                              ,self.n_cmpl,self.dmp_par,self.psi,self.psix,self.psiy,self.phi,self.phix,self.phiy)
            
            #update
            ptr=self.Up
            self.Up=self.U
            self.U=ptr
            
            ptrx=self.Upx
            self.Upx=self.Ux
            self.Ux=ptrx
            
            ptry=self.Upy
            self.Upy=self.Uy
            self.Uy=ptry
            self.disp(self.i_show)
            self.n_iter += 1
            
            self.record[i] = self.U[self.n_cmpl+self.ext_cells,self.n_cmpl+1:-self.n_cmpl-1]
        end = timer()
        
        print("Runtime ",end-start ,"s")
