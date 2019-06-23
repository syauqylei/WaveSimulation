#include <iostream>
#include <iomanip>
#include <fstream>
#include "wvesim.h"
#include <cmath>

wvesim::wvesim(){
	velmod=NULL;
	ext_velmod = NULL;
	Un = NULL;
	Wn = NULL;
	Vn = NULL;
	Wn_p = NULL;
	Vn_p = NULL;
	Un_p = NULL;}

wvesim::wvesim( input fwd_input) {
	nx = fwd_input.Nx;
	ny = fwd_input.Ny;
	
	h = fwd_input.h; 
	
	Nx_ext = 2*nb+2+fwd_input.Nx;
	Ny_ext = 2*nb+2+fwd_input.Ny;
	
	dim = Nx_ext*Ny_ext;
	velmod=fwd_input.velmod;
	ext_velmod = new double[dim];
	Un = new double[dim];
	Wn = new double[dim];
	Vn = new double[dim];
	Wn_p = new double[dim];
	Vn_p = new double[dim];
	Un_p = new double[dim];
	extend_velmod();
	}

wvesim::~wvesim(){
	delete [] Un;
	delete [] Wn;
	delete [] Vn;
	delete [] Wn_p;
	delete [] Vn_p;	
	delete [] Un_p;
	delete [] ext_velmod;
	}

void wvesim::extend_velmod(){
	for (int i=0;i<ny;i++){
		for (int j=0;j<nx;j++){
				ext_velmod[(i+nb+1)*Nx_ext+(j+nb+1)]=velmod[i*nx+j];
			}
		}
	
	double *y= new double[25];
	double *x= new double[25];
	double *f= new double[25];
	
	//this block of codes fills top region of velmod 
	for(int i=0;i<nb+1;i++){
		for(int j=0;j<nx;j++){
			double yi=double(i);
			double xi=double(j);
			
			for(int k=0;k<5;k++){
				for(int l=0;l<5;l++){
					int some_var;
					if(j==0){some_var=0;}
					else if(j==1){some_var=0;}
					else if(j==nx-1){some_var=nx-5;}
					else if(j==nx-2){some_var=nx-5;}
					else {some_var=j-2;}
					
					x[k*5+l]=double((some_var+l));
					y[k*5+l]=double((k+nb+1));
					f[k*5+l]=ext_velmod[(k+nb+1)*Nx_ext+(some_var+l+nb+1)];
				}
			}		
			ext_velmod[i*Nx_ext+(j+nb+1)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}

	//this block of codes fills bottom region of velmod 
	for(int i=nb+1+ny;i<Ny_ext;i++){
		for(int j=0;j<nx;j++){
			double yi=double(i);
			double xi=double(j);
			
			for(int k=0;k<5;k++){
				for(int l=0;l<5;l++){
					
					int some_var;
			
					if(j==0){some_var=0;}
					else if(j==1){some_var=0;}
					else if(j==nx-1){some_var=nx-5;}
					else if(j==nx-2){some_var=nx-5;}
					else {some_var=j-2;}
					
					x[k*5+l]=double((some_var+l));
					y[k*5+l]=double((nb+1+ny-1-k));
					f[k*5+l]=ext_velmod[(nb+1+ny-1-k)*Nx_ext+(some_var+l+nb+1)];
				}
			}		
			ext_velmod[i*Nx_ext+(j+nb+1)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}
	
	//this block of codes fills left region of velmod 
	for(int i=0;i<Ny_ext;i++){
		for(int j=0;j<nb+1;j++){
			double yi=double(i);
			double xi=double(j);
			
			for(int k=0;k<5;k++){
				for(int l=0;l<5;l++){
					
					int some_var;
			
					if(i==0){some_var=0;}
					else if(i==1){some_var=0;}
					else if(i==Ny_ext-1){some_var=Ny_ext-5;}
					else if(i==Ny_ext-2){some_var=Ny_ext-5;}
					else {some_var=i-2;}
					
					x[k*5+l]=double((nb+1+l));
					y[k*5+l]=double((some_var+k));
					f[k*5+l]=ext_velmod[(some_var+k)*Nx_ext+(l+nb+1)];
				}
			}		
			ext_velmod[i*Nx_ext+(j)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}
	
	//this block of codes fills right region of velmod 
	for(int i=0;i<Ny_ext;i++){
		for(int j=nb+1+nx;j<Nx_ext;j++){
			double yi=double(i);
			double xi=double(j);
			
			for(int k=0;k<5;k++){
				for(int l=0;l<5;l++){
					
					int some_var;
			
					if(i==0){some_var=0;}
					else if(i==1){some_var=0;}
					else if(i==Ny_ext-1){some_var=Ny_ext-5;}
					else if(i==Ny_ext-2){some_var=Ny_ext-5;}
					else {some_var=i-2;}
					
					x[k*5+l]=double((nb+1+nx-l-1));
					y[k*5+l]=double((some_var+k));
					f[k*5+l]=ext_velmod[(some_var+k)*Nx_ext+(nb+1+nx-l-1)];
				}
			}		
			ext_velmod[i*Nx_ext+(j)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}

	delete [] y;
	delete [] f;
	delete [] x;
	}

double wvesim::iwd_interp(double xi, double yi ,double *y, double *x, double *f, int len){
	double fx;
	double dist;
	double up_var=0;
	double bot_var=0;
	
	for (int i=0;i<len;i++){
		dist=sqrt((xi-x[i])*(xi-x[i])+(yi-y[i])*(yi-y[i]));
		up_var+=1.0/dist*f[i];
		bot_var+=1.0/dist;
		}

	fx=up_var/bot_var;
	return fx;
	}
