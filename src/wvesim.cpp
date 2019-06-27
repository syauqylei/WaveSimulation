#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "wvesim.h"

wvesim::wvesim(){
	velmod=NULL;
	ext_velmod = NULL;
	Un = NULL;
	Wn = NULL;
	Vn = NULL;}

wvesim::wvesim( input fwd_input) {
	nx = fwd_input.Nx;
	ny = fwd_input.Ny;
	
	h = fwd_input.h;
	
	Nx = 2*nb+2+fwd_input.Nx;
	Ny = 2*nb+2+fwd_input.Ny;
	NX = 2*nb+fwd_input.Nx;
	NY = 2*nb+fwd_input.Nx;
	dim = Nx*Ny;
	velmod=fwd_input.velmod;
	
	vmax=*std::max_element(velmod,velmod+nx*ny);		
	dt = cfl*h/vmax;
	nt = int(fwd_input.T/dt);
	std::cout<<" dt = "<<dt <<std::endl;
	ext_velmod = new double[dim];
	Un = alloc_array(3, dim);
	Wn = alloc_array(3, dim);
	Vn = alloc_array(3, dim);
	
	//auxilaries var for CPML
	//left
	psix_u_lf = alloc_array(Ny, nb);
	phix_u_lf = alloc_array(Ny, (nb+1));
	psix_w_lf = alloc_array(Ny, nb);
	phix_w_lf = alloc_array(Ny, (nb+1));
	psix_v_lf = alloc_array(Ny, nb);
	phix_v_lf = alloc_array(Ny, (nb+1));
	
	//right
	psix_u_rt = alloc_array(Ny, nb);
	phix_u_rt = alloc_array(Ny, (nb+1));
	psix_w_rt = alloc_array(Ny, nb);
	phix_w_rt = alloc_array(Ny, (nb+1));
	psix_v_rt = alloc_array(Ny, nb);
	phix_v_rt = alloc_array(Ny, (nb+1));
	
	
	//top
	psix_u_tp = alloc_array(Nx, nb);
	phix_u_tp = alloc_array(Nx, (nb+1));
	psix_w_tp = alloc_array(Nx, nb);
	phix_w_tp = alloc_array(Nx, (nb+1));
	psix_v_tp = alloc_array(Nx, nb);
	phix_v_tp = alloc_array(Nx, (nb+1));
	
	//bottom
	psix_u_bt = alloc_array(Nx, nb);
	phix_u_bt = alloc_array(Nx, (nb+1));
	psix_w_bt = alloc_array(Nx, nb);
	phix_w_bt = alloc_array(Nx, (nb+1));
	psix_v_bt = alloc_array(Nx, nb);
	phix_v_bt = alloc_array(Nx, (nb+1));
	
	//extend_velmod
	extend_velmod();
	init_cpml_parm();
	
	//initiate values to array container
	for(int i=0;i<Ny;i++){
		for(int j=0;j<nb;j++){
			psix_u_lf[i][j]=0.0;
			psix_w_lf[i][j]=0.0;
			psix_v_lf[i][j]=0.0;
			
			psix_u_rt[i][j]=0.0;
			psix_w_rt[i][j]=0.0;
			psix_v_rt[i][j]=0.0;
		}
		for(int j=0;j<nb+1;j++){
			phix_v_rt[i][j]=0.0;
			phix_w_rt[i][j]=0.0;
			phix_u_rt[i][j]=0.0;
			
			phix_u_lf[i][j]=0.0;
			phix_w_lf[i][j]=0.0;
			phix_v_lf[i][j]=0.0;
		}
	}
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<nb;j++){
			psix_u_tp[i][j]=0.0;
			psix_w_tp[i][j]=0.0;
			psix_v_tp[i][j]=0.0;
			
			psix_u_bt[i][j]=0.0;
			psix_w_bt[i][j]=0.0;
			psix_v_bt[i][j]=0.0;
			}
		for (int j=0;j<nb+1;j++){
			phix_v_tp[i][j]=0.0;
			phix_u_tp[i][j]=0.0;
			phix_w_tp[i][j]=0.0;
			
			phix_v_bt[i][j]=0.0;
			phix_u_bt[i][j]=0.0;
			phix_w_bt[i][j]=0.0;
			}
		}
		
	for (int i=0;i<dim;i++){
		Un[0][i]=0;
		Wn[0][i]=0;
		Vn[0][i]=0;
		Wn[1][i]=0;
		Vn[1][i]=0;
		Un[1][i]=0;
		Wn[2][i]=0;
		Vn[2][i]=0;
		Un[2][i]=0;
		}
	
	//creating index to calculate from
	idx =  new int[(NX)*(NY)];
	for (int i=1;i<(Ny-1);i++){
		for (int j=1;j<(Nx-1);j++){
			idx[(i-1)*(NX)+(j-1)]=i*(Nx)+j;
			}
		}
		
	//creating approriate source function by interp1d
	src_func = new double[nt];
	src_t = new double[nt];
	
	interp1d src_interp(3,fwd_input.t_src,fwd_input.srcfunc);
	for (int i=0 ;i<nt;i++){
		src_t[i]=i*dt;
		src_func[i]=src_interp.interp(i*dt);
		}
	
	init_cpml_parm();
	
	src_loc =(nb+1+fwd_input.yloc)*(Nx)+(nb+1+fwd_input.xloc);
	record= alloc_array(nt,nx);
	for (int i=0;i<nt-1;i++){
		if((i+1)%(nt/10) == 0) {std::cout<<"Calculating Wavefields .... "<< std::setw(3)<<std::setprecision(1) << double(i+2)/double(nt)*100 <<"%\n";}
		Un[1][src_loc]=src_func[i];
		for (int j=0;j<nx;j++){
			record[i][j]=Un[1][(nb+1)*Nx+nb+1+j];
		}
		wve_calc();
		wve_update();
		/*
		std::string name="output";
		std::stringstream ss;
		ss<<i;
		std::string fname=name+ss.str()+".txt";
		write_txt(fname);*/
		}
	
	for (int j=0;j<nx;j++){
			record[nt-1][j]=Un[2][(nb+1)*Nx+nb+1+j];
		}
		
	write_rec(fwd_input.f_out);
	}

wvesim::~wvesim(){
	delete [] idx;
	delete [] ext_velmod;
	delete [] velmod;
	delete [] b;
	delete [] src_func;
	delete [] src_t;
	
	free_array_mem(Un);
	free_array_mem(Wn);
	free_array_mem(Vn);
	free_array_mem(psix_u_lf);
	free_array_mem(psix_w_lf);
	free_array_mem(psix_v_lf);
	free_array_mem(phix_u_lf);
	free_array_mem(phix_w_lf);
	free_array_mem(phix_v_lf);
	free_array_mem(psix_u_rt);
	free_array_mem(psix_w_rt);
	free_array_mem(psix_v_rt);
	free_array_mem(phix_u_rt);
	free_array_mem(phix_w_rt);
	free_array_mem(phix_v_rt);
	free_array_mem(psix_u_tp);
	free_array_mem(psix_w_tp);
	free_array_mem(psix_v_tp);
	free_array_mem(phix_u_tp);
	free_array_mem(phix_w_tp);
	free_array_mem(phix_v_tp);
	free_array_mem(psix_u_bt);
	free_array_mem(psix_w_bt);
	free_array_mem(psix_v_bt);
	free_array_mem(phix_u_bt);
	free_array_mem(phix_w_bt);
	free_array_mem(phix_v_bt);
	}

void wvesim::init_cpml_parm(){
	b=new double[nb+1];
	
	double L=nb*h;

	for (int i=0;i<nb+1;i++){
		double du;
		du=-3.0*vmax/2/L*log(R)*(L-i*h)*(L-i*h)/L/L;
		b[i]=exp(-du*dt);
		}
	}

void wvesim::wve_calc(){
	for (int i=0;i<(NX)*(NY);i++){
			int index=idx[i];
			double c1=ext_velmod[index]*ext_velmod[index]*dt*dt;
			double c2=(ext_velmod[index]*ext_velmod[index]
						*h*h*dt*dt
						-ext_velmod[index]*ext_velmod[index]
						*ext_velmod[index]*ext_velmod[index]
						*dt*dt*dt*dt)/12.0;
			double c3=ext_velmod[index]*ext_velmod[index]
						*ext_velmod[index]*ext_velmod[index]
						*dt*dt*dt*dt/6.0;
						
			Un[2][index]=2.0*Un[1][index]-Un[0][index]
						+c1*d2xd2y(Un[1],index)
						-c2*d4x(Un[1],Wn[1],index)-c2*d4y(Un[1],Vn[1],index)
						+c3*d2x2y(Un[1],Wn[1],Vn[1],index);
			
			Wn[2][index]=2.0*Wn[1][index]-Wn[0][index]
						+c1*d2xd2y(Wn[1],index)
						-c2*d5x(Un[1],Wn[1],index)-c2*dx4y(Un[1],Vn[1],index)
						+c3*d3x2y(Un[1],Wn[1],index);
				
			Vn[2][index]=2.0*Vn[1][index]-Vn[0][index]
						+c1*d2xd2y(Vn[1],index)
						-c2*d4xy(Un[1],Wn[1],index)-c2*d5y(Un[1],Vn[1],index)
						+c3*d2x3y(Un[1],Vn[1],index);
		}
	
	left_cpml();
	right_cpml();
	top_cpml();
	bottom_cpml();
	}

void wvesim::left_cpml(){
	//left_cmpl
	for (int i=1;i<Ny-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_u_lf[i][j]=b[j]*phix_u_lf[i][j]+(b[j]-1.0)*(Un[1][i*Nx+j+2]-Un[1][i*Nx+j+1])/h;
				}
			for(int j=0;j<nb;j++){
				psix_u_lf[i][j]=b[j]*psix_u_lf[i][j]+(b[j]-1.0)
							*((Un[1][i*Nx+j+2]-2.0*Un[1][i*Nx+j+1]+Un[1][i*Nx+j])/h/h
								+(phix_u_lf[i][j+1]-phix_u_lf[i][j])/h);
								
				Un[2][i*Nx+j+1]+=ext_velmod[i*Nx+j+1]*ext_velmod[i*Nx+j+1]*dt*dt
									*(psix_u_lf[i][j]+(phix_u_lf[i][j+1]-phix_u_lf[i][j])/h);
				}
		}
	for (int i=1;i<Ny-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_w_lf[i][j]=b[j]*phix_w_lf[i][j]+(b[j]-1.0)*(Wn[1][i*Nx+j+2]-Wn[1][i*Nx+j+1])/h;
				}
			for(int j=0;j<nb;j++){
				psix_w_lf[i][j]=b[j]*psix_w_lf[i][j]+(b[j]-1.0)
							*((Wn[1][i*Nx+j+2]-2.0*Wn[1][i*Nx+j+1]+Wn[1][i*Nx+j])/h/h
								+(phix_w_lf[i][j+1]-phix_w_lf[i][j])/h);
								
				Wn[2][i*Nx+j+1]+=ext_velmod[i*Nx+j+1]*ext_velmod[i*Nx+j+1]*dt*dt
									*(psix_w_lf[i][j]+(phix_w_lf[i][j+1]-phix_w_lf[i][j])/h);
				}
		}
	for (int i=1;i<Ny-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_v_lf[i][j]=b[j]*phix_v_lf[i][j]+(b[j]-1.0)*(Vn[1][i*Nx+j+2]-Vn[1][i*Nx+j+1])/h;
				}
			for(int j=0;j<nb;j++){
				psix_v_lf[i][j]=b[j]*psix_v_lf[i][j]+(b[j]-1.0)
							*((Vn[1][i*Nx+j+2]-2.0*Vn[1][i*Nx+j+1]+Vn[1][i*Nx+j])/h/h
								+(phix_v_lf[i][j+1]-phix_v_lf[i][j])/h);
								
				Vn[2][i*Nx+j+1]+=ext_velmod[i*Nx+j+1]*ext_velmod[i*Nx+j+1]*dt*dt
									*(psix_v_lf[i][j]+(phix_v_lf[i][j+1]-phix_v_lf[i][j])/h);
				}
		}
	}

void wvesim::right_cpml(){
	//right_cmpl
	for (int i=1;i<Ny-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_u_rt[i][nb-j]=b[j]*phix_u_rt[i][nb-j]+(b[j]-1.0)*(Un[1][i*Nx+Nx-j-2]-Un[1][i*Nx+Nx-j-3])/h;
				}
			for(int j=0;j<nb;j++){
				psix_u_rt[i][nb-1-j]=b[j]*psix_u_rt[i][nb-1-j]+(b[j]-1.0)
							*((Un[1][i*Nx+Nx-j-3]-2.0*Un[1][i*Nx+Nx-j-2]+Un[1][i*Nx+Nx-1-j])/h/h
								+(phix_u_rt[i][nb-j]-phix_u_rt[i][nb-1-j])/h);
								
				Un[2][i*Nx+Nx-j-2]+=ext_velmod[i*Nx+Nx-j-2]*ext_velmod[i*Nx+Nx-j-2]*dt*dt
									*(psix_u_rt[i][nb-1-j]+(phix_u_rt[i][nb-j]-phix_u_rt[i][nb-1-j])/h);
				}
		}
	for (int i=1;i<Ny-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_w_rt[i][nb-j]=b[j]*phix_w_rt[i][nb-j]+(b[j]-1.0)*(Wn[1][i*Nx+Nx-j-2]-Wn[1][i*Nx+Nx-j-3])/h;
				}
			for(int j=0;j<nb;j++){
				psix_w_rt[i][nb-1-j]=b[j]*psix_w_rt[i][nb-1-j]+(b[j]-1.0)
							*((Wn[1][i*Nx+Nx-j-3]-2.0*Wn[1][i*Nx+Nx-j-2]+Wn[1][i*Nx+Nx-1-j])/h/h
								+(phix_w_rt[i][nb-j]-phix_w_rt[i][nb-1-j])/h);
								
				Un[2][i*Nx+Nx-j-2]+=ext_velmod[i*Nx+Nx-j-2]*ext_velmod[i*Nx+Nx-j-2]*dt*dt
									*(psix_w_rt[i][nb-1-j]+(phix_w_rt[i][nb-j]-phix_w_rt[i][nb-1-j])/h);
				}
		}
	for (int i=1;i<Ny-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_v_rt[i][nb-j]=b[j]*phix_v_rt[i][nb-j]+(b[j]-1.0)*(Vn[1][i*Nx+Nx-j-2]-Vn[1][i*Nx+Nx-j-3])/h;
				}
			for(int j=0;j<nb;j++){
				psix_v_rt[i][nb-1-j]=b[j]*psix_v_rt[i][nb-1-j]+(b[j]-1.0)
							*((Vn[1][i*Nx+Nx-j-3]-2.0*Vn[1][i*Nx+Nx-j-2]+Vn[1][i*Nx+Nx-1-j])/h/h
								+(phix_v_rt[i][nb-j]-phix_v_rt[i][nb-1-j])/h);
								
				Vn[2][i*Nx+Nx-j-2]+=ext_velmod[i*Nx+Nx-j-2]*ext_velmod[i*Nx+Nx-j-2]*dt*dt
									*(psix_v_rt[i][nb-1-j]+(phix_v_rt[i][nb-j]-phix_v_rt[i][nb-1-j])/h);
				}
		}
	}

void wvesim::top_cpml(){
	//top
	for (int i=1;i<Nx-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_u_tp[i][j]=b[j]*phix_u_tp[i][j]+(b[j]-1)*(Un[1][(j+2)*Nx+i]-Un[1][(j+1)*Nx+i])/h;
				}
			for(int j=0;j<nb;j++){
				psix_u_tp[i][j]=b[j]*psix_u_tp[i][j]+(b[j]-1.0)
							*((Un[1][(j+2)*Nx+i]-2.0*Un[1][(j+1)*Nx+i]+Un[1][(j)*Nx+i])/h/h
								+(phix_u_tp[i][j+1]-phix_u_tp[i][j])/h);
								
				Un[2][(j+1)*Nx+i]+=ext_velmod[(j+1)*Nx+i]*ext_velmod[(j+1)*Nx+i]*dt*dt
									*(psix_u_tp[i][j]+(phix_u_tp[i][j+1]-phix_u_tp[i][j])/h);
				}
		}
	
	for (int i=1;i<Nx-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_w_tp[i][j]=b[j]*phix_w_tp[i][j]+(b[j]-1)*(Wn[1][(j+2)*Nx+i]-Wn[1][(j+1)*Nx+i])/h;
				}
			for(int j=0;j<nb;j++){
				psix_w_tp[i][j]=b[j]*psix_w_tp[i][j]+(b[j]-1.0)
							*((Wn[1][(j+2)*Nx+i]-2.0*Wn[1][(j+1)*Nx+i]+Wn[1][(j)*Nx+i])/h/h
								+(phix_w_tp[i][j+1]-phix_w_tp[i][j])/h);
								
				Wn[2][(j+1)*Nx+i]+=ext_velmod[(j+1)*Nx+i]*ext_velmod[(j+1)*Nx+i]*dt*dt
									*(psix_w_tp[i][j]+(phix_w_tp[i][j+1]-phix_w_tp[i][j])/h);
				}
		}
	
	for (int i=1;i<Nx-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_v_tp[i][j]=b[j]*phix_v_tp[i][j]+(b[j]-1)*(Vn[1][(j+2)*Nx+i]-Vn[1][(j+1)*Nx+i])/h;
				}
			for(int j=0;j<nb;j++){
				psix_v_tp[i][j]=b[j]*psix_v_tp[i][j]+(b[j]-1.0)
							*((Vn[1][(j+2)*Nx+i]-2.0*Vn[1][(j+1)*Nx+i]+Vn[1][(j)*Nx+i])/h/h
								+(phix_v_tp[i][j+1]-phix_v_tp[i][j])/h);
								
				Vn[2][(j+1)*Nx+i]+=ext_velmod[(j+1)*Nx+i]*ext_velmod[(j+1)*Nx+i]*dt*dt
									*(psix_v_tp[i][j]+(phix_v_tp[i][j+1]-phix_v_tp[i][j])/h);
				}
		}
	}
void wvesim::bottom_cpml(){
	//bottom
	for (int i=1;i<Nx-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_u_bt[i][nb-j]=b[j]*phix_u_bt[i][nb-j]+(b[j]-1)*(Un[1][Nx*Ny-(j+2)*Nx+i]-Un[1][Nx*Ny-(j+3)*Nx+i])/h;
				}
			for(int j=0;j<nb;j++){
				psix_u_bt[i][nb-1-j]=b[j]*psix_u_bt[i][nb-1-j]+(b[j]-1.0)
							*((Un[1][Nx*Ny-(j+3)*Nx+i]-2.0*Un[1][Nx*Ny-(j+2)*Nx+i]+Un[1][Nx*Ny-(j+1)*Nx+i])/h/h
							+(phix_u_bt[i][nb-j]-phix_u_bt[i][nb-1-j])/h);
								
				Un[2][Nx*Ny-(j+2)*Nx+i]+=ext_velmod[Nx*Ny-(j+2)*Nx+i]*ext_velmod[Nx*Ny-(j+2)*Nx+i]*dt*dt
									*(psix_u_bt[i][nb-1-j]+(phix_u_bt[i][nb-j]-phix_u_bt[i][nb-j-1])/h);
				}
		}
	
	for (int i=1;i<Nx-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_w_bt[i][nb-j]=b[j]*phix_w_bt[i][nb-j]+(b[j]-1)*(Wn[1][Nx*Ny-(j+2)*Nx+i]-Wn[1][Nx*Ny-(j+3)*Nx+i])/h;
				}
			for(int j=0;j<nb;j++){
				psix_w_bt[i][nb-1-j]=b[j]*psix_w_bt[i][nb-1-j]+(b[j]-1.0)
							*((Wn[1][Nx*Ny-(j+3)*Nx+i]-2.0*Wn[1][Nx*Ny-(j+2)*Nx+i]+Wn[1][Nx*Ny-(j+1)*Nx+i])/h/h
								+(phix_w_bt[i][nb-j]-phix_w_bt[i][nb-1-j])/h);
								
				Wn[2][Nx*Ny-(j+2)*Nx+i]+=ext_velmod[Nx*Ny-(j+2)*Nx+i]*ext_velmod[Nx*Ny-(j+2)*Nx+i]*dt*dt
									*(psix_w_bt[i][nb-1-j]+(phix_w_bt[i][nb-j]-phix_w_bt[i][nb-j-1])/h);
				}
		}
	
	for (int i=1;i<Nx-1;i++){
			for(int j=0;j<nb+1;j++){
				phix_v_bt[i][nb-j]=b[j]*phix_v_bt[i][nb-j]+(b[j]-1)*(Vn[1][Nx*Ny-(j+2)*Nx+i]-Vn[1][Nx*Ny-(j+3)*Nx+i])/h;
				}
			for(int j=0;j<nb;j++){
				psix_v_bt[i][nb-1-j]=b[j]*psix_v_bt[i][nb-1-j]+(b[j]-1.0)
							*((Vn[1][Nx*Ny-(j+3)*Nx+i]-2.0*Vn[1][Nx*Ny-(j+2)*Nx+i]+Vn[1][Nx*Ny-(j+1)*Nx+i])/h/h
								+(phix_v_bt[i][nb-j]-phix_v_bt[i][nb-1-j])/h);
								
				Vn[2][Nx*Ny-(j+2)*Nx+i]+=ext_velmod[Nx*Ny-(j+2)*Nx+i]*ext_velmod[Nx*Ny-(j+2)*Nx+i]*dt*dt
									*(psix_v_bt[i][nb-1-j]+(phix_v_bt[i][nb-j]-phix_v_bt[i][nb-j-1])/h);
				}
		}
	
	}

void wvesim::wve_update(){
		for (int i=0;i<dim;i++){
		Un[0][i]=Un[1][i];
		Wn[0][i]=Wn[1][i];
		Vn[0][i]=Vn[1][i];
		Un[1][i]=Un[2][i];
		Wn[1][i]=Wn[2][i];
		Vn[1][i]=Vn[2][i];
		}
	}

void wvesim::extend_velmod(){
	for (int i=0;i<ny;i++){
		for (int j=0;j<nx;j++){
				ext_velmod[(i+nb+1)*Nx+(j+nb+1)]=velmod[i*nx+j];
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
					f[k*5+l]=ext_velmod[(k+nb+1)*Nx+(some_var+l+nb+1)];
				}
			}		
			ext_velmod[i*Nx+(j+nb+1)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}

	//this block of codes fills bottom region of velmod 
	for(int i=nb+1+ny;i<Ny;i++){
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
					f[k*5+l]=ext_velmod[(nb+1+ny-1-k)*Nx+(some_var+l+nb+1)];
				}
			}		
			ext_velmod[i*Nx+(j+nb+1)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}
	
	//this block of codes fills left region of velmod 
	for(int i=0;i<Ny;i++){
		for(int j=0;j<nb+1;j++){
			double yi=double(i);
			double xi=double(j);
			
			for(int k=0;k<5;k++){
				for(int l=0;l<5;l++){
					
					int some_var;
			
					if(i==0){some_var=0;}
					else if(i==1){some_var=0;}
					else if(i==Ny-1){some_var=Ny-5;}
					else if(i==Ny-2){some_var=Ny-5;}
					else {some_var=i-2;}
					
					x[k*5+l]=double((nb+1+l));
					y[k*5+l]=double((some_var+k));
					f[k*5+l]=ext_velmod[(some_var+k)*Nx+(l+nb+1)];
				}
			}		
			ext_velmod[i*Nx+(j)]=iwd_interp(xi,yi,y,x,f,25);
			}
		}
	
	//this block of codes fills right region of velmod 
	for(int i=0;i<Ny;i++){
		for(int j=nb+1+nx;j<Nx;j++){
			double yi=double(i);
			double xi=double(j);
			
			for(int k=0;k<5;k++){
				for(int l=0;l<5;l++){
					
					int some_var;
			
					if(i==0){some_var=0;}
					else if(i==1){some_var=0;}
					else if(i==Ny-1){some_var=Ny-5;}
					else if(i==Ny-2){some_var=Ny-5;}
					else {some_var=i-2;}
					
					x[k*5+l]=double((nb+1+nx-l-1));
					y[k*5+l]=double((some_var+k));
					f[k*5+l]=ext_velmod[(some_var+k)*Nx+(nb+1+nx-l-1)];
				}
			}		
			ext_velmod[i*Nx+(j)]=iwd_interp(xi,yi,y,x,f,25);
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

double wvesim::dxf(double *u,int i){
	return
	(u[i+1]-u[i])/h;}

double wvesim::dxb(double *u,int i){
	return
	(u[i]-u[i-1])/h;}

double wvesim::dyf(double *u,int i){
	return
	(u[i+Nx]-u[i])/h;} 

double wvesim::dyb(double *u,int i){
	return
	(u[i]-u[i-Nx])/h;}

double wvesim::d2x(double *u,int i){
	return
	(u[i-1]-2.0*u[i]+u[i+1])/h/h;
	}

double wvesim::d2y(double *u,int i){
	return
	(u[i-Nx]-2.0*u[i]+u[i+Nx])/h/h;
	}

double wvesim::d2xd2y(double *u,int i){
	return 
	(u[i+1]+u[i-1]+u[i+Nx]+u[i-Nx]-4.0*u[i])/h/h;
	}
	
double wvesim::d2x2y(double *u,double *u_x,double *u_y,int i){
	return
	(2.0*(u[i+1]+u[i-1]+u[i+Nx]+u[i-Nx]-2.0*u[i])
		-u[i+1+Nx]-u[i-1-Nx]-u[i-1+Nx]-u[i+1-Nx])/h/h/h/h
	+(u_x[i+1+Nx]+u_x[i+1-Nx]-u_x[i-1-Nx]-u_x[i-1+Nx]
		-2.0*u_x[i+1]+2.0*u_x[i-1])/h/h/h/2.0
	+(u_y[i+1+Nx]+u_y[i-1+Nx]-u_y[i-1+Nx]-u_y[i-1-Nx]
		-2.0*u_y[i+Nx]+2.0*u_y[i-Nx])/h/h/h/2.0;
	}

double wvesim::d4x(double *u, double *u_x,int i){
	return
	-12.0/h/h/h/h*(u[i+1]+u[i-1]-2.0*u[i])
	+6.0/h/h/h*(u_x[i+1]-u_x[i-1]);
	}

double wvesim::d4y(double *u, double *u_y,int i){
	return
	-12.0/h/h/h/h*(u[i+Nx]+u[i-Nx]-2.0*u[i])
	+6.0/h/h/h*(u_y[i+Nx]-u_y[i-Nx]);
	}

double wvesim::d3x2y(double *u, double *u_x,int i){
	return
	-3.0/2.0/h/h/h/h/h*(u[i+Nx+1]-u[i-1-Nx]+u[i+1-Nx]
						-u[i-1+Nx]+2.0*u[i-1]-2.0*u[i+1])
	+3.0/2.0/h/h/h/h*(u_x[i+1+Nx]+u_x[i-1-Nx]+u_x[i-1+Nx]
						+u_x[i+1-Nx]-2.0*u_x[i+1]-2.0*u_x[i-1]);
	}

double wvesim::d5x(double *u, double *u_x,int i){
	return
	-90.0/h/h/h/h/h*(u[i+1]-u[i-1])
	+30.0/h/h/h/h*(u_x[i+1]+4.0*u_x[i]+u_x[i-1]);
	}

double wvesim::dx4y(double *u, double *u_y,int i){
	return
	-6.0/h/h/h/h/h*(u[i+Nx+1]-u[i-1-Nx]-u[i-1+Nx]
					+u[i+1-Nx]+2.0*u[i-1]-2.0*u[i+1])
	+3.0/h/h/h/h*(u_y[i+1+Nx]+u_y[i-1-Nx]-u_y[i-1+Nx]-u_y[i+1-Nx]);
}
	
double wvesim::d2x3y(double *u, double *u_y,int i){
	return
	-3.0/2.0/h/h/h/h/h*(u[i+Nx+1]-u[i-1-Nx]+u[i-1+Nx]
						-u[i+1-Nx]+2.0*u[i-Nx]-2.0*u[i+Nx])
	+3.0/2.0/h/h/h/h*(u_y[i+1+Nx]+u_y[i-1-Nx]+u_y[i-1+Nx]
						+u_y[i+1-Nx]-2.0*u_y[i+Nx]-2.0*u_y[i-Nx]);
	}

double wvesim::d4xy(double *u,double *u_x,int i){
	return
	-6.0/h/h/h/h/h*(u[i+Nx+1]-u[i-1-Nx]+u[i-1+Nx]
					-u[i+1-Nx]+2.0*u[i-Nx]-2.0*u[i+Nx])
	+3.0/h/h/h/h*(u_x[i+1+Nx]+u_x[i-1-Nx]-u_x[i-1+Nx]-u_x[i+1-Nx]);
	}

double wvesim::d5y(double *u,double *u_y,int i){
	return
	-90.0/h/h/h/h/h*(u[i+Nx]-u[i-Nx])
	+30.0/h/h/h/h*(u_y[i+Nx]+4.0*u_y[i]+u_y[i-Nx]);
	}

void wvesim::write_txt(std::string f_name){
	std::ofstream file(f_name);
	for(int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			file << std::setprecision(3)<<std::setw(5)<< Un[2][(nb+1+i)*Nx+(nb+1+j)] << "\t";
			}
			file << std::endl;
		}
	file.close();
	}

void wvesim::write_rec(std::string f_name){
	std::ofstream file(f_name);
	for(int i=0;i<nt;i++){
		for(int j=0;j<nx;j++){
			file << std::setprecision(3)<<std::setw(5)<< record[i][j] << "\t";
			}
			file << std::endl;
		}
	file.close();
	}

double **wvesim::alloc_array(const int nrows, const int ncols){
	double **mat=new double*[nrows];
	
	mat[0]= new double[nrows*ncols];
	for (int i=1; i<nrows;i++){
		mat[i]=&mat[0][i*ncols];
	}	
	return mat;
}

void wvesim::free_array_mem(double **mat){
	delete [] mat[0];
	delete [] mat;
}
