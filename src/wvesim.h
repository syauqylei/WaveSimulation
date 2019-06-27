#ifndef WVESIM_H
#define WVESIM_H

#include "input.h"
#include <vector>
#include <algorithm>
const int nb = 40;
const double R = 0.001;

class wvesim {
	
	public:
	double **record;
	wvesim();
	wvesim(input fwd_input);
	~wvesim();

	private:
	int dim;
	int NX;
	int NY;
	int Nx;
	int Ny;
	int nx;
	int ny;
	int nt;
	int src_loc;
	
	double h;
	double dt;
	
	int *idx;
	
	double *ext_velmod;
	double *velmod;
	
	void extend_velmod();
	double iwd_interp(double xi, double yi ,double *y, double *x, double *f, int len);
	
	double d2x(double *u,int i);
	double d2y(double *u,int i);
	double dxf(double *u,int i);
	double dxb(double *u,int i);
	double dyf(double *u,int i);
	double dyb(double *u,int i);
	
	double d2xd2y(double *u,int i);
	
	double d2x2y(double *u,double *u_x,double *u_y,int i);
	double d4x(double *u, double *u_x,int i);
	double d4y(double *u, double *u_y,int i);
	
	double d3x2y(double *u, double *u_x,int i);
	double d5x(double *u, double *u_x,int i);
	double dx4y(double *u, double *u_y,int i);
	
	double d2x3y(double *u, double *u_y,int i);
	double d4xy(double *u,double *u_x,int i);
	double d5y(double *u,double *u_y,int i);
	
	//cpml paramters
	double *b;
	double vmax;
	void init_cpml_parm();
	
	double **Un;
	double **Wn;
	double **Vn;
	void wve_calc();
	void wve_update();
	
	double **psix_u_lf;
	double **phix_u_lf;
	double **psix_w_lf;
	double **phix_w_lf;
	double **psix_v_lf;
	double **phix_v_lf;
	void left_cpml();
	
	double **psix_u_rt;
	double **phix_u_rt;
	double **psix_w_rt;
	double **phix_w_rt;
	double **psix_v_rt;
	double **phix_v_rt;
	void right_cpml();
	
	double **psix_u_tp;
	double **phix_u_tp;
	double **psix_w_tp;
	double **phix_w_tp;
	double **psix_v_tp;
	double **phix_v_tp;
	void top_cpml();
	
	
	double **psix_u_bt;
	double **phix_u_bt;
	double **psix_w_bt;
	double **phix_w_bt;
	double **psix_v_bt;
	double **phix_v_bt;
	void bottom_cpml();
	
	void write_txt(std::string f_name);
	void write_rec(std::string f_name);
	
	double **alloc_array(const int nrows, const int ncols);
	void free_array_mem(double **mat);
	};

#endif
