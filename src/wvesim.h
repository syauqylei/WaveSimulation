#ifndef WVESIM_H
#define WVESIM_H

#include "input.h"
#include <vector>
const int nb = 32;

class wvesim {
	
	public:
	wvesim();
	wvesim(input fwd_input);
	~wvesim();

	private:
	int dim;
	int Nx_ext;
	int Ny_ext;
	int nx;
	int ny;
	
	double h;
	
	double *ext_velmod;
	double *velmod;
	double *Un;
	double *Wn;
	double *Vn;
	double *Wn_p;
	double *Vn_p;
	double *Un_p;
		
	void extend_velmod();
	double iwd_interp(double xi, double yi ,double *y, double *x, double *f, int len);
	};

#endif
