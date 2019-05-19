#include <cmath>

interp1d::interp1d(double *x, double *f, int nx, int order){
	x_dat = x;
	f_dat = f;
	n = nx;
	loc = 0;
	order = ord;

	}

void interp1d::src_loc(double xi){
	while (xi < x_dat[loc] ){
	loc++;

	neigs = new int [order] ;
	int start;
	if ( loc <= order/2) {start =0 ;}
	else if (loc+order/2 >= n-1) { start = n-1-order;}
	else { start = loc-order/2;}

	for (int i=0;i<order;i++){
		neigs[i]=start+i;}
	
	}

double interp(double xi){
	double fx=1;
	double Lu=1;
	double Lb=1;

	for (int i=0; i<order ; i++){
		int ii=neigs[i];
		for (int j =0 ; j<order;j++){
			int jj=neigs[jj];
			Lu*=(xi-x_dat[jj]);
			Lb*=(x_dat[ii]-x_dat[jj]);
			}
		fx*=Lu/Lb*f_dat[ii];
		}

	return fx;
	}

interp1d::~interp1d(){
	delete [] x_dat;
	delete [] f_dat;
	delete [] neigs;
	}
