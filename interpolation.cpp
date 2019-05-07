#include <cmath>

interp1d::interp1d(double *x, double *f, int nx){
	x_dat = x;
	f_dat = f;
	n = nx;
	loc = 0;
	}

double interp1d::interp_lgr(double xi, int order){
	
	}

double interp1d::weight_func_lgr(double xi, int order){
	index_xi(xi);//give loc a value
	int mid_index = order/2;
	double L;
	for (int k = 0; k < order; k++){
		int i=loc-mid_index+k;
		int j=loc-mid_index+1+k;
		L*= xi-x_dat[j]/(x_dat[i]-x_dat[j]);
	return L;
	}

//search the appropriate point to interpolate from
int interp1d::index_xi(double xi){
	loc = 0;
	while (xi < x_dat[i]){
		loc++;
		}
	}

interp1d::~interp1d(){
	delete [] x_dat;
	delete [] f_dat;
	}
