#ifndef INTERPOLATION_H
#define INTERPOLATION_H
class interp1d {
	public :
	
	interp1d(double *x, double *f, int nx, int order);
	interp1d(double x);
	interp1d(double *x);
	~interp1d();

	private :
	
	double *x_dat;
	double *f_dat;
	
	int n;
	int loc;
	int order;
	
	double weight_func_lgr(double xi, int order);
	double interp_lgr(double xi, int order);
	int index_xi(double xi);
	
	};

class interp2d {
	
	};
#endif
