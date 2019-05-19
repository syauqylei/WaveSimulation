#ifndef INTERPOLATION_H
#define INTERPOLATION_H

class interp1d {
	public :
	double *fi;
	
	interp1d(double *xi, double *x, double *f, int nx, int ord);
	~interp1d();

	private :
	
	double *x_dat;
	double *f_dat;
	
	int n;
	int loc;
	int order;
	int *neigs; //arrays of neighbor points

	double L;
	
	double weg(double xi);
	double interp(double xi);
	
	//function to search loc in array
	void src_loc(double xi, double *x, int order);
	};

#endif
