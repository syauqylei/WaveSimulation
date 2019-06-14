#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <algorithm>
#include <utility>
#include "error.h"

class interp1d {
	public :
	double *fi;
	
	interp1d();
	interp1d(std::vector<double> xi, int ord, std::vector<double> x, std::vector<double> f);
	~interp1d();

	private :
	int order;
	std::vector<double> x_raw;
	std::vector<double> f_raw;
	
	std::vector< std::pair<double,double>> x_dist;
	std::vector< std::pair<double,double>> f_dist;
	
	double interp(double xx);
	void cal_dist(double xx);
	void sort_dist();
	void clear_vect();
	
	};

#endif
