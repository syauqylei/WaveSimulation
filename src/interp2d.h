#ifndef INTERPOLATION2_H
#define INTERPOLATION2_H

#include <vector>
#include <algorithm>
#include <utility>
#include "error.h"

class interp2d {
	public :
	
	interp2d();
	interp2d(std::vector<double> x, std::vector<double> y, std::vector<double> f, int ord);
	~interp2d();
	
	double interp(double xx, double yy);

	private :
	int order;
	std::vector<double> x_raw;
	std::vector<double> y_raw;
	std::vector<double> f_raw;
	
	std::vector< std::pair<double,double>> x_dist;
	std::vector< std::pair<double,double>> y_dist;
	std::vector< std::pair<double,double>> f_dist;
	
	double euc_dist(double x, double y, double xx, double yy);
	void cal_dist(double xx, double yy);
	void sort_dist();
	void clear_vect();
	
	};

#endif
