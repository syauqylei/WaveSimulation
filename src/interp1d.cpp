#include "interp1d.h"
#include <iostream>
#include <cmath>

interp1d::interp1d(){
	}
	
interp1d::interp1d(std::vector<double> xi, int ord, std::vector<double> x, std::vector<double> f){
	x_raw = x;
	f_raw = f;
	order = ord;
	int len = xi.size();
	
	fi = new double[len];	
	for (int i=0;i<len;i++){
		fi[i]=interp(xi[i]);
		}
	}


double interp1d::interp(double xx){
	double fx=0;
	
	//searchin nearest neighbors
	cal_dist(xx);
	sort_dist();
	
	for (int i=0;i<order;i++){
		double L=1;
		for (int j=0;j<order;j++){
			
			if (i==j){continue;}
			
			L*=(xx-x_dist[j].second)/(x_dist[i].second-x_dist[j].second);
			}
			
			fx+=L*f_dist[i].second;
		}
	
	clear_vect();
	
	return fx;
	}

void interp1d::cal_dist(double xx){
	int len = x_raw.size();
	
	for (int i=0; i<len;i++){
		x_dist.push_back(std::make_pair(sqrt((xx-x_raw[i])*(xx-x_raw[i])),x_raw[i]));
		f_dist.push_back(std::make_pair(sqrt((xx-x_raw[i])*(xx-x_raw[i])),f_raw[i]));
		}
		
	}

void interp1d::sort_dist(){
	
	if(x_dist.empty() == true){throw(Exception("Vector Container of x_distance","Vector is empty"));}
	if(f_dist.empty() == true){throw(Exception("Vector Container of f_distance","Vector is empty"));}
	
	std::sort(x_dist.begin(),x_dist.end());
	std::sort(f_dist.begin(),f_dist.end());
	}

void interp1d::clear_vect(){
	x_dist.clear();
	f_dist.clear();
	}

  
interp1d::~interp1d(){
	delete [] fi;
	}
