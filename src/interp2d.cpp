#include "interp2d.h"
#include <iostream>
#include <cmath>

interp2d::interp2d(){
	}
	
interp2d::interp2d(std::vector<double> xi,std::vector<double> yi, int ord, 
					std::vector<double> x, std::vector<double> y, std::vector<double> f){
	x_raw = x;
	y_raw = y;
	f_raw = f;
	order = ord;
	int len = xi.size();
	double check_dist;
	fi = new double[len];	
	std::cout<<"size of Velmod"<< len<<std::endl;
	for (int i=0;i<len;i++){
		std::cout<<i <<" ";
		check_dist=euc_dist(xi[i],yi[i],x[i],y[i]);
		if(check_dist  == 0){continue;}
		fi[i]=interp(xi[i],yi[i]);
		}
	}


double interp2d::interp(double xx, double yy){
	double fx;
	double Upr,Bot;
	//searchin nearest neighbors
	cal_dist(xx,yy);
	sort_dist();
	
	for (int i=0;i<order;i++){
		Upr+=f_dist[i].second/euc_dist(xx,yy,x_dist[i].second,y_dist[i].second);
		Bot+=1/euc_dist(xx,yy,x_dist[i].second,y_dist[i].second);
		}
	clear_vect();
	fx=Upr/Bot;
	return fx;
	}

double interp2d::euc_dist(double x, double y, double xx, double yy){
	return sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy));
	}

void interp2d::cal_dist(double xx, double yy){
	int len = x_raw.size();
	
	for (int i=0; i<len;i++){
		x_dist.push_back(std::make_pair(sqrt((xx-x_raw[i])*(xx-x_raw[i])+(yy-y_raw[i])*(yy-y_raw[i])),x_raw[i]));
		y_dist.push_back(std::make_pair(sqrt((xx-x_raw[i])*(xx-x_raw[i])+(yy-y_raw[i])*(yy-y_raw[i])),y_raw[i]));
		f_dist.push_back(std::make_pair(sqrt((xx-x_raw[i])*(xx-x_raw[i])+(yy-y_raw[i])*(yy-y_raw[i])),f_raw[i]));
		}
		
	}

void interp2d::sort_dist(){
	
	if(x_dist.empty() == true){throw(Exception("Vector Container of x_distance","Vector is empty"));}
	if(y_dist.empty() == true){throw(Exception("Vector Container of f_distance","Vector is empty"));}
	if(f_dist.empty() == true){throw(Exception("Vector Container of f_distance","Vector is empty"));}
	
	std::sort(x_dist.begin(),x_dist.end());
	std::sort(y_dist.begin(),y_dist.end());
	std::sort(f_dist.begin(),f_dist.end());
	}

void interp2d::clear_vect(){
	x_dist.clear();
	y_dist.clear();
	f_dist.clear();
	}

  
interp2d::~interp2d(){
	delete [] fi;
	}
