#ifndef INPUT_H
#define INPUT_H

#include "interp1d.h"
#include "error.h"
class input: public interp1d {
	public :
	int Nx;
	int Ny;
	
	double dt;
	double Ts;
	double h;
	double T;
	
	//soure location
	int xloc;
	int yloc;
	
	double *velmod;
	std::vector<double> srcfunc;
	std::vector<double> t_src;
	
	std::string f_out;
	
	input();
	input(std::string filename);
	~input();
	
	private:
	
	std::string f_velmod;
	std::string f_srcfunc;
	
	void read_parms(std::string fname);
	void read_srcfunc();
	void read_velmod();
};

#endif
