#ifndef INPUT_H
#define INPUT_H

#include "interp1d.h"
#include "interp2d.h"
#include "error.h"

const double nb = 32;

class input: public interp1d {
	public :
	int Nx;
	int Ny;
	int Nx_ext;
	int Ny_ext;
	int Nt;
	
	double dt;
	double h;
	
	//soure location
	int xloc;
	int yloc;
	
	
	double *velmod;
	double *srcfunc;
	
	input();
	input(std::string filename);
	~input();
	
	private:
	std::vector<double> t_src;
	std::vector<double> t_src_raw;
	std::vector<double> srcfunc_raw;
	std::vector<double> velmod_raw;
	
	std::vector<double> xx_vel;
	std::vector<double> yy_vel;
	std::vector<double> xx_vel_ext;
	std::vector<double> yy_vel_ext;
	
	std::string f_velmod;
	std::string f_srcfunc;
	std::string f_out;
	
	void read_parms(std::string fname);
	void read_srcfunc();
	void read_velmod();
};

#endif
