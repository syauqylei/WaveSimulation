#ifndef INPUT_H
#define INPUT_H

class input {
	public :
	int Nx;
	int Ny;
	int Nt;
	
	double dt;
	double h;
	
	//soure location
	int xloc;
	int yloc;
	
	
	double *velmod;
	double *srcfunc;
	double *t_src;
	
	input();
	input(std::string filename);
	~input();
	
	private:
	std::string f_velmod;
	std::string f_srcfunc;
	std::string f_out;
	
	void read_parms(std::string fname);
	void read_srcfunc();
	void read_velmod();
};

#endif
