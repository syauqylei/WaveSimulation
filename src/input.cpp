#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "error.h"
#include "input.h"

input::input(){
	Nx=0;
	Ny=0;
	Nt=0;
	xloc=0;
	yloc=0;
	dt=0;
	h=0;
	srcfunc=NULL;
	velmod=NULL;
	}

input::input(std::string filename){
	
	read_parms(filename);//read parameters non-array
	srcfunc = new double[Nt];
	for (int i=0;i<Nt;i++){ t_src.push_back( dt*i);}
	read_srcfunc();
	
	velmod = new double[Nx*Ny];
	read_velmod(); //read velmod
	
	}

void input::read_parms(std::string filename){
	std::ifstream f(filename);
	
	if (f.is_open() == false) {throw (Exception("FILE","File doesn't exist or File can't be opened"));}
	
	int i=0;
	std::string content;
	std::string::size_type sz; 
	while (f >> content){
		if (i==1){f_velmod = content ;}
		if (i==3){Nx = std::stoi(content,&sz);}
		if (i==5){Ny = std::stoi(content,&sz);}
		if (i==7){Nt = std::stoi(content,&sz);}
		if (i==9){dt = std::stod(content,&sz);}
		if (i==11){h = std::stod(content,&sz);}
		if (i==13){f_srcfunc = content ;}
		if (i==15){xloc = std::stoi(content,&sz);}
		if (i==16){yloc = std::stoi(content,&sz);}
		if (i==18){f_out = content ;}
		i++;
		}
	f.close();
	}

void input::read_srcfunc(){
	std::ifstream file_src(f_srcfunc);	
	
	if (file_src.is_open() == false) {throw (Exception("FILE","File doesn't exist or File can't be opened"));}
	
	int i=0;
	double a;
	double b;
	while (file_src >> a >> b){
		t_src_raw.push_back(a);
		srcfunc_raw.push_back(b);
		i++;
		}
	file_src.close();
	
	//interpolate
	interp1d intrp(t_src, 5, t_src_raw, srcfunc_raw);
	
	for (int i=0;i<Nt;i++){
		srcfunc[i]=intrp.fi[i];
		}
	}

void input::read_velmod(){
	std::ifstream v_file(f_velmod);
	
	if (v_file.is_open() == false) {throw (Exception("FILE","File doesn't exist or File can't be opened"));}
	
	int i=0;
	double val;
	while (v_file >> val){velmod[i]=val;i++;}
	v_file.close();
}

input::~input(){
	delete [] srcfunc;
	delete [] velmod;
	}
