#include <string>
#include <fstream>
#include <iostream>
#include <vector>
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
	
	velmod = new double[Nx*Ny];
	srcfunc = new double[Nt];
	t_src = new double[Nt];
	
	read_srcfunc();
	read_velmod(); //read velmod
	
	}

void input::read_parms(std::string filename){
	std::ifstream f(filename);
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
	std::ifstream f(f_srcfunc);
	int i=0;
	while (f >> t_src[i] >> srcfunc[i]){
		i++;
		}
	f.close();
	}

void input::read_velmod(){
	std::ifstream f(f_velmod);
	int i=0;
	while (f >> velmod[i]){i++;}
	f.close();
	}

input::~input(){
	delete [] srcfunc;
	delete [] velmod;
	delete [] t_src;
	}
