#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "error.h"
#include "input.h"

input::input(){
	Nx=0;
	Ny=0;
	Nx_ext=0;
	Ny_ext=0;
	Nt=0;
	xloc=0;
	yloc=0;
	dt=0;
	h=0;
	srcfunc=NULL;
	velmod=NULL;
	}

input::input(std::string filename){
	
	std::cout<< "Reading Forward Parameters \n";
	read_parms(filename);//read parameters non-array
	
	std::cout<< "Reading Source Function File \n";
	srcfunc = new double[Nt];
	for (int i=0;i<Nt;i++){ t_src.push_back( dt*i);}
	read_srcfunc();
	std::cout<< "Reading Source Function File  is done\n";
	
	std::cout<< "Extending model \n";
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			xx_vel.push_back(j);
			yy_vel.push_back(i);
			}
		}	
	for (int i=-nb;i<Ny+nb;i++){
		for(int j=-nb;j<Nx+nb;j++){
			xx_vel_ext.push_back(j);
			yy_vel_ext.push_back(i);
			}
		}
	//extending domain
	Nx_ext=2*nb+Nx;
	Ny_ext=2*nb+Ny;
	
	velmod = new double[Nx_ext*Ny_ext];
	std::cout<< "Velocity Model is being read \n";
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
	std::ifstream v_file;
	
	v_file.open(f_velmod);
	
	if (v_file.is_open() == false) {throw (Exception("FILE","File doesn't exist or File can't be opened"));}
	
	int i=0;
	double val;
	while (v_file >> val){velmod_raw.push_back(val);i++;}
	v_file.close();
	interp2d intrp_vel(xx_vel_ext,yy_vel_ext, 13, 
						xx_vel,yy_vel, velmod_raw);
	
	for (int i=0;i<Ny_ext;i++){
		for(int j=0 ;j<Nx_ext;j++){
			velmod[i*Nx_ext+j]=intrp_vel.fi[i*Nx_ext+j];
			}
		}
}

input::~input(){
	delete [] srcfunc;
	delete [] velmod;
	}
