#include <iostream>
#include "input.h"

wvesim::wvesim(){
	Un=NULL;
	Wn=NULL;
	Vn=NULL;
	Wn_2=NULL;
	Wn_1=NULL;
	Vn_2=NULL;
	Vn_1=NULL;
	Un_2=NULL;
	Un_1=NULL;
	}

wvesim::wvesim( input parms_obj) const {
	
	dim = parms_obj.Nx*parms_obj.Ny;
	
	Un = new double[dim];
	Wn = new double[dim];
	Vn = new double[dim];
	Wn_2 = new double[dim];
	Wn_1 = new double[dim];
	Vn_2 = new double[dim];
	Vn_1 = new double[dim];
	Un_2 = new double[dim];
	Un_1 = new double[dim];
	}


wvesim::~wvesim(){
	delete [] Un;
	delete [] Wn;
	delete [] Vn;
	delete [] Wn_2;
	delete [] Wn_1;
	delete [] Vn_2;
	delete [] Vn_1;
	delete [] Un_2;
	delete [] Un_1;
	}
