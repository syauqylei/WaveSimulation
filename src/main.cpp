#include <iostream>
#include <cassert>
#include "input.h"
#include "wvesim.h"

int main(int* argc, char** argv){
	std::string filename;
	filename = argv[1];
	input parms(filename);
	std::cout<<"is this error1?"<<std::endl;
	wvesim propag(parms);
	std::cout<<"is this error2?"<<std::endl;
	return 0;
}
