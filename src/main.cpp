#include <iostream>
#include <cassert>
#include "input.h"
#include "wvesim.h"

int main(int* argc, char** argv){
	std::string filename;
	filename = argv[1];
	input parms(filename);
	wvesim propag(parms);
	return 0;
}
