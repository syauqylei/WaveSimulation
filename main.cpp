#include <iostream>
#include <cassert>
#include "input.h"

int main(int* argc, char** argv){
	std::string filename;
	filename = argv[1];
	input parms(filename);
	return 0;
}
