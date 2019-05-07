#ifndef WVESIM_H
#define WVESIM_H

class wvesim {
	
	public:

	double *Un;
	
	wvesim();
	~wvesim();

	private:
	
	double *Wn;
	double *Vn;
	double *Wn_2;
	double *Wn_1;
	double *Vn_2;
	double *Vn_1;
	double *Un_2;
	double *Un_1;
	
	extend_velmod();
	
	};

#endif
