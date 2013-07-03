/*
 *  WRF.h
 *
 *  Created by Michael Bell on 5/29/13.
 *  Based on example code from netCDF library
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef WRF_H
#define WRF_H

#include <netcdfcpp.h>

class WRF  
{

public:
	WRF();
	~WRF();
	
	bool readWRF(const char* filename);
	bool getData(const double &lat, const double &lon,const double &alt, double &u_out, double &v_out, double &w_out, double &dbz_out);
	
private:
	int NDIMS, NLVL, NLAT, NLON, NREC;

    float* lats;
	float* lons;	
	float* u;
	float* v;
	float* w;
	float* dbz;
	float* z;
	
	float* uwind_in;
	float* vwind_in;
	float* wwind_in;
	float* dbz_in;
	float* z_in;
	float* ph_in;
	float* phb_in;
};


#endif
