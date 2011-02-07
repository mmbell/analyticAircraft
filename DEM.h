/*
 *  DEM.h
 *  analyticAircraft
 *
 *  Created by Michael Bell on 2/4/11.
 *  Copyright 2011. All rights reserved.
 *
 */

#ifndef DEM_H
#define DEM_H

#include "geotiff.h"
#include "geo_normalize.h"

class DEM  
{

public:
	DEM();
	~DEM();
	
	bool readDem(char* fname);
	
private:
	int GTIFReportACorner( GTIF *gtif, GTIFDefn *defn, FILE * fp_out,
						const char * corner_name,
						double x, double y, int inv_flag, int dec_flag );
	void GTIFPrintCorners( GTIF *, GTIFDefn *, FILE *, int, int, int, int );
};


#endif
