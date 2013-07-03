/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is an example which reads some 4D pressure and
   temperatures. The data file read by this program is produced by the
   companion program pres_temp_4D_wr.cpp. It is intended to illustrate
   the use of the netCDF C++ API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C++ API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-cxx

   $Id: pres_temp_4D_rd.cpp,v 1.11 2006/08/22 19:22:06 ed Exp $
*/

#include "WRF.h"
#include <iostream>
#include <cmath>

WRF::WRF() {
	NDIMS = 4;
	NLVL = 27;
	NLAT = 249;
	NLON = 249;
	NREC = 1;
	
	lats = new float[NLAT*NLON];
	lons= new float[NLAT*NLON];	
	u= new float[NLON*NLAT*NLVL];
	v= new float[NLON*NLAT*NLVL];
	w= new float[NLON*NLAT*NLVL];
	dbz= new float[NLON*NLAT*NLVL];
	z= new float[NLON*NLAT*NLVL];

	uwind_in= new float[NLVL*NLAT*(NLON+1)];
	vwind_in= new float[NLVL*(NLAT+1)*NLON];
	wwind_in= new float[(NLVL+1)*NLAT*NLON];
	dbz_in= new float[NLVL*NLAT*NLON];
	ph_in= new float[(NLVL+1)*NLAT*NLON];
	phb_in= new float[(NLVL+1)*NLAT*NLON]; 

}

WRF::~WRF() {

	delete[] u;
	delete[] v;
	delete[] w;
	delete[] dbz;
	delete[] z;
	delete[] uwind_in;
	delete[] vwind_in;
	delete[] wwind_in;
	delete[] dbz_in;
	delete[] z_in;
	delete[] ph_in;
	delete[] phb_in;
	delete[] lats;
	delete[] lons;

}

bool WRF::readWRF(const char* filename) {

	// We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
	// timesteps of data.

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;

   // These arrays will hold the data we will read in. We will only
   // need enough space to hold one timestep of data; one record.
   /* float uwind_in[NLVL][NLAT][NLON+1];
   float vwind_in[NLVL][NLAT+1][NLON];
   float wwind_in[NLVL+1][NLAT][NLON];
   float dbz_in[NLVL][NLAT][NLON];
   float ph_in[NLVL+1][NLAT][NLON];
   float phb_in[NLVL+1][NLAT][NLON]; */
   
   // Change the error behavior of the netCDF C++ API by creating an
   // NcError object. Until it is destroyed, this NcError object will
   // ensure that the netCDF C++ API returns error codes on any
   // failure, prints an error message, and leaves any other error
   // handling to the calling program. In the case of this example, we
   // just exit with an NC_ERR error code.
   NcError err(NcError::verbose_nonfatal);

   // Open the file.
   NcFile dataFile(filename, NcFile::ReadOnly);

   // Check to see if the file was opened.
   if(!dataFile.is_valid())
      return NC_ERR;

   // Get pointers to the latitude and longitude variables.
   NcVar *latVar, *lonVar;
   if (!(latVar = dataFile.get_var("XLAT")))
      return NC_ERR;
   if (!(lonVar = dataFile.get_var("XLONG")))
      return NC_ERR;

   // Get pointers to the pressure and temperature variables.
   NcVar *uwindVar, *vwindVar, *wwindVar, *dbzVar, *phVar, *phbVar;
   if (!(uwindVar = dataFile.get_var("U")))
      return NC_ERR;
   if (!(vwindVar  = dataFile.get_var("V")))
      return NC_ERR;
   if (!(wwindVar = dataFile.get_var("W")))
      return NC_ERR;
   if (!(dbzVar  = dataFile.get_var("REFL_10CM")))
      return NC_ERR;
   if (!(phVar  = dataFile.get_var("PH")))
      return NC_ERR;
   if (!(phbVar  = dataFile.get_var("PHB")))
      return NC_ERR;
   
   // Read the data. Since we know the contents of the file we know
   // that the data arrays in this program are the correct size to
   // hold one timestep. 
   for (int rec = 0; rec < NREC; rec++)
   {
      // Read the data one record at a time.
      if (!latVar->set_cur(rec, 0, 0))
 	 return NC_ERR;
      if (!lonVar->set_cur(rec, 0, 0))
 	 return NC_ERR;
      if (!uwindVar->set_cur(rec, 0, 0, 0))
	 return NC_ERR;
      if (!vwindVar->set_cur(rec, 0, 0, 0))
	 return NC_ERR;
      if (!wwindVar->set_cur(rec, 0, 0, 0))
	 return NC_ERR;
      if (!dbzVar->set_cur(rec, 0, 0, 0))
	 return NC_ERR;
      if (!phVar->set_cur(rec, 0, 0, 0))
	 return NC_ERR;
      if (!phbVar->set_cur(rec, 0, 0, 0))
	 return NC_ERR;
	 
      // Get the lat/lon data from the file.
      if (!latVar->get(lats, 1, NLAT, NLON))
	   return NC_ERR;
      if (!lonVar->get(lons, 1, NLAT, NLON))
       return NC_ERR;
	 
      // Get 1 record of NLVL by NLAT by NLON values for each variable.
      if (!uwindVar->get(uwind_in, 1, NLVL, NLAT, NLON+1))
	 return NC_ERR;
      if (!vwindVar->get(vwind_in, 1, NLVL, NLAT+1, NLON))
	 return NC_ERR;
      if (!wwindVar->get(wwind_in, 1, NLVL+1, NLAT, NLON))
	 return NC_ERR;
      if (!dbzVar->get(dbz_in, 1, NLVL, NLAT, NLON))
	 return NC_ERR;
      if (!phVar->get(ph_in, 1, NLVL+1, NLAT, NLON))
	 return NC_ERR;
      if (!phbVar->get(phb_in, 1, NLVL+1, NLAT, NLON))
	 return NC_ERR;
   } // next record 

   for (int i = 0; i < NLON; i++) {
	   for (int j = 0; j < NLAT; j++) {
	      for (int k = 0; k < NLVL; k++) {
			  u[k*NLON*NLAT + j*NLON + i] = (uwind_in[k*(NLON*1)*NLAT + j*(NLON+1)+ i] + uwind_in[k*(NLON+1)*NLAT + j*(NLON+1) + i+1])/2.0;
			  v[k*NLON*NLAT + j*NLON + i] = (vwind_in[k*NLON*(NLAT+1) + j*NLON + i] + vwind_in[k*NLON*(NLAT+1) + (j+1)*NLON + i])/2.0;
			  w[k*NLON*NLAT + j*NLON + i] = (wwind_in[k*NLON*NLAT + j*NLON + i] + wwind_in[(k+1)*NLON*NLAT + j*NLON + i])/2.0;
			  dbz[k*NLON*NLAT + j*NLON + i] = pow(10.0,(dbz_in[k*NLON*NLAT + j*NLON + i]*0.1));
			  z[k*NLON*NLAT + j*NLON + i] = (ph_in[k*NLON*NLAT + j*NLON + i] + ph_in[(k+1)*NLON*NLAT + j*NLON + i]
				  + phb_in[k*NLON*NLAT + j*NLON + i] + phb_in[(k+1)*NLON*NLAT + j*NLON + i])/(2.0*9.81);
		  }
	  }
  }
   // The file is automatically closed by the destructor. This frees
   // up any internal netCDF resources associated with the file, and
   // flushes any buffers.

   std::cout << "*** SUCCESS reading netCDF!" << std::endl;
   return 0;
}

bool WRF::getData(const double &lat,const double &lon,const double &alt, double &u_out, double &v_out, double &w_out, double &dbz_out)
{
	int altup, altdown, latup, latdown, lonup, londown;
	double altwgt, latwgt, lonwgt;
	
	// Find the closest point
	float mindist = 1e34;
	float minalt = 1e34;
	
	int loni, lati, alti;
    for (int i = 0; i < NLON; i++) {
		for (int j = 0; j < NLAT; j++) {
			float dist = sqrt((lon - lons[j*NLON+i])*(lon - lons[j*NLON+i]) + (lat - lats[j*NLON+i])*(lat - lats[j*NLON+i]));
			if (dist < mindist) {
				mindist = dist;
				loni = i;
				lati = j;
			}
		}
	}
 	for (int k = 0; k < NLVL; k++) {
		if (fabs(alt - z[k*NLON*NLAT + lati*NLON + loni]) < minalt) {
			minalt = fabs(alt - z[k*NLON*NLAT + lati*NLON + loni]);
			alti = k;
		}
	}
	
	// Trilinear interpolate
	if ((mindist != 1e34) and (minalt != 1e34)) {
		double altdiff = alt - z[alti*NLON*NLAT + lati*NLON + loni];
		if (altdiff > 0) {
			if (alti != NLVL-1) {
				altup = alti+1;
				altdown = alti;
				altwgt = altdiff / (z[altup*NLON*NLAT + lati*NLON + loni]-z[altdown*NLON*NLAT + lati*NLON + loni]);
			} else {
				altup = alti;
				altdown = alti;
				altwgt = 1.0;
			}
		} else {
			if (alti != 0) {
				altup = alti;
				altdown = alti-1;
				altwgt = (alt - z[altdown*NLON*NLAT + lati*NLON + loni]) / (z[altup*NLON*NLAT + lati*NLON + loni]-z[altdown*NLON*NLAT + lati*NLON + loni]);
				
			} else {
				altup = alti;
				altdown = alti;
				altwgt = 1.0;
			}
		}
		
		double latdiff = lat - lats[lati*NLON+loni];
		if (latdiff > 0) {
			if (lati != NLAT-1) {
				latup = lati+1;
				latdown = lati;
				latwgt = latdiff / (lats[latup*NLON+loni] - lats[latdown*NLON+loni]);
			} else {
				latup = lati;
				latdown = lati;
				latwgt = 1.0;
			}
		} else {
			if (lati != 0) {
				latup = lati;
				latdown = lati-1;
				latwgt = (lat - lats[latdown*NLON+loni]) / (lats[latup*NLON+loni] - lats[latdown*NLON+loni]);
			} else {
				latup = lati;
				latdown = lati;
				latwgt = 1.0;
			}
		}
		
		double londiff = lon - lons[lati*NLON+loni];
		if (londiff > 0) {
			if (loni != NLAT-1) {
				lonup = loni+1;
				londown = loni;
				lonwgt = londiff / (lons[lati*NLON+lonup] - lons[lati*NLON+londown]);
			} else {
				lonup = loni;
				londown = loni;
				lonwgt = 1.0;
			}
		} else {
			if (loni != 0) {
				lonup = loni;
				londown = loni-1;
				lonwgt = (lon - lons[lati*NLON+londown]) / (lons[lati*NLON+lonup] - lons[lati*NLON+londown]);
			} else {
				lonup = loni;
				londown = loni;
				lonwgt = 1.0;
			}
		}
		double c00 = (1.0-altwgt)*u[altdown*NLON*NLAT + latdown*NLON + londown] + altwgt*u[altup*NLON*NLAT + latdown*NLON + londown];
		double c10 = (1.0-altwgt)*u[altdown*NLON*NLAT + latup*NLON + londown] + altwgt*u[altup*NLON*NLAT + latup*NLON + londown];
		double c01 = (1.0-altwgt)*u[altdown*NLON*NLAT + latdown*NLON + lonup] + altwgt*u[altup*NLON*NLAT + latdown*NLON + lonup];
		double c11 = (1.0-altwgt)*u[altdown*NLON*NLAT + latup*NLON + lonup] + altwgt*u[altup*NLON*NLAT + latup*NLON + lonup];
		double c0 = (1.0-latwgt)*c00 + latwgt*c10;
		double c1 = (1.0-latwgt)*c01 + latwgt*c11;
		u_out = (1.0-lonwgt)*c0 + lonwgt*c1;
		
		c00 = (1.0-altwgt)*v[altdown*NLON*NLAT + latdown*NLON + londown] + altwgt*v[altup*NLON*NLAT + latdown*NLON + londown];
		c10 = (1.0-altwgt)*v[altdown*NLON*NLAT + latup*NLON + londown] + altwgt*v[altup*NLON*NLAT + latup*NLON + londown];
		c01 = (1.0-altwgt)*v[altdown*NLON*NLAT + latdown*NLON + lonup] + altwgt*v[altup*NLON*NLAT + latdown*NLON + lonup];
		c11 = (1.0-altwgt)*v[altdown*NLON*NLAT + latup*NLON + lonup] + altwgt*v[altup*NLON*NLAT + latup*NLON + lonup];
		c0 = (1.0-latwgt)*c00 + latwgt*c10;
		c1 = (1.0-latwgt)*c01 + latwgt*c11;
		v_out = (1.0-lonwgt)*c0 + lonwgt*c1;
		
		c00 = (1.0-altwgt)*w[altdown*NLON*NLAT + latdown*NLON + londown] + altwgt*w[altup*NLON*NLAT + latdown*NLON + londown];
		c10 = (1.0-altwgt)*w[altdown*NLON*NLAT + latup*NLON + londown] + altwgt*w[altup*NLON*NLAT + latup*NLON + londown];
		c01 = (1.0-altwgt)*w[altdown*NLON*NLAT + latdown*NLON + lonup] + altwgt*w[altup*NLON*NLAT + latdown*NLON + lonup];
		c11 = (1.0-altwgt)*w[altdown*NLON*NLAT + latup*NLON + lonup] + altwgt*w[altup*NLON*NLAT + latup*NLON + lonup];
		c0 = (1.0-latwgt)*c00 + latwgt*c10;
		c1 = (1.0-latwgt)*c01 + latwgt*c11;
		w_out = (1.0-lonwgt)*c0 + lonwgt*c1;
		
		c00 = (1.0-altwgt)*dbz[altdown*NLON*NLAT + latdown*NLON + londown] + altwgt*dbz[altup*NLON*NLAT + latdown*NLON + londown];
		c10 = (1.0-altwgt)*dbz[altdown*NLON*NLAT + latup*NLON + londown] + altwgt*dbz[altup*NLON*NLAT + latup*NLON + londown];
		c01 = (1.0-altwgt)*dbz[altdown*NLON*NLAT + latdown*NLON + lonup] + altwgt*dbz[altup*NLON*NLAT + latdown*NLON + lonup];
		c11 = (1.0-altwgt)*dbz[altdown*NLON*NLAT + latup*NLON + lonup] + altwgt*dbz[altup*NLON*NLAT + latup*NLON + lonup];
		c0 = (1.0-latwgt)*c00 + latwgt*c10;
		c1 = (1.0-latwgt)*c01 + latwgt*c11;
		dbz_out = (1.0-lonwgt)*c0 + lonwgt*c1;
						
		return true;
	}
	return false;

}

