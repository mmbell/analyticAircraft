/*
 *  DEM.cpp
 *  analyticAircraft
 *
 *  Created by Michael Bell on 2/4/11.
 *  Copyright 2011. All rights reserved.
 *
 */

#include "DEM.h"
#include "xtiffio.h"
#include "geo_simpletags.h"
#include "geovalues.h"
#include "cpl_serv.h"
#include <stdio.h>
#include <string.h>
#include <GeographicLib/TransverseMercatorExact.hpp>

DEM::DEM()
{
}

DEM::~DEM()
{
}

int DEM::getElevation(const double& lat, const double& lon)
{
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	double refX, refY, pointX, pointY;
	tm.Forward(refLon, refLat, refLon, refX, refY);	
	tm.Forward(refLon, lat, lon, pointX, pointY);
	int xIndex = (int)(pointX - refX)/dx;
	int yIndex = (int)(pointY - refY)/dy;
	int pixel = yIndex*xsize + xIndex;
	if ((pixel >= 0) and (pixel < npixels)) {
		return elevations[pixel];
	} 
	
	return -999;
	
}

bool DEM::readDem(char* fname) 
{

    char    *outfile = NULL;
    TIFF 	*tif=(TIFF*)0;  /* TIFF-level descriptor */
    GTIF	*gtif=(GTIF*)0; /* GeoKey-level descriptor */
    int	norm_print_flag = 0, proj4_print_flag = 0;
    int	inv_flag = 0, dec_flag = 1;
    int      st_test_flag = 0;
    unsigned int x, y, skip;
    GTIFDefn	defn;
    FILE* out;

    skip = 5;
	
    /*
     * Open the file, read the GeoTIFF information, and print some info to stdout. 
     */
    int flength = strlen(fname);
    outfile = (char *) malloc(flength);
    strncpy(outfile,fname,flength-4);
    strcat(outfile, ".asc\0");
    out = fopen(outfile, "w");
    printf("Writing to %s\n\n",outfile);
    free(outfile);
	
    tif=XTIFFOpen(fname,"r");
    if (!tif) goto failure;
    
    gtif = GTIFNew(tif);
    if (!gtif)
	{
		fprintf(stderr,"failed in GTIFNew\n");
		goto failure;
	}
	
    /* dump the GeoTIFF metadata to std out */
	
    GTIFPrint(gtif,0,0);
	
    if( GTIFGetDefn( gtif, &defn ) )
	{
		printf( "\n" );
		GTIFPrintDefn( &defn, stdout );
		
		if( proj4_print_flag )
		{
			printf( "\n" );
			printf( "PROJ.4 Definition: %s\n", GTIFGetProj4Defn(&defn));
		}
		int count, orient;
		double* data;
		TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &xsize );
		TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &ysize );
		TIFFGetField( tif, TIFFTAG_GEOPIXELSCALE, &count, &data);
		TIFFGetField( tif, TIFFTAG_ORIENTATION, &orient );
		printf("Orientation:%d\n",orient);
		
		GTIFPrintCorners( gtif, &defn, stdout, xsize, ysize, inv_flag, dec_flag );
		const char* project = "ASTERGDEM";
		const char* yymmdd = "090629";
		double originx, originy;
		GTIFImageToPCS( gtif, &originx, &originy);
		int lat = (int)originy*1000;
		int lon = (int)originx*1000;
		refLat = originx;
		refLon = originy;
		int xmin = 0;
		int ymin = 0;
		int nx = xsize/skip + 1;
		int ny = ysize/skip + 1;
		dx = 30 * skip;
		dy = 30 * skip;
		fprintf(out, "%12s%12s%7d%7d%7d%7d%7d%7d%7d%7d\n",
				project, yymmdd, lat, lon, xmin, ymin, nx, ny, dx,dy);
		npixels = xsize * ysize;
		int16* buf;
		tsample_t sample;
		int nbytes;
		nbytes = TIFFScanlineSize(tif);
		buf = (int16*) _TIFFmalloc(nbytes);
		elevations = (int16*) _TIFFmalloc(npixels * sizeof (int16));
		if (elevations != NULL) {
			for( y = 0; y < ysize; y++ )
			{
				uint32 row = ysize - 1 - y;
				TIFFReadScanline(tif, buf, row, sample);
				for( x = 0; x < xsize; x++ ) 
				{
					elevations[y*xsize + x] = buf[x];
				}
			}	    
			_TIFFfree(buf);
		}
		for( y = 0; y < ysize; y+=skip ) {
			for( x = 0; x < xsize; x+=skip ) {
				fprintf(out, "%6d", elevations[y*xsize + x]);
			}
			fprintf(out, "\n");
		}
		_TIFFfree(elevations);
		
	}
    

    GTIFFree(gtif);
    if( st_test_flag )
        ST_Destroy( (ST_TIFF *) tif );
    else
        XTIFFClose(tif);
    GTIFDeaccessCSV();
    return 0;
	
failure:
    fprintf(stderr,"failure in listgeo\n");
    if (tif) XTIFFClose(tif);
    if (gtif) GTIFFree(gtif);
    GTIFDeaccessCSV();
    return 1;
}

int DEM::GTIFReportACorner( GTIF *gtif, GTIFDefn *defn, FILE * fp_out,
							 const char * corner_name,
							 double x, double y, int inv_flag, int dec_flag )

{
    double	x_saved, y_saved;
	
    /* Try to transform the coordinate into PCS space */
    if( !GTIFImageToPCS( gtif, &x, &y ) )
        return FALSE;
    
    x_saved = x;
    y_saved = y;
	
    fprintf( fp_out, "%-13s ", corner_name );
	
    if( defn->Model == ModelTypeGeographic )
    {
		if (dec_flag) 
		{
			fprintf( fp_out, "(%.7f,", x );
			fprintf( fp_out, "%.7f)\n", y );
		} 
		else 
		{
			fprintf( fp_out, "(%s,", GTIFDecToDMS( x, "Long", 2 ) );
			fprintf( fp_out, "%s)\n", GTIFDecToDMS( y, "Lat", 2 ) );
		}
    }
    else
    {
        fprintf( fp_out, "(%12.3f,%12.3f)", x, y );
		
        if( GTIFProj4ToLatLong( defn, 1, &x, &y ) )
        {
			if (dec_flag) 
			{
                fprintf( fp_out, "  (%.7f,", x );
                fprintf( fp_out, "%.7f)", y );
			} 
			else 
			{
				fprintf( fp_out, "  (%s,", GTIFDecToDMS( x, "Long", 2 ) );
				fprintf( fp_out, "%s)", GTIFDecToDMS( y, "Lat", 2 ) );
			}
        }
		
        fprintf( fp_out, "\n" );
    }
	
    if( inv_flag && GTIFPCSToImage( gtif, &x_saved, &y_saved ) )
    {
        fprintf( fp_out, "      inverse (%11.3f,%11.3f)\n", x_saved, y_saved );
    }
    
    return TRUE;
}

void DEM::GTIFPrintCorners( GTIF *gtif, GTIFDefn *defn, FILE * fp_out,
							 int xsize, int ysize, int inv_flag, int dec_flag )

{
    printf( "\nCorner Coordinates:\n" );
    if( !GTIFReportACorner( gtif, defn, fp_out,
						   "Upper Left", 0.0, 0.0, inv_flag, dec_flag ) )
    {
        printf( " ... unable to transform points between pixel/line and PCS space\n" );
        return;
    }
	
    GTIFReportACorner( gtif, defn, fp_out, "Lower Left", 0.0, ysize, 
					  inv_flag, dec_flag );
    GTIFReportACorner( gtif, defn, fp_out, "Upper Right", xsize, 0.0,
					  inv_flag, dec_flag );
    GTIFReportACorner( gtif, defn, fp_out, "Lower Right", xsize, ysize,
					  inv_flag, dec_flag );
    GTIFReportACorner( gtif, defn, fp_out, "Center", xsize/2.0, ysize/2.0,
					  inv_flag, dec_flag );
}
