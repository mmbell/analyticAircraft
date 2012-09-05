/*
 *  AnalyticAircraft.cpp
 *
 *  Created by Michael Bell
 *  Copyright 2010. All rights reserved.
 *
 */

#include "AnalyticAircraft.h"
#include <iterator>
#include <fstream>
#include <iostream>
#include <QFile>
#include <QTextStream>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include "read_dorade.h"

AnalyticAircraft::AnalyticAircraft(const QString& in, const QString& out, const QString& suffix, const int& analytic)
{
	
	// Setup the data path
	dataPath = QDir(in);
	outPath = QDir(out);
	swpSuffix = suffix;
	readSwpDir();
	
	analyticType = analytic;
	Pi = acos(-1);
	beamwidth= -999.;
	
}

AnalyticAircraft::~AnalyticAircraft()
{
}

bool AnalyticAircraft::readSwpDir()
{
	dataPath.setNameFilters(QStringList("swp.*"));
	dataPath.setFilter(QDir::Files);
	dataPath.setSorting(QDir::Name);
	QStringList filenames = dataPath.entryList();
	
	// Read in the list sweepfiles
	for (int i = 0; i < filenames.size(); ++i) {
		QStringList fileparts = filenames.at(i).split(".");
		if (fileparts.size() == 6) swpfileList.append(filenames.at(i));
	}
	
	if (swpfileList.isEmpty()) { 
		return false; 
	} else {
		return true;
	}
	
}

bool AnalyticAircraft::processSweeps()
{

	// Set these parameters!
	double refLat = 16.5;
	double refLon = 148.0;
	QTime refTime(23,48);
	QString demfile = "ASTGTM_N16E146_dem.tif";
	
	// Load the DEM
	if(!asterDEM.readDem(demfile.toAscii().data())) return false;
	
	// Resample an analytic field
	if (getfileListsize()) {
		for (int f = 0; f < getfileListsize(); ++f) {
			if (load(f)) {
				printf("\n\nProcessing file %d\n", f);
				// Create an analytic track, then add some navigation errors
				analyticTrack(refLat, refLon, refTime, beltrami);
				recalculateAirborneAngles();
				resample_wind(refLat, refLon, beltrami);
				addNavError(refLat, refLon, refTime, beltrami);
				saveQCedSwp(f);
			} else {
				printf("\n\nError loading file %d\n", f);
				return false;
			}
		} 	
	} else {
		std::cout << "No swp files exist in " << dataPath.dirName().toStdString() << "\n"; 
		return false;
	}
		
	return true;
}

/****************************************************************************************
 ** load : This function loads the information from an individual swp file given the index.
 ****************************************************************************************/
bool AnalyticAircraft::load(const int& swpIndex)
{
	QString filename = dataPath.absolutePath() + "/" + getswpfileName(swpIndex);
	swpfile.setFilename(filename);
	// Read in the swp file
	if(swpfile.readSwpfile()) 
		return true;
	
	return false;
	
}

/****************************************************************************************
 ** save : This function saves a modified swp file with suffix appended to the end.
 ****************************************************************************************/
bool AnalyticAircraft::saveQCedSwp(const int& swpIndex)
{
	QString qcfilename = outPath.absolutePath() + "/" + getswpfileName(swpIndex) + "." + swpSuffix;
	if (swpfile.writeSwpfile(qcfilename)) return true;
	
	return false;
}

void AnalyticAircraft::recalculateAirborneAngles()
{
	swpfile.recalculateAirborneAngles();
}	

void AnalyticAircraft::analyticTrack(double refLat, double refLon, QTime refTime, int analytic)
{
	
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	
	// Clear the cfac block or manually enter values here
	cfac_info* cfptr = swpfile.getCfacBlock();
	cfptr->c_azimuth = 0.0;
	cfptr->c_elevation = 0.0;
	cfptr->c_range_delay = 0.0;
	cfptr->c_rad_lon = 0.0;
	cfptr->c_rad_lat = 0.0;
	cfptr->c_alt_msl = 0.0;
	cfptr->c_alt_agl = 0.0;
	cfptr->c_ew_grspeed = 0.0;
	cfptr->c_ns_grspeed = 0.0;
	cfptr->c_vert_vel = 0.0;
	cfptr->c_head = 0.0;
	cfptr->c_roll = 0.0;
	cfptr->c_pitch = 0.0;
	cfptr->c_drift = 0.0;
	cfptr->c_rotang = 0.0;
	cfptr->c_tiltang = 0.0;

	double ns_gspeed = 120.0;
	double ew_gspeed = 0.0;
	double refX, refY;
	tm.Forward(refLon, refLat, refLon, refX, refY);
	for (int i=0; i < swpfile.getNumRays(); i++) {		
		asib_info* aptr = swpfile.getAircraftBlock(i);
		ryib_info* ryptr = swpfile.getRyibBlock(i);

		QTime rayTime(ryptr->hour, ryptr->min, ryptr->sec, ryptr->msec);
		int msecElapsed = refTime.msecsTo(rayTime);

		double radarLat, radarLon,radarAlt;
		double radarX = refX + ew_gspeed * msecElapsed/1000.0;
		double radarY = refY + ns_gspeed * msecElapsed/1000.0;
		radarAlt = 3.0;
		
		// If the aircraft is below the ground there is a problem
		tm.Reverse(refLon, radarX, radarY, radarLat, radarLon);
		int h = asterDEM.getElevation(radarLat, radarLon);
		if (radarAlt*1000 < h) {
			std::cout << "Problem with heights! Aircraft below ground\n";
		}

		double x = radarX - refX;
		double y = radarY - refY;
		double z = radarAlt*1000;
		double t = 0.;
		double u, v, w, dz;
		double hwavelength = 10000.;
		double vwavelength = 32000.;
		if (analytic == beltrami) {
			BeltramiFlow(hwavelength, vwavelength, x, y, z, t, h, u, v, w, dz);
		} else if (analytic == wrf) {
			WrfResample(x, y, z, t, h, u, v, w, dz);
		}
		aptr->lon = radarLon;
		aptr->lat = radarLat;
		aptr->alt_msl= radarAlt;
		aptr->alt_agl= radarAlt;
		aptr->ew_gspeed= ew_gspeed;
		aptr->ns_gspeed = ns_gspeed;
		aptr->vert_vel= 0.;
		aptr->head= 0.;
		aptr->roll= 0.;
		aptr->pitch= 0.;
		aptr->drift= 0.;
		aptr->ew_horiz_wind= u;
		aptr->ns_horiz_wind= v;
		aptr->vert_wind= w;
		aptr->head_change= 0.;
		aptr->pitch_change= 0.;
		aptr->tilt_ang=15.6;
	}
	
}

void AnalyticAircraft::addNavError(double refLat, double refLon, QTime refTime, int analytic)
{
	
	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
		
	double ns_gspeed = 120.0;
	double ew_gspeed = 0.0;
	double refX, refY;
	tm.Forward(refLon, refLat, refLon, refX, refY);
	for (int i=0; i < swpfile.getNumRays(); i++) {		
		asib_info* aptr = swpfile.getAircraftBlock(i);
		ryib_info* ryptr = swpfile.getRyibBlock(i);
		
		QTime rayTime(ryptr->hour, ryptr->min, ryptr->sec, ryptr->msec);
		int msecElapsed = refTime.msecsTo(rayTime);
		
		double radarLat, radarLon,radarAlt;
		double radarX = ew_gspeed * msecElapsed/1000.0;
		double radarY = ns_gspeed * msecElapsed/1000.0;
		radarAlt = 3.0;
		
		// If the aircraft is below the ground there is a problem
		tm.Reverse(refLon, refX + radarX, refY + radarY, radarLat, radarLon);
		int h = asterDEM.getElevation(radarLat, radarLon);
		if (radarAlt*1000 < h) {
			std::cout << "Problem with heights! Aircraft below ground\n";
		}
		
		aptr->lon = radarLon;
		aptr->lat = radarLat;
		aptr->alt_msl= radarAlt+0.2;
		aptr->alt_agl= radarAlt+0.2;
		aptr->ew_gspeed= ew_gspeed;
		aptr->ns_gspeed = ns_gspeed+1.0;
		aptr->vert_vel= 0.;
		aptr->head= 0.;
		aptr->roll= 1.0;
		aptr->pitch= 1.5;
		aptr->drift= 0.2;
		//aptr->tilt_ang -= 0.3;
		aptr->head_change= 0.;
		aptr->pitch_change= 0.;
	}
	
}


void AnalyticAircraft::resample_wind(double refLat, double refLon, int analytic)
{

	// Resample the Doppler field with a Beltrami flow from Shapiro et al. 2009

	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM;
	double refX, refY;
	tm.Forward(refLon, refLat, refLon, refX, refY);
	int maxElevation = asterDEM.getMaxElevation();
	float nyquist = swpfile.getNyquistVelocity();
	
	for (int i=0; i < swpfile.getNumRays(); i++) {
		float az = swpfile.getAzimuth(i)*Pi/180.;
		float el = swpfile.getElevation(i)*Pi/180.;
		float radarLat = swpfile.getRadarLat(i);
		float radarLon = swpfile.getRadarLon(i);
		float radarAlt = swpfile.getRadarAlt(i);		
		float* refdata = swpfile.getRayData(i, "ZZ");
		float* veldata = swpfile.getRayData(i, "VR");
		float* velcorr = swpfile.getRayData(i, "VV");
		float* swdata = swpfile.getRayData(i, "SW");
		float* ncpdata = swpfile.getRayData(i, "NCP");
		QDateTime rayTime = swpfile.getRayTime(i);
		float* gatesp = swpfile.getGateSpacing();
		double radarX, radarY;
		tm.Forward(refLon, radarLat, radarLon, radarX, radarY);

		for (int n=0; n < swpfile.getNumGates(); n++) {
			float range = gatesp[n];
			double relX = range*sin(az)*cos(el);
			double relY = range*cos(az)*cos(el);
			double rEarth = 6371000;
			double relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(el)) - rEarth;
			
			double x = radarX + relX - refX;
			double y = radarY + relY - refY;
			double z = relZ + radarAlt*1000;
			double t = 0.;
						
			double u, v, w, dz;
			if ((z > -5000.0) and (z <= 20000.0)) {
				
				// Check the altitude if the beam is pointing downward
				double absLat, absLon, h;
				if ((el <= 0) and (z < maxElevation)) {
					tm.Reverse(refLon, radarX + relX, radarY + relY, absLat, absLon);
					h = asterDEM.getElevation(absLat, absLon);
					if (h < 0) h=0;
				} else {
					h = 0;
				}
				if (analytic == beltrami) {
					double vwavelength = 32000.;
					double utmp, vtmp, wtmp;
					u = v = w = 0.0;
					for (int wl = 16; wl < 17; wl = wl*2) {
						double hwavelength = wl*1000.;
						BeltramiFlow(hwavelength, vwavelength, x, y, z, t, h, utmp, vtmp, wtmp, dz);
						u+=utmp;
						v+=vtmp;
						w+=wtmp;
					}
					/* double k = 2*Pi/(4000.);
					u += 5*(sin(k*x)*cos(k*y));
					v += -5*(cos(k*x)*sin(k*y)); */
					
				} else if (analytic == wrf) {
					WrfResample(x, y, z, t, h, u, v, w, dz);
				}
				
				/* The default beamwidth is set to -999 which assumes the beam is infinitely small.
				    Increasing it to realistic values and/or changing the beam pattern to include sidelobes increases
				    the calculation time significantly but gives a more realistic representation of the winds */
				beamwidth = 1.8;
				
				if (beamwidth < 0) {
					refdata[n] = 10*log10(dz);
					swdata[n] = 1.;
					ncpdata[n] = 1.;
					double vr = u*sin(az)*cos(el) + v*cos(az)*cos(el) + w*sin(el);
					velcorr[n] = vr;
					
					// Add in the aircraft motion
					double aircraft_vr = swpfile.getAircraftVelocity(i);
					vr -= aircraft_vr;
					if (fabs(vr) > nyquist) {
						// Fold data back into Nyquist range
						while (vr > nyquist) {
							vr -= 2*nyquist;
						}
						while (vr < -nyquist) {
							vr += 2*nyquist;
						}
					}
					
					// Add some random noise
					double dznoise = rand() % ((int)refdata[n] + 20) + 1;
					double noise = 1/(dznoise) * ((rand() % 2) - 1);
					veldata[n] = vr + noise;
				} else {
					// Loop over the width of the beam
					double maxbeam = (beamwidth*3.)*Pi/180.;
					double beamincr = maxbeam/10.;
					// Circle in spherical plane to radar beam
					double reftmp, veltmp, velcorrtmp, swtmp, ncptmp, weight;
					reftmp = dz;
					double vr = u*sin(az)*cos(el) + v*cos(az)*cos(el) + w*sin(el);
					velcorrtmp = vr;
					double aircraft_vr = swpfile.getAircraftVelocity(i);
					vr -= aircraft_vr;
					if (fabs(vr) > nyquist) {
						// Fold data back into Nyquist range
						while (vr > nyquist) {
							vr -= 2*nyquist;
						}
						while (vr < -nyquist) {
							vr += 2*nyquist;
						}
					}					
					veltmp = vr;
					swtmp = 0.0;
					ncptmp = 1.0;
					weight = 1.0;
					for (double r=beamincr; r <= maxbeam; r += beamincr) {
						for (double theta=0; theta < 360; theta += 10) {
							double azmod = az + r*cos(theta * Pi / 180.);
							double elmod = el + r*sin(theta * Pi / 180.);
							relX = range*sin(azmod)*cos(elmod);
							relY = range*cos(azmod)*cos(elmod);
							relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(elmod)) - rEarth;
							double beamaxis = r/(beamwidth*Pi/180.); 
							// Gaussian beam
							//double power = exp(-(beamaxis*beamaxis)/1.443695);
							
							// Rectangular beam
							double power = (sin(2*Pi*beamaxis)/2*Pi*beamaxis)*(sin(2*Pi*beamaxis)/2*Pi*beamaxis);
							
							// Gnarly sidelobes
							//double power = (sin(10*beamaxis)/sin(beamaxis));
							//power = (power > -0.1) ? fabs(power) : 0.1;
							
							x = radarX + relX - refX;
							y = radarY + relY - refY;
							z = relZ + radarAlt*1000;
							// Check the altitude if the beam is pointing downward
							double absLat, absLon, h;
							if (el <= 0) {
								tm.Reverse(refLon, radarX + relX, radarY + relY, absLat, absLon);
								h = asterDEM.getElevation(absLat, absLon);
								if (h < 0) h=0;
							} else {
								h = 0;
							}
							double hwavelength = 16000.;
							double vwavelength = 32000.;
							if (analytic == beltrami) {
								BeltramiFlow(hwavelength, vwavelength, x, y, z, t, h, u, v, w, dz);
							} else if (analytic == wrf) {
								WrfResample(x, y, z, t, h, u, v, w, dz);
							}
							reftmp += dz*power;
							vr = u*sin(azmod)*cos(elmod) + v*cos(azmod)*cos(elmod) + w*sin(elmod);
							velcorrtmp += vr*power;
							// Add in the aircraft motion
							double aircraft_vr = swpfile.getAircraftVelocity(i);
							vr -= aircraft_vr;
							if (fabs(vr) > nyquist) {
								// Fold data back into Nyquist range
								while (vr > nyquist) {
									vr -= 2*nyquist;
								}
								while (vr < -nyquist) {
									vr += 2*nyquist;
								}
							}
							// Add some random noise
							double dznoise = 10*log10(dz) + 40.;
							if (dznoise > 1.) dznoise = 1.;
							double noise = (rand() % (int)dznoise) * ((rand() % 2) - 1);
							vr += noise;
							veltmp += vr*power;
							weight += power;	
							double vrmean = veltmp/weight;
							swtmp += power*(vr - vrmean)*(vr - vrmean);
							//ncptmp += 1/swtmp;
						}
					}
					refdata[n] = 10*log10(reftmp/weight);
					swdata[n] = sqrt(swtmp/weight);
					ncpdata[n] = (1 - swdata[n] / 8.) + refdata[n] / 50.;
					if (ncpdata[n] < 0.0) ncpdata[n] = 0.01;
					if (ncpdata[n] > 1.0) ncpdata[n] = 1.0;
					veldata[n] = veltmp/weight;
					// Add some flecks of bad data
					if (ncpdata[n] < 0.4) {
						if ((rand() % 100) < 10) veldata[n] += (rand() % 50) *((rand() % 2) - 1);
					}
					if (ncpdata[n] < 0.2) veldata[n] += (rand() % 50) *((rand() % 2) - 1);
					
					velcorr[n] = velcorrtmp/weight;
				}
			} else {
				refdata[n] = -32768.;
				swdata[n] = -32768.;
				ncpdata[n] = -32768.;
				veldata[n] = -32768.;
				velcorr[n] = -32768.;
			}
			

		}
	}
}

void AnalyticAircraft::BeltramiFlow(double hwavelength, double vwavelength, double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz)
{
	double k = 2*Pi / hwavelength; // Horizontal wavelengths
	double l = k;
	double m = 2*Pi / vwavelength;
	double A = 10.; // Peak Vertical velocity
	double amp = A / (k*k + l*l);
	double wavenum = sqrt(k*k + l*l + m*m);
	double U = 10.;
	double V = 10.;
	double nu = 15.11e-6;
	
	u = U - amp*(wavenum*l*cos(k*(x - U*t))*sin(l*(y-V*t))*sin(m*z) + 
						m*k*sin(k*(x-U*t))*cos(l*(y-V*t))*cos(m*z))
	*exp(-nu*wavenum*wavenum*t);
	v = V + amp*(wavenum*k*sin(k*(x - U*t))*cos(l*(y-V*t))*sin(m*z) - 
						m*l*cos(k*(x-U*t))*sin(l*(y-V*t))*cos(m*z))
	*exp(-nu*wavenum*wavenum*t);
	w = A*cos(k*(x-U*t))*cos(l*(y-V*t))*sin(m*z)*exp(-nu*wavenum*wavenum*t);
	//dz = 1.0;
	double dbz = 45.*cos(k*(x-U*t))*cos(l*(y-V*t))*cos(m*(z-1000)/2);
	if (dbz < -25.) 	dbz = -25.;
	dz = pow(10.0,(dbz*0.1));
	if (z < (h+1)) {
		u = v = w = 0.0;
		dz = 100000.0;
	}
}

void AnalyticAircraft::WrfResample(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz)
{

	
	
//height(i,k,j) = 0.5*(phb(i,k,j)+phb(i,k+1,j)+ph(i,k,j)+ph(i,k+1,j))/9.81

}
