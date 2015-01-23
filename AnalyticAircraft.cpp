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
#include "ReferenceState.h"

AnalyticAircraft::AnalyticAircraft(const QString& in, const QString& out, const QString& suffix, const QDomElement& config)
{

	// Setup the data path
	dataPath = QDir(in);
	outPath = QDir(out);
	swpSuffix = suffix;
    parseXMLconfig(config);
	readSwpDir();

	//analyticType = analytic;
	Pi = acos(-1);
	beamwidth= -999.;
	msecOffset = 0;
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
	double refLat = configHash.value("ref_lat").toFloat();
	double refLon = configHash.value("ref_lon").toFloat();
	QTime refTime(configHash.value("ref_hr").toFloat(),
                  configHash.value("ref_min").toFloat(),
                  configHash.value("ref_sec").toFloat());
	QDate refDate(configHash.value("ref_year").toFloat(),
                  configHash.value("ref_mon").toFloat(),
                  configHash.value("ref_day").toFloat());
	QDateTime refDateTime(refDate, refTime, Qt::UTC);
	QString demfile = configHash.value("dem_file");

	// Load the DEM
	if(!asterDEM.readDem(demfile.toAscii().data())) return false;

	QString mode = configHash.value("analytic");
	int analytic = 0;
	if (mode == "beltrami") {
		analytic = beltrami;
	} else if (mode == "wrf") {
		analytic = wrf;
		wrfFile.readWRF(configHash.value("wrf_file").toAscii().data());
	} else if (mode == "constant") {
		analytic = constant;
	} else if (mode == "cylindrical") {
		analytic = cylindrical;
	}

	// Resample an analytic field
	if (getfileListsize()) {
		for (int f = 0; f < getfileListsize(); ++f) {
			if (load(f)) {
				printf("\n\nProcessing file %d\n", f);
				// Create an analytic track, then add some navigation errors
				clearCfacs();
				if (configHash.value("track") == "analytic") {
					if (f == 0) {
						ryib_info* ryptr = swpfile.getRyibBlock(0);
						QTime rayTime(ryptr->hour, ryptr->min, ryptr->sec, ryptr->msec);
						QDateTime rayDateTime(refDateTime.date(), rayTime);
						msecOffset = refDateTime.time().msecsTo(rayTime);
					}
					analyticTrack(refLat, refLon, refDateTime, analytic);
				}
				recalculateAirborneAngles();
				resample_wind(refLat, refLon, analytic);
				addNavError(refLat, refLon, refDateTime, analytic);
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

void AnalyticAircraft::analyticTrack(double refLat, double refLon, QDateTime refDateTime, int analytic)
{

	GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();

	double ns_gspeed = configHash.value("ns_gspeed").toFloat();
	double ew_gspeed = configHash.value("ew_gspeed").toFloat();
	double refX, refY;
	tm.Forward(refLon, refLat, refLon, refX, refY);

	QFile insituFile("insitudata.txt");
	if (!insituFile.open(QIODevice::Append | QIODevice::Text))
		std::cerr << "Problem opening in situ file\n";
	QTextStream outstream(&insituFile);
	outstream.setRealNumberPrecision(4);
	outstream.setFieldWidth(10);
    outstream.setRealNumberNotation(QTextStream::FixedNotation);
	vold_info* vptr = swpfile.getVolumeBlock();
	vptr->year = refDateTime.date().year();
	vptr->mon = refDateTime.date().month();
	vptr->day = refDateTime.date().day();
	vptr->hour = refDateTime.time().hour();
	vptr->min = refDateTime.time().minute();
	vptr->sec = refDateTime.time().second();

	for (int i=0; i < swpfile.getNumRays(); i++) {
		asib_info* aptr = swpfile.getAircraftBlock(i);
		ryib_info* ryptr = swpfile.getRyibBlock(i);
		ryptr->julian_day = refDateTime.date().dayOfYear();
		QTime rayTime(ryptr->hour, ryptr->min, ryptr->sec, ryptr->msec);
		QDateTime rayDateTime(refDateTime.date(), rayTime);
		int msecElapsed = refDateTime.time().msecsTo(rayTime);
		msecElapsed -= msecOffset;
		rayTime = refDateTime.time().addMSecs(msecElapsed);
		ryptr->hour = rayTime.hour();
		ryptr->min = rayTime.minute();
		ryptr->sec = rayTime.second();
		ryptr->msec = rayTime.msec();
		double radarLat, radarLon,radarAlt;
		double radarX = refX + ew_gspeed * msecElapsed/1000.0;
		double radarY = refY + ns_gspeed * msecElapsed/1000.0;
		radarAlt = configHash.value("radar_alt").toFloat() / 1000.0;

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
		double hwavelength = configHash.value("hwavelength").toFloat();
		double vwavelength = configHash.value("vwavelength").toFloat();
		if (analytic == beltrami) {
			BeltramiFlow(hwavelength, vwavelength, x, y, z, t, h, u, v, w, dz);
		} else if (analytic == wrf) {
			WrfResample(radarLat, radarLon, z, t, h, u, v, w, dz);
		} else if (analytic == constant) {
			ConstantWind(x, y, z, t, h, u, v, w, dz);
		} else if (analytic == cylindrical) {
			CylindricalWind(x, y, z, t, h, u, v, w, dz);
		}
		aptr->lon = radarLon;
		aptr->lat = radarLat;
		aptr->alt_msl= radarAlt;
		aptr->alt_agl= radarAlt;
		aptr->ew_gspeed= ew_gspeed;
		aptr->ns_gspeed = ns_gspeed;
		aptr->vert_vel= 0.;
		aptr->head= configHash.value("heading").toFloat();
		aptr->roll= 0.;
		aptr->pitch= 0.;
		aptr->drift= 0.;
		aptr->ew_horiz_wind= u;
		aptr->ns_horiz_wind= v;
		aptr->vert_wind= w;
		aptr->head_change= 0.;
		aptr->pitch_change= 0.;
		aptr->tilt_ang=configHash.value("tilt_angle").toFloat();
		outstream << ryptr->hour << ryptr->min << ryptr->sec << ryptr->msec << msecElapsed << radarLat << radarLon << radarAlt << u << v << w << endl;
	}
	insituFile.close();

}

void AnalyticAircraft::clearCfacs()
{

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
}


void AnalyticAircraft::addNavError(double refLat, double refLon, QDateTime refDateTime, int analytic)
{

	for (int i=0; i < swpfile.getNumRays(); i++) {
		asib_info* aptr = swpfile.getAircraftBlock(i);
		//ryib_info* ryptr = swpfile.getRyibBlock(i);

		aptr->lon += configHash.value("lon_error").toFloat();
		aptr->lat += configHash.value("lat_error").toFloat();
		aptr->alt_msl += configHash.value("alt_error").toFloat();
		aptr->alt_agl += configHash.value("alt_error").toFloat();
		aptr->ew_gspeed += configHash.value("ew_error").toFloat();
		aptr->ns_gspeed += configHash.value("ns_error").toFloat();
		aptr->vert_vel += configHash.value("vv_error").toFloat();
		aptr->head += configHash.value("heading_error").toFloat();
		aptr->roll += configHash.value("roll_error").toFloat();
		aptr->pitch += configHash.value("pitch_error").toFloat();
		aptr->drift += configHash.value("drift_error").toFloat();
		aptr->tilt_ang += configHash.value("tilt_error").toFloat();
		aptr->head_change= 0.;
		aptr->pitch_change= 0.;
	}

}


void AnalyticAircraft::resample_wind(double refLat, double refLon, int analytic)
{
	ReferenceState* refstate = new ReferenceState("dunion_mt.snd");
	// Resample the Doppler field with a Beltrami flow from Shapiro et al. 2009
#pragma omp parallel for
	for (int i=0; i < swpfile.getNumRays(); i++) {
		GeographicLib::TransverseMercatorExact tm = GeographicLib::TransverseMercatorExact::UTM();
		double refX, refY;
		tm.Forward(refLon, refLat, refLon, refX, refY);
		int maxElevation = asterDEM.getMaxElevation();
		float nyquist = swpfile.getNyquistVelocity();

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
			float range = gatesp[n]+configHash.value("range_delay_error").toFloat();
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
				tm.Reverse(refLon, radarX + relX, radarY + relY, absLat, absLon);
				if ((el <= 0) and (z < maxElevation)) {
					h = asterDEM.getElevation(absLat, absLon);
					if (h < 0) h=0;
				} else {
					h = 0;
				}
				if (analytic == beltrami) {
					double vwavelength = configHash.value("vwavelength").toFloat();
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
					WrfResample(absLat, absLon, z, t, h, u, v, w, dz);
				} else if (analytic == constant) {
					ConstantWind(x, y, z, t, h, u, v, w, dz);
				} else if (analytic == cylindrical) {
					CylindricalWind(x, y, z, t, h, u, v, w, dz);
				}
				// Background reflectivity noise floor
				//if (dz == 0) dz = pow(10.0,((-89.339 + 19.295 * log10(range))/10.0));

				/* The default beamwidth is set to -999 which assumes the beam is infinitely small.
				    Increasing it to realistic values and/or changing the beam pattern to include sidelobes increases
				    the calculation time significantly but gives a more realistic representation of the winds */
				beamwidth = configHash.value("beamwidth").toFloat();

				if (beamwidth < 0) {
                    /* if ((z-h) < 75) {
                        refdata[n] = 50.;
                        swdata[n] = 0.;
                        ncpdata[n] = 1.;
                        veldata[n] = 0.;
                        velcorr[n] = 0.;
                    } else { */
                        refdata[n] = 10*log10(dz);
                        swdata[n] = 0.;
                        ncpdata[n] = 1.;

	                    // Fall speed
	                    double Z = refdata[n];
	                    double H = z;
	                    double ZZ=pow(10.0,(Z*0.1));
	                    double melting_zone = 1000.0;
	                    double hlow= 5000.0;
	                    double hhi= hlow + melting_zone;

	                    /* density correction term (rhoo/rho)*0.45
	                     0.45 density correction from Beard (1985, JOAT pp 468-471) */
	                    double rho = refstate->getReferenceVariable(ReferenceVariable::rhoref, H);
	                    double rhosfc = refstate->getReferenceVariable(ReferenceVariable::rhoref, 0.);
	                    double DCOR = pow((rhosfc/rho),(double)0.45);

	                    // The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s)
	                    double VTS=-DCOR * (0.817*pow(ZZ,(double)0.063));

	                    // The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s) */
	                    double VTR=-DCOR * (2.6*pow(ZZ,(double).107));

	                    /* Test if height is in the transition region between SNOW and RAIN
	                     defined as hlow in km < H < hhi in km
	                     if in the transition region do a linear weight of VTR and VTS */
	                    double mixed_dbz = 20.0;
	                    double rain_dbz = 30.0;
	                    if ((Z > mixed_dbz) and
	                        (Z <= rain_dbz)) {
	                        double WEIGHTR=(Z-mixed_dbz)/(rain_dbz - mixed_dbz);
	                        double WEIGHTS=1.-WEIGHTR;
	                        VTS=(VTR*WEIGHTR+VTS*WEIGHTS)/(WEIGHTR+WEIGHTS);
	                    } else if (Z > rain_dbz) {
	                        VTS=VTR;
	                    }
	                    double w_term=VTR*(hhi-H)/melting_zone + VTS*(H-hlow)/melting_zone;
	                    if (H < hlow) w_term=VTR;
	                    if (H > hhi) w_term=VTS;

                        double vr = u*sin(az)*cos(el) + v*cos(az)*cos(el) + (w+w_term)*sin(el);
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
                        double noise = configHash.value("noise").toFloat();
                        if (noise > 0) {
														double u1, u2, g;
														g = 2.0;
														while ( g >= 1 ) {
															u1 = 2 * float(rand() % 1000)/1000.0 - 1;
															u2 = 2 * float(rand() % 1000)/1000.0 - 1;
															g = u1*u1 + u2*u2;
															//std::cout << g << "\t" << u1 << "\t" << u2 << "\n";
														}

														g = sqrt( (-2 * log(g))  / g );
														double g1 = u2 * g;
                            //noise = (rand() % int(10*noise) + 1) * ((rand() % 3) - 1) / 10.0;
                            vr += g1*noise;
                        }
                        //double dznoise = rand() % ((int)refdata[n] + 20) + 1;
                        //double noise = configHash.value("noise").toFloat() * 1/(dznoise) * ((rand() % 2) - 1);
                        veldata[n] = vr;
                    //}
				} else {
					// Loop over the width of the beam
					double maxbeam = Pi;
					double beamincr = maxbeam/90.;
                    int beamtype = configHash.value("beamtype").toFloat();

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
						double beamaxis = r; //(beamwidth*Pi/180.);

						for (double theta=0; theta < 360; theta += 5) {
							double azmod = az + r*cos(theta * Pi / 180.);
							double elmod = el + r*sin(theta * Pi / 180.);
							relX = range*sin(azmod)*cos(elmod);
							relY = range*cos(azmod)*cos(elmod);
							relZ = sqrt(range*range + rEarth*rEarth + 2.0 * range * rEarth * sin(elmod)) - rEarth;
                            double power;
                            if (beamtype == 0) {
                                // Gaussian beam
                                power = exp(-(beamaxis*beamaxis)/1.443695);
							} else if (beamtype == 1) {
                                // Rectangular beam (ELDORA-like)
                                power = fabs(sin(27*sin(beamaxis))/(27*sin(beamaxis)));
                                if (beamaxis > Pi/2.) beamaxis = Pi/2.; // power = power * 0.001;
                                power = pow(10.0,4.5*log10(power));
                                //if (beamaxis > Pi/2.) power = 0.000000316227766;
							} else if (beamtype == 2) {
                                // Gnarly sidelobes
                                power = (sin(10*beamaxis)/sin(beamaxis));
                                power = (power > -0.1) ? fabs(power) : 0.1;
                            }
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
                            double hwavelength = configHash.value("hwavelength").toFloat();
                            double vwavelength = configHash.value("vwavelength").toFloat();
							if (analytic == beltrami) {
								BeltramiFlow(hwavelength, vwavelength, x, y, z, t, h, u, v, w, dz);
							} else if (analytic == wrf) {
								WrfResample(x, y, z, t, h, u, v, w, dz);
							} else if (analytic == constant) {
								ConstantWind(x, y, z, t, h, u, v, w, dz);
							} else if (analytic == cylindrical) {
								CylindricalWind(x, y, z, t, h, u, v, w, dz);
							}
							if (dz == 100000.0) dz = 100000.0 * sin(elmod)*sin(elmod);

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
							//double dznoise = 10*log10(dz) + 40.;
							//if (dznoise > 1.) dznoise = 1.;
							double noise = configHash.value("noise").toFloat();
							if (noise > 0) {
								double u1, u2, g;
								g = 2.0;
								while ( g >= 1 ) {
									u1 = 2 * float(rand() % 1000)/1000.0 - 1;
									u2 = 2 * float(rand() % 1000)/1000.0 - 1;
									g = u1*u1 + u2*u2;
									//std::cout << g << "\t" << u1 << "\t" << u2 << "\n";
								}

								g = sqrt( (-2 * log(g))  / g );
								double g1 = u2 * g;
								//noise = (rand() % int(10*noise) + 1) * ((rand() % 3) - 1) / 10.0;
								vr += g1*noise;
							}
							veltmp += vr*dz*power;
							weight += dz*power;
							double vrmean = veltmp/weight;
							swtmp += dz*power*(vr - vrmean)*(vr - vrmean);
							//ncptmp += 1/swtmp;
						}
					}
					double bgdz = pow(10.0,((-89.339 + 19.295 * log10(range))/10.0));
					if (reftmp < bgdz) reftmp = bgdz;
					refdata[n] = 10*log10(reftmp);
					swdata[n] = sqrt(swtmp/weight);
					ncpdata[n] = (1 - swdata[n] / 8.) + refdata[n] / 50.;
					if (ncpdata[n] < 0.0) ncpdata[n] = 0.01;
					if (ncpdata[n] > 1.0) ncpdata[n] = 1.0;
					veldata[n] = veltmp/weight;
					// Add some flecks of bad data
					if (ncpdata[n] < 0.4) {
						if ((rand() % 100) < 10) veldata[n] += configHash.value("noise").toFloat() *
                            (rand() % 50) *((rand() % 2) - 1);
					}
					if (ncpdata[n] < 0.2) veldata[n] += configHash.value("noise").toFloat() *
                        (rand() % 50) *((rand() % 2) - 1);

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
	delete refstate;
}

void AnalyticAircraft::BeltramiFlow(double hwavelength, double vwavelength, double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz)
{
	double k = 2*Pi / hwavelength; // Horizontal wavelengths
	double l = k;
	double m = 2*Pi / vwavelength;
	double A = configHash.value("peak_w").toFloat(); // Peak Vertical velocity
	double amp = A / (k*k + l*l);
	double wavenum = sqrt(k*k + l*l + m*m);
	double U = configHash.value("mean_u").toFloat();
	double V = configHash.value("mean_v").toFloat();
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
	//if (fabs(z-h) < 75.0) {
    if (z-h < 75.0) {
		u = v = w = 0.0;
		dz = 100000.0;
	}
}

void AnalyticAircraft::WrfResample(double lat, double lon, double z, double t, double h, double &u, double &v, double &w, double &dz)
{
    if ((z-h) < 0.0) {
		u = v = w = 0.0;
		dz = 100000.0;
	} else {
		double dbz = 0.0;
		if (!wrfFile.getData(lat,lon,z,u,v,w,dbz)) {
			u = v = w = dbz = 0.0;
		}
		dz = dbz;
		//dz = pow(10.0,(dbz*0.1));
	}

//height(i,k,j) = 0.5*(phb(i,k,j)+phb(i,k+1,j)+ph(i,k,j)+ph(i,k+1,j))/9.81

}

void AnalyticAircraft::ConstantWind(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz)
{
	u = configHash.value("mean_u").toFloat();
	v = configHash.value("mean_v").toFloat();
	w = 0.0;
	dz = 0.01;
    if (fabs(z-h) < 75.0) {
		u = v = w = 0.0;
		dz = 100000.0;
	}
}

void AnalyticAircraft::CylindricalWind(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz)
{
	double vr = configHash.value("mean_u").toFloat();
	double vt = configHash.value("mean_v").toFloat();
	double az = atan2(y,x);
	u = vr*cos(az) - vt*sin(az);
	v = vr*sin(az) + vt*cos(az);
	w = configHash.value("mean_w").toFloat();
	dz = 0.01;
    if (fabs(z-h) < 75.0) {
		u = v = w = 0.0;
		dz = 100000.0;
	}
}

bool AnalyticAircraft::parseXMLconfig(const QDomElement& config)
{

    std::cout << "Parsing configuration file...\n";

	// Parse the nodes to a hash
	QDomNodeList nodeList = config.childNodes();
	for (int i = 0; i < nodeList.count(); i++) {
		QDomNode currNode = nodeList.item(i);
		// Check to see if this is a set of pass parameters
		QString iter = 0;
		if (currNode.hasAttributes() and currNode.attributes().contains("iter")) {
			iter = currNode.toElement().attribute("iter");
		}
		QDomNodeList configList = currNode.childNodes();
		for (int j = 0; j < configList.count(); j++) {
			QDomNode configItem = configList.item(j);
			QDomElement group = configItem.toElement();
			QString tag = group.tagName();
			if (iter.toInt() > 1) {
				// Append the pass number to the tagName
				tag += "_" + iter;
			}
			if (!group.text().isEmpty()) {
				configHash.insert(tag, group.text());
                std::cout << tag.toStdString() << " => " << configHash.value(tag).toStdString() << std::endl;
			}
		}
	}

	// Validate the hash -- multiple passes are not validated currently
	QStringList configKeys;
	configKeys << "ref_lat" << "ref_lon" << "ref_hr" << "ref_min" << "ref_sec"
		<< "dem_file" << "analytic" << "ns_gspeed" << "radar_alt" << "beamtype"
		<< "beamwidth" << "tilt_angle" << "lon_error" << "lat_error" << "alt_error"
		<< "ew_error" << "ns_error" << "vv_error" << "heading_error" << "roll_error"
		<< "pitch_error" << "drift_error" << "tilt_error" << "range_delay_error"
		<< "hwavelength" << "vwavelength" << "mean_u" << "mean_v" << "peak_w" << "noise"
		<< "heading" << "track";
 	for (int i = 0; i < configKeys.count(); i++) {
		if (!configHash.contains(configKeys.at(i))) {
            std::cout <<	"No configuration found for <" << configKeys.at(i).toStdString() << "> aborting..." << std::endl;
			return false;
		}
	}
	return true;

}
