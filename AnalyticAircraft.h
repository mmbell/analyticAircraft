/*
 *  AnalyticAircraft.h
 *  analyticAircraft
 *
 *  Created by Michael Bell
 *  Copyright 2010. All rights reserved.
 *
 */

#ifndef ANALYTICAIRCRAFT_H
#define ANALYTICAIRCRAFT_H

#include "Dorade.h"
#include <QList>
#include <QDir>
#include "DEM.h"

class AnalyticAircraft
{
	
public:
	AnalyticAircraft(const QString& in, const QString& out, const QString& suffix, const int& analytic);
	~AnalyticAircraft();
	
	bool readSwpDir();
	bool load(const int& swpIndex);
	bool saveQCedSwp(const int& swpIndex);
	
	int getfileListsize() { return swpfileList.size(); }
	QString getswpfileName(int n) { return swpfileList[n]; }

	void recalculateAirborneAngles();

	void analyticTrack(double refLat, double refLon, QTime refTime, int analytic);
	void resample_wind(double refLat, double refLon, int analytic);

	bool processSweeps();
	
	enum analytics {
		beltrami,
		wrf
	};
	
private:
	QList<QString> swpfileList;
	
	QDir dataPath;
	QDir outPath;
	QString swpSuffix;
	Dorade swpfile;
	double Pi;
	int analyticType;
	int beamwidth;
	
	DEM asterDEM;
	
	void BeltramiFlow(double hwavelength, double vwavelength, double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz);
	void WrfResample(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz);

};

#endif