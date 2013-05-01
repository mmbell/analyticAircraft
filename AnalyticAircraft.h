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
#include <QDomDocument>
#include <QHash>
#include "DEM.h"

class AnalyticAircraft
{
	
public:
	AnalyticAircraft(const QString& in, const QString& out, const QString& suffix, const QDomElement& config);
	~AnalyticAircraft();
	
	bool readSwpDir();
	bool load(const int& swpIndex);
	bool saveQCedSwp(const int& swpIndex);
	
	int getfileListsize() { return swpfileList.size(); }
	QString getswpfileName(int n) { return swpfileList[n]; }

	void recalculateAirborneAngles();

	void analyticTrack(double refLat, double refLon, QTime refTime, int analytic);
	void addNavError(double refLat, double refLon, QTime refTime, int analytic);
	void resample_wind(double refLat, double refLon, int analytic);

	bool processSweeps();
	bool parseXMLconfig(const QDomElement& config);
    
	enum analytics {
		beltrami,
		wrf,
		constant,
		cylind
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
    QHash<QString, QString> configHash;
	void BeltramiFlow(double hwavelength, double vwavelength, double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz);
	void WrfResample(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz);
	void ConstantWind(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz);
	void CylindricalWind(double x, double y, double z, double t, double h, double &u, double &v, double &w, double &dz);
};

#endif