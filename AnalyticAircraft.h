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

class AnalyticAircraft
{
	
public:
	AnalyticAircraft(QString path);
	~AnalyticAircraft();
	
	bool readSwpDir();
	bool load(const int& swpIndex);
	bool saveQCedSwp(const QString& suffix);
	
	int getfileListsize() { return swpfileList.size(); }
	QString getswpfileName(int n) { return swpfileList[n]; }

	void recalculateAirborneAngles();

	void analyticTrack(double refLat, double refLon, QTime refTime, int analytic);
	void resample_wind(double refLat, double refLon, int analytic);
	
private:
	QList<QString> swpfileList;
	QDir dataPath;
	
	Dorade swpfile;
	double Pi;

	void BeltramiFlow(double x, double y, double z, double t, double &u, double &v, double &w);
	void WrfResample(double x, double y, double z, double t, double &u, double &v, double &w);

};

#endif