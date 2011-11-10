/*
 *  AnalyticAircraft
 *  Software to generate analytic airborne Doppler radar data
 *
 *  Created by Michael Bell
 *  Copyright 2011. All rights reserved.
 *
 */

#include <iostream>
#include <QApplication>
#include <QFile>
#include <QDateTime>
#include <cmath>
#include "Dorade.h"
#include "AnalyticAircraft.h"

using namespace std;

int main (int argc, char *argv[]) {
	
	
	// Get the arguments	
	if (argc < 3) {
		cout << "Usage: analyticAircraft /path/to/sweepfiles /path/to/output\n";
                cout << "\t This program requires base sweepfiles, but the internal data\n";
		cout << "\t is overwritten with analytic wind and reflectivity data\n";
		exit(1);
	}

	QString inpath = argv[1];
	QString outpath = argv[2];	
	QString suffix = "AC";
	AnalyticAircraft AC(inpath, outpath, suffix, AnalyticAircraft::beltrami);
	if (AC.processSweeps()) return 0;

	return 1;
	
}
