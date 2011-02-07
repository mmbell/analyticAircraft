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
	if (argc == 1) {
		cout << "Usage: eldoraqc /path/to/sweepfiles\n";
		exit(1);
	}
	// For now, just take a directory, we can add more arguments later
	QString path = argv[1];
	

	AnalyticAircraft AC(path);
	
	QString suffix = "AC";
	double refLat = 16.5;
	double refLon = 148.;
	QTime refTime(23,0);
	enum analytics {
		beltrami,
		wrf
	};
	
	// Resample an analytic field
	if (AC.getfileListsize()) {
		for (int f = 0; f < AC.getfileListsize(); ++f) {
			if (AC.load(f)) {
				printf("\n\nProcessing file %d\n", f);
				AC.analyticTrack(refLat, refLon, refTime, beltrami);
				AC.resample_wind(refLat, refLon, beltrami);
				AC.recalculateAirborneAngles();
				AC.resample_wind(refLat, refLon, beltrami);
				AC.saveQCedSwp(suffix);
			}
		}	
	} else printf("No swp files exist in %s\n\n", argv[1]); 
	
	return 0;
	
}
