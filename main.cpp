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
		cout << "Usage: eldoraqc /path/to/sweepfiles /path/to/output\n";
		exit(1);
	}

	QString inpath = argv[1];
	QString outpath = argv[2];	
	QString suffix = "AC";
	AnalyticAircraft AC(inpath, outpath, suffix, AnalyticAircraft::beltrami);
	if (AC.processSweeps()) return 0;

	return 1;
	
}
