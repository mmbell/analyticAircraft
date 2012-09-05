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
#include <QtXml>
#include "Dorade.h"
#include "AnalyticAircraft.h"

using namespace std;

int main (int argc, char *argv[]) {
	
	
	// Get the arguments	
	if (argc < 4) {
		cout << "Usage: analyticAircraft /path/to/sweepfiles /path/to/output configuration.xml\n";
                cout << "\t This program requires base sweepfiles, but the internal data\n";
		cout << "\t is overwritten with analytic wind and reflectivity data\n";
		exit(1);
	}

	QString inpath = argv[1];
	QString outpath = argv[2];
    
    // Check to make sure the last argument has the right suffix
    QString xmlfile(argv[3]);
    if (xmlfile.right(3) != "xml") {
        std::cout << xmlfile.toStdString() << " does not look like an XML file\n";
        return EXIT_FAILURE;
    }
    
    // Open the file
    QFile file(xmlfile);
    if (!file.open(QIODevice::ReadOnly)) {
        std::cout << "Error Opening Configuration File, Check Permissions on " << xmlfile.toStdString() << "\n";
        return EXIT_FAILURE;
    }
    
    // Create a DOM document with contents from the configuration file
    QDomDocument domDoc;
    QString errorStr;
    int errorLine;
    int errorColumn;
    if (!domDoc.setContent(&file, true, &errorStr, &errorLine, &errorColumn)) {
        // Exit on malformed XML
        QString errorReport = QString("XML Parse Error in "+xmlfile+" at Line %1, Column %2:\n%3")
        .arg(errorLine)
        .arg(errorColumn)
        .arg(errorStr);
        std::cout << errorReport.toStdString() << "\n";
        file.close();
        return EXIT_FAILURE;
    }
    
    // Successful file read
    file.close();
    
    // Check the root node to make sure this is really a configuration file
    QDomElement root = domDoc.documentElement();
    if (root.tagName() != "analyticaircraft") {
        std::cout << "The XML file " << xmlfile.toStdString() << " is not a configuration file\n.";
        return EXIT_FAILURE;
    }
    
	QString suffix = "AC";
	AnalyticAircraft AC(inpath, outpath, suffix, root);
	if (AC.processSweeps()) return 0;

	return 1;
	
}
