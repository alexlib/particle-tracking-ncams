/*
 *  Calibration.h
 *  
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 */

#ifndef CALIBRATION_H
#define	CALIBRATION_H

#include <string>
#include <deque>
#include <stdexcept>
#include <utility>

#include <Camera.h>
#include <Frame.h>
#include <Position.h>

class Calibration {
public:
	// build a Calibration object from a file
	Calibration(std::string& fname);
	~Calibration() {};
    
    Position m_pos();
    
    int mcam;
	
    void writeGDFHeader(std::string filename);
	void fixHeader(int nr, int cols);
	// do the stereomatching
	Frame Stereomatch(const std::deque<Frame>& iframes, int framenumber) throw(std::runtime_error);
		
private:

	std::string outname;
	std::ofstream outfile;

	int ncams;
    
	std::deque<Camera> cams;
	// threshold distance between a line of sight and real particles *on each image
	// plane* (in mm)
	double mindist_2D;
	// threshold distance between nearby lines of sight (in 3D, in mm) to match a particle
	double mindist_3D;
	
	// create a 3D world position from multiple positions on image planes
	std::pair<double,Position> WorldPosition(std::deque<Position> ipos) throw(std::runtime_error);
	
};
inline Position Calibration::m_pos()
{
    return Position(mcam, mcam, mcam, mcam);
}    

#endif // CALIBRATION_H
