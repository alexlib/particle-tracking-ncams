/*
 *  GDF.h
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 */

#ifndef GDF_H
#define GDF_H

#include <string>
#include <deque>
#include <stdexcept>

#include <Frame.h>

class GDF{

public:
	// constructor: process the given pixel array
	GDF(std::string filename) throw(std::out_of_range);
	// destructor
	~GDF();
	
	// read GDF files
	int readGDF2D(int frame);

    int seekGDF(int start);
	
	// fix header information
	void fixHeader(int nr, int cols);
	
	// make a Frame object with the particle positions.
	Frame CreateFrame();
	// return the number of particles found
	int NumParticles() const;;
	
private:

	std::string outname;
	std::ifstream infile;
	std::ofstream outfile;
	double filePos[3];

	int magic,tmpi, cols, rows;
	double xi, yi, fi, orii;
    
	int prevFrameNum;
	int currentFrameNum;
	int nextFrameNum;
    int startFrameNum;
	int missedFrame;
	
	int waiting_to_be_written;
	
	// store vectors of the x and y coordinates of the particles
	std::deque<double> x;
	std::deque<double> y;
	std::deque<double> ori;

};

inline int GDF::NumParticles() const
{
  return x.size();
}

inline GDF::~GDF()
{
  infile.close();
  outfile.close();
}

#endif // GDF_H
