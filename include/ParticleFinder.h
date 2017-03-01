/*
 *  ParticleFinder.h
 * 
 *  A ParticleFinder object knows how to find particles given a pixel array.
 *  Currently, it does this by using the 1D Gaussian Estimator procedure 
 *  described in:
 *
 *  N.T. Ouellette, H. Xu, and E. Bodenschatz, "A quantitative study of three-
 *  dimensional Lagrangian particle tracking algorithms," Exp. Fluids 40, 301-313 (2006).
 *
 *  To change the particle-finding method, re-implement the constructor.
 *
 *  Last update: 10/12/11 by NTO
 *
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 */

#ifndef PARTICLEFINDER_H
#define PARTICLEFINDER_H

#include <string>
#include <deque>
#include <stdexcept>

#include <Frame.h>

class ParticleFinder {

public:
  // constructor: process the given pixel array
  ParticleFinder(int**& p, int rows, int cols, int depth, int threshold /*int thresh*/) throw(std::out_of_range);
  // destructor
  ~ParticleFinder() {};

  // write the particle centers to a file
  void WriteToFile(std::string filename);
  // make a Frame object with the particle positions.
  Frame CreateFrame();
   
  // return the number of particles found
  int NumParticles() const;
  
  // remove large clusters
  void Squash(double rad);

private:
	int** pixels;
	
  // store vectors of the x and y coordinates of the particles
  std::deque<double> x;
  std::deque<double> y;

	// helper function
  bool IsLocalMax(int r, int c);

};

inline int ParticleFinder::NumParticles() const
{
  return x.size();
}

#endif // PARTICLEFINDER_H
