/*
 *  Trackfile.h
 *  
 *
 *  Created by Nicholas Ouellette on 7/28/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 *  Added third coordinate 5/18/09
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 */

#ifndef TRACKFILE_H
#define TRACKFILE_H

#include <string>
#include <stdexcept>
#include <fstream>
#include <deque>

class Trackfile {

public: 
  // constructor: takes a filename
  Trackfile(std::string filename) throw(std::invalid_argument);
  // destructor
  ~Trackfile();
  
  int NumTracks();
  bool eof() const;
  void SkipNextTrack();
  void GetNextTrack(std::deque<float>* x, std::deque<float>* y, std::deque<float>* z, 
										std::deque<float>* t, std::deque<float>* fake);
  
private:

  int npoints;
  int current_index;
  int ntracks;
  
  std::ifstream file;
};

inline Trackfile::~Trackfile()
{
  file.close();
}

inline bool Trackfile::eof() const
{
  return file.eof();
}

#endif // TRACKFILE_H
