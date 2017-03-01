/*
 *  Trackfile.cpp
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

#include <cstring>

#include <Trackfile.h>

using namespace std;

Trackfile::Trackfile(string filename) throw(invalid_argument)
{
  ntracks = -1;
	
  file.open(filename.c_str(), ios::in | ios::binary);
  
  int magic = 82991;
  int tmpi;
  file.read(reinterpret_cast<char*>(&tmpi), 4);
  if (tmpi != magic) {
    throw invalid_argument("Not a GDF file!");
  }
  
  int ndim;
  file.read(reinterpret_cast<char*>(&ndim), 4);
  
  int ncol;
  file.read(reinterpret_cast<char*>(&ncol), 4);
  
  file.read(reinterpret_cast<char*>(&npoints), 4);
  
  int code;
  file.read(reinterpret_cast<char*>(&code), 4);
  
  int total;
  file.read(reinterpret_cast<char*>(&total), 4);
  
  // read the track index of the first track
  float cur;
  file.read(reinterpret_cast<char*>(&cur), 4);
  current_index = static_cast<int>(cur);
  // reset the get pointer
  file.seekg(-4, ios::cur);
}

int Trackfile::NumTracks()
{
  // have we already computed this?
  if (ntracks != -1) {
    return ntracks;
  }
  
  // no, so compute it now
  ntracks = 1;
  streampos current_position = file.tellg();
  int old_current_index = current_index;
  
  // loop over the whole file
  while (!file.eof()) {
    float index;
		file.read(reinterpret_cast<char*>(&index), 4);
		if (file.eof()) {
			break;
		}
		if (current_index != static_cast<int>(index)) {
			// this is a new track
			++ntracks;
			current_index = static_cast<int>(index);
		}
		// read the rest of the values for this position
		char buffer[20];
		file.read(buffer, 20);
  }
  
  // now reset everything
  file.clear();
  current_index = old_current_index;
  file.seekg(current_position);
  
  return ntracks;
}

void Trackfile::SkipNextTrack()
{
  int index = current_index;
  while (current_index == index) {
    float tmpf;
    file.read(reinterpret_cast<char*>(&tmpf), 4);
		if (file.eof()) {
			break;
		}
		index = static_cast<int>(tmpf);
    char buffer[20];
		file.read(buffer, 20);
  }
  
  if (!file.eof()) {
    file.seekg(-4, ios::cur);
    current_index = static_cast<int>(index);
  }
}

void Trackfile::GetNextTrack(deque<float>* x, deque<float>* y, deque<float>* z, deque<float>* t, deque<float>* fake)
{  
  int index = current_index;
  while (current_index == index) {
    float tmpf;
    file.read(reinterpret_cast<char*>(&tmpf), 4);
		if (file.eof()) {
			break;
		}
		index = static_cast<int>(tmpf);
		
    float buffer[5];
		file.read(reinterpret_cast<char*>(buffer), 20);
		x->push_back(buffer[0]);
		y->push_back(buffer[1]);
		z->push_back(buffer[2]);
		t->push_back(buffer[3]);
		fake->push_back(buffer[4]);
  }
  
  float tmpf;
  file.read(reinterpret_cast<char*>(&tmpf), 4);
  if (!file.eof()) {
    file.seekg(-4, ios::cur);
    current_index = static_cast<int>(tmpf);
  }
}
