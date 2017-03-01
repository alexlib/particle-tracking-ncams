/*
 * ----------------------------------------------------------------------------
 * File: Tracker.h
 * ----------------------------------------------------------------------------
 * This is the header file for Tracker objects.  A Tracker object stores a 
 * collection of Frames, and uses them to calculate particle tracks.  It then
 * stores the tracks in one or more files.  
 * ----------------------------------------------------------------------------
 * Created 7/17/03
 * updated 6/8/04
 * ____________________________________________________________________________
 * 
 * added the calculation of the statistics of track length
 * updated 08/01/2010 by HX
 *
 * ----------------------------------------------------------------------------
 * Modified for Yale tracking package by NTO on 10/12/11
 * New output format:
 
 // FORMAT:
 //   HEADER:
 //   magic number: 82991										(4-byte int)
 //   number of array dimensions: 2					(4-byte int)
 //   number of data fields per point: 5     (4-byte int)
 //   number of data points									(4-byte int)
 //   data type code: 4											(4-btye int)
 //   total number of data fields						(4-byte int)
 // 
 //   EACH POINT:
 //   track index (unique for each track)		(4-byte float)
 //   X position                             (4-byte float)
 //   Y position                             (4-byte float)
 //   Z Position                             (4-byte float)
 //   time                                   (4-byte float)
 //   estimate bit (0 if the point is real,	(4-byte float)
 //  							 1 if extrapolated)
 
 *
 * ----------------------------------------------------------------------------
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 */

#ifndef TRACKER_H
#define TRACKER_H

#include <vector>
#include <deque>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <utility>

#include <Frame.h>
#include <Track.h>
#include <Position.h>

class Tracker {

public:
  // data types
  
  // this enum defines the mode of the Tracker: how many frames to use, 
  // the criterion for determining the correct track, etc.
  enum TrackMode {
    FRAME2,
    FRAME3,
    FRAME4
  };

  // constructor
  Tracker(TrackMode m, double md, int mem, double ifps, 
          std::string name = std::string("track"));
  // destructor: nothing to do
  ~Tracker();

  // do the work of making all the tracks
  void MakeTracks(std::vector<Frame>& f);

private:
  typedef std::map<int, Track*> TrackMap;
  
  static const int UNLINKED = -1;
  static const int MINTRACK = 10;
  static const int EMPTY = INT_MAX;
  
  int too_short;
  double ntracks;
  int ntotalpoints;
  
  TrackMap tracks;
  TrackMode mode;
  std::string outname;
  std::ofstream outfile;
  
  double max_disp;
  int memory;
  double fps;
  
  // helper functions

  // extend the tracks that weren't added to in this frame pair by 
  // extrapolation
  void PadTracks(std::deque<int>& activelist, int framenum);
  
  // generate the proper frame-to-frame links (replacing class LinkMatrix)
  void MakeLinks(const std::deque<int>& activelist, Frame& fr1, Frame& fr2, 
                 float*& costs, int*& links);
  
  // compute the cost function for a possible link
  std::pair<int, float> ComputeCost(Frame& fr1, Frame& fr2, 
                                    const Position& estimate, 
                                    const Position& velocity,
                                    const Position& now, bool stopflag);

};

// Inline Function Definitions
inline Tracker::~Tracker()
{
  outfile.close();
}

#endif /* TRACKER_H */
