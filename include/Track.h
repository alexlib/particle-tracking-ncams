/*
 * ----------------------------------------------------------------------------
 * File: Track.h
 * ----------------------------------------------------------------------------
 * This is the header file for Track objects.  A Track object is essentially
 * a fancy frame: it inherits from Frame, but contains more information.
 * ----------------------------------------------------------------------------
 * Created 7/17/03
 * Last updated 3/2/04
 * ----------------------------------------------------------------------------
 * Modified for Yale by NTO on 10/28/11
 *	No longer inherits from Frame
 * ----------------------------------------------------------------------------
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 */

#ifndef TRACK_H
#define TRACK_H

#include <deque>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <Position.h>

class Track {

public:

  // default constructor
  Track();
  // constructor: give one Position argument, and a time (defaults to 0)
  Track(const Position& p, int t = 0);
  // copy-constructor
  Track(const Track& t);
  // destructor: nothing to do
  ~Track() {};

  // Add a new point to the Track
  void Add(const Position& p, int t);
  // Add another Track onto the end of this one
  void Add(const Track& t);

  // return the last point on the Track
  const Position Last() const;
  // return the penultimate point on the Track
  const Position Penultimate() const;
  // return the antepenultimate point on the Track
  const Position Antepenultimate() const;

  // get the length of the Track (not including any ending extrapolated points)
  int Length() const;
  // get the time of a particular element
  int GetTime(int index) const throw(std::out_of_range);

  // check the occlusion counter
  int OcclusionCount() const;
  // increment the occlusion counter
  void Occluded();
  // reset the occlusion counter
  void ResetCounter();

  // find the number of usable fake points on the track
  int NumFake() const;

  // write the track as part of a GDF file (see Tracker.h for format)
  void WriteGDF(std::ofstream& output, double index, double fps = 1) const;
  
  // member operators
  Track& operator=(const Track& t);

  // non-member operators
  friend std::ostream& operator<<(std::ostream& os, const Track& t);

  // print only the estimated points
  void PrintEstimates(std::ostream& os) const;
  
  // free memory
  void Clear();

private:
	std::deque<Position> pos;
  std::deque<int> time;	// the time (as an integer frame number)
  int occluded; // a counter keeping track of the number of frames this 
                // track hasn't had a real particle added to it.
  int npoints;

};

// Inline Function Definitions

inline Track::Track() : occluded(0), npoints(0) {}

inline Track::Track(const Position& p, int t /* = 0 */) 
: occluded(0), npoints(1)
{
	pos.push_back(p);
  time.push_back(t);
}

inline Track::Track(const Track& t)
: pos(t.pos), time(t.time), occluded(t.occluded), npoints(t.npoints)
{}

inline void Track::Add(const Position& p, int t)
{
  pos.push_back(p);
  time.push_back(t);
  ++npoints;
}

inline const Position Track::Last() const
{
  return pos[npoints - 1];
}

inline const Position Track::Penultimate() const
{
  return pos[npoints - 2];
}

inline const Position Track::Antepenultimate() const
{
  return pos[npoints - 3];
}

inline int Track::GetTime(int index) const throw(std::out_of_range)
{
  try {
    return time.at(index);
  } catch (std::out_of_range& e) {
    std::cerr << e.what() << std::endl;
    throw std::out_of_range("Caught out_of_range in Track::GetTime()");
  }
}

inline int Track::OcclusionCount() const
{
  return occluded;
}

inline void Track::Occluded() 
{
  ++occluded;
}

inline void Track::ResetCounter()
{
  occluded = 0;
}

inline Track& Track::operator=(const Track& t)
{
  pos = t.pos;
  time = t.time;
  occluded = t.occluded;
  npoints = t.npoints;

  return *this;
}

inline void Track::Clear()
{
  pos.clear();
  time.clear();
  npoints = 0;
  occluded = 0;
}

#endif /* TRACK_H */
