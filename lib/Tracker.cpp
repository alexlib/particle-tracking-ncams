/*
 * ----------------------------------------------------------------------------
 * File: Tracker.cpp
 * ----------------------------------------------------------------------------
 * This is the implementation file for Tracker objects.
 * ----------------------------------------------------------------------------
 * Created 7/17/03
 * updated 7/22/04
 * ----------------------------------------------------------------------------
 * 
 * added the calculation of the statistics of track length
 * updated 08/01/2010 by HX
 * ----------------------------------------------------------------------------
 * Modified for Yale by NTO on 10/12/11
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

#include <iostream>
#include <fstream>
#include <climits>
#include <list>
#include <algorithm>
#include <iterator>

#ifdef TIME
#include <sys/time.h>
#include <time.h>
#endif

#include <Tracker.h>

using namespace std;

Tracker::Tracker(TrackMode m, double md, int mem, double ifps, string name /* = "track" */)
: too_short(0), ntracks(0), ntotalpoints(0), mode(m), outname(name),
max_disp(md), memory(mem), fps(ifps)
{
	// open the output file with a temporary header
	outfile.open(outname.c_str(), ios::out | ios::binary);
	int magic = 82991;
	outfile.write(reinterpret_cast<const char*>(&magic), 4);
	// number of dimensions
  int tmpi = 2;
  outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
  // number of columns
  tmpi = 19;
  outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
  // number of rows: we don't know this yet!
  tmpi = 0;
  outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
  // a 4 means floating point numbers
  tmpi = 5;
  outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
  // number of total points: we don't know this yet!
  tmpi = 0;
  outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
}

void Tracker::MakeTracks(vector<Frame>& f)
{  
#ifdef TIME
  struct timeval t1;
  gettimeofday(&t1, NULL);
#endif
  
  // keep a list of the indices of active tracks: i.e., tracks which can 
  // still have points added to them
  deque<int> activelist;
  
  long trackindex = 0;
  
  // initialize the tracks vector with the values in the first frame
  for (int i = 0; i < f[0].NumParticles(); ++i) {
    tracks[i] = new Track(f[0][i]);
    // also add these tracks to the active track list
    activelist.push_back(i);
    ++trackindex;
  }
  
  // the frame number we're on: it will be one when we start the loop
  int framenum = 1;
  
  // loop over all the frames
  
  // precalculate the end of the iterator...
  vector<Frame>::iterator fr_end = f.end();
  // ...but we don't want to go all the way to the end
  switch (mode) {
    case FRAME3:
      fr_end -= 1;
      break;
    case FRAME4:
      fr_end -= 2;
      break;
    default:
      break;
  }
  for (vector<Frame>::iterator fr0 = f.begin(); fr0 != fr_end; ++fr0, ++framenum) {
#ifdef DEBUG
    cout << "Processing frame number " << framenum << endl;
#endif
    // this is the next frame:
    vector<Frame>::iterator fr1 = fr0 + 1;
    
    // make sure there are particles in the next frame!
    if (fr1->NumParticles() <= 0) {
      // there aren't; skip to the next timestep
      PadTracks(activelist, framenum);
      continue;
    }
    
    // allocate a 1D array that will store the best matches to the given
    // particle in the next frame
    float* costs = new float[fr1->NumParticles()];
    int* links = new int[fr1->NumParticles()];
    for (int i = 0; i < fr1->NumParticles(); ++i) {
      links[i] = UNLINKED;
    }
    
    // make the links!
    MakeLinks(activelist, *fr1, *(fr1+1), costs, links);
    
    // now add the best matches to the appropriate Tracks
    long n_new_tracks = 0;
    size_t n_ended_tracks = activelist.size();
    for (int i = 0; i < fr1->NumParticles(); ++i) {
      if (links[i] == UNLINKED) {
        // this particle should start a new track
        tracks[trackindex] = new Track((*fr1)[i], framenum);
        activelist.push_back(trackindex);
				++trackindex;
        ++n_new_tracks;
      } else {
        // add this particle to an existing track
        tracks.find(links[i])->second->Add((*fr1)[i], framenum);
        tracks.find(links[i])->second->ResetCounter();
        --n_ended_tracks;
      }
    }
    
    // deal with tracks that didn't get added to: if the track has not been 
    // added to for a time greater than memory, remove the track
    // from the active track list.  otherwise, make an estimate of the next 
    // point on the track and store that, updating the track's occlusion 
    // counter.
    PadTracks(activelist, framenum);
    
    // free memory
    delete []costs;
    delete []links;
    
    // print diagnostic information
    cout << "Processed frame " << framenum << endl;
    cout << "\tNumber of particles found: " << fr1->NumParticles() << endl;
		cout << "\t Number of active tracks: " << activelist.size() << endl;
		cout << "\t Number of new tracks started here: " << n_new_tracks << endl;
		cout << "\t Number of tracks that found no match: " << n_ended_tracks << endl;
		cout << "\t Total number of tracks: " << tracks.size() << endl;
  }
  
  // Write the rest of the tracks out, freeing their memory as we go
	deque<int>::const_iterator tr_end = activelist.end();
	for (deque<int>::const_iterator tr = activelist.begin(); tr != tr_end; ++tr) {
		Track* t = tracks.find(*tr)->second;
		if (t->Length() >= MINTRACK) {
		  t->WriteGDF(outfile, ntracks, fps);
		  ++ntracks;
		  ntotalpoints += t->Length();
		}
		// free memory
		delete t;
	}
	// now fix up the header with the proper sizes
	outfile.seekp(12, ios::beg);
	outfile.write(reinterpret_cast<const char*>(&ntotalpoints), 4);
	outfile.seekp(4, ios::cur);
	int tmpi = 19 * ntotalpoints;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
  
#ifdef TIME
  struct timeval t2;
  gettimeofday(&t2, NULL);
  cerr << "Time for tracking: " << (t2.tv_sec - t1.tv_sec) + 1e-6 * 
  (t2.tv_usec - t1.tv_usec) << endl;
#endif
}

void Tracker::PadTracks(deque<int>& activelist, int framenum)
{
	list<int> writelist;
	deque<int> stillactive;
	
  for (deque<int>::iterator tr = activelist.begin(); tr != activelist.end(); ++tr) {
		Track* t = tracks.find(*tr)->second;
		int len = t->Length();
		
    // if this track was just added to, skip it
		if (t->GetTime(len - 1) == framenum) {
			stillactive.push_back(*tr);
      continue;
    }
    // check the occlusion counter of this track.
		if (t->OcclusionCount() >= memory) {
			
			// if this track is too short, free its memory
			if (len < MINTRACK) {
				delete t;
				tracks.erase(*tr);
				++too_short;
			} else {
				// otherwise, we'll want to write it out
				writelist.push_back(*tr);
			}
    } else {
      // add an estimated position to this track, and increment its 
      // occlusion counter
			
      // make one more check: if this track was only one or two points long
      // and hasn't found another one, drop it
      if (len <= 2) {
				delete t;
				tracks.erase(*tr);
				++too_short;
				continue;
      }
			
      // don't bother with the times in our estimate -- they'll cancel out.
			const Position last = t->Last();
			const Position penultimate = t->Penultimate();
			const Position antepenultimate = t->Antepenultimate();
			
			Position velocity = last - penultimate;
			Position acceleration = 0.5 * (last - 2.0 * penultimate + antepenultimate);
			Position estimate = last + velocity + 0.5 * acceleration;
			if (Distance(estimate, last) > (max_disp * max_disp)) {
				estimate = last;
			}
			
      // keep track of the fact that this position is an estimate
      estimate.SetFake();
      t->Add(estimate, framenum);
			
      t->Occluded();
			
			// and remember that this is still an active track
			stillactive.push_back(*tr);
    }
  }
	
	// write out the tracks that have ended and are long enough, then free their memory
	list<int>::const_iterator tr_end = writelist.end();
	for (list<int>::const_iterator tr = writelist.begin(); tr != tr_end; ++tr) {
		Track* t = tracks.find(*tr)->second;
		t->WriteGDF(outfile, ntracks, fps);
		++ntracks;
		ntotalpoints += t->Length();
		// free memory
		delete t;
		tracks.erase(*tr);
	}
	
	activelist = stillactive;
}

void Tracker::MakeLinks(const deque<int>& activelist, Frame& fr1, Frame& fr2, 
                        float*& costs, int*& links)
{
  // loop over all the tracks on the activelist
  deque<int>::const_iterator tr_end = activelist.end();
  int listcount = 0;
  for (deque<int>::const_iterator tr = activelist.begin(); tr != tr_end; ++tr, ++listcount) {
    // the current track
		Track* t = tracks.find(*tr)->second;
    // the current position
    Position now = t->Last();
		
		int len = t->Length();
		
    // we'll need a velocity and an estimated future position
    Position velocity;
    Position estimate;
    
    // does the current track have more than one point? Or are we using 
    // nearest neighbor search?
    if (len == 1 || mode == FRAME2) {
      estimate = now;
    } else {
			// this track was at least two points long; use it to get an estimate of the velocity
			velocity = now - t->Penultimate();
			// if the track contains more than two particles, also calculate
			// an acceleration to help with the estimate
			if (len > 2) {
				// this track contains multiple particles
				Position acceleration = now - 2.0 * t->Penultimate() + t->Antepenultimate();
				estimate = now + velocity + 0.5 * acceleration;
			} else {
				// this track doesn't contain multiple particles, so just use the 
				// velocity to estimate a position
				estimate = now + velocity;
			}
    }
    
    pair<int, float> cost;
    if (mode == FRAME4) {
      cost = ComputeCost(fr1, fr2, estimate, velocity, now, false);
    } else {
      cost = ComputeCost(fr1, fr1, estimate, velocity, now, true);
    }
    if (cost.first == -1) {
			// no matches found!
			continue;
		}
		if (links[cost.first] == UNLINKED || costs[cost.first] > cost.second) {
		  costs[cost.first] = cost.second;
		  links[cost.first] = *tr;
		}
  }
}

pair<int, float> Tracker::ComputeCost(Frame& fr1, Frame& fr2, 
                                      const Position& estimate, 
                                      const Position& velocity, 
                                      const Position& now, bool stopflag)
{
  // find possible continuations in the next frame.
  
  // find possible matches
  deque<int> matches;
  deque<float> match_costs;
  float mincost = 1e6;
  
	Frame::const_iterator fitend = fr1.end();
	for (Frame::const_iterator fit = fr1.begin(); fit != fitend; ++fit) {
    float mag = Distance(estimate, *fit);
    if (mag > max_disp * max_disp) {
      continue;
    }
    matches.push_back(fit.where());
    if (stopflag == true) {
      // don't look into the future any more
      match_costs.push_back(mag);
      if (mincost > mag) {
        mincost = mag;
      }
    } else {
      // project again!
      Position new_velocity = *fit - now;
			Position acceleration = new_velocity - velocity;
			Position new_estimate = *fit + new_velocity + 0.5 * acceleration;
      pair<int, float> cost = ComputeCost(fr2, fr2, new_estimate, new_velocity, *fit, true);
      match_costs.push_back(cost.second);
      if (mincost > cost.second) {
        mincost = cost.second;
      }
    }
  }
  
  pair<int, float> retval = make_pair(-1, 1e6);
	deque<int>::const_iterator m_end = matches.end();
	deque<float>::const_iterator mc = match_costs.begin();
	for (deque<int>::const_iterator m = matches.begin(); m != m_end; ++m, ++mc) {
		if (*mc > mincost) {
			continue;
		}
		retval = make_pair(*m, *mc);
	}
	
	return retval;
}
