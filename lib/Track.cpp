/*
 * ----------------------------------------------------------------------------
 * File: Track.cpp
 * ----------------------------------------------------------------------------
 * This is the implementation file for Track objects.
 * ----------------------------------------------------------------------------
 * Created 7/17/03
 * Last updated 7/22/04
 * ----------------------------------------------------------------------------
 * Modified for Yale by NTO on 10/12/11
 * ----------------------------------------------------------------------------
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 */

#include <vector>
#include <Track.h>

using namespace std;

void Track::Add(const Track& t)
{
  pos.insert(pos.end(), t.pos.begin(), t.pos.end());
  time.insert(time.end(), t.time.begin(), t.time.end());
  npoints = pos.size();
}

int Track::Length() const
{
  if (pos.empty()) {
    return 0;
  }
  // calculate effective size of the track: if the track ends with one or
  // more estimated positions, don't count them
  deque<Position>::const_reverse_iterator p_end = pos.rend();
  unsigned int adjust = 0;
  for (deque<Position>::const_reverse_iterator p = pos.rbegin(); 
       p != p_end; ++p, ++adjust) {
    if (!p->IsFake()) {
      break;
    }
  }
	
  return (npoints - adjust);
}

ostream& operator<<(ostream& os, const Track& t)
{
  for (int i = 0; i < t.Length(); ++i) {
    os << t.time[i] << "\t" << t.pos[i] << "\n";
  }

  return os;
}

void Track::PrintEstimates(ostream& os) const
{
  for (int i = 0; i < Length(); ++i) {
    if (pos[i].IsFake()) {
      os << time[i] << "\t" << pos[i] << "\n";
    }
  }
}

int Track::NumFake() const
{
  int count = 0;
  for (int i = 0; i < Length(); ++i) {
    if (pos[i].IsFake()) {
      ++count;
    }
  }

  return count;
}

void Track::WriteGDF(ofstream& output, double index, double fps /* = 1 */) const
{
int len = Length();
	for (int i = 0; i < len; ++i) {
		// Format:
		// Track Index
		// X
		// Y
		// Z
		// Framenumber
		// Camera 1-4 X and Y and Orientation
		// Info
		// Fake bit
		output.write(reinterpret_cast<const char*>(&index), 8);
		double tmp = pos[i].X();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Y();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Z();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = static_cast<double>(time[i]) / fps;
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].X1();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Y1();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Ori1();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].X2();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Y2();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Ori2();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].X3();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Y3();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Ori3();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].X4();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Y4();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Ori4();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = pos[i].Info();
		output.write(reinterpret_cast<const char*>(&tmp), 8);
			if (pos[i].IsFake()) {
				tmp = 1;
			} 
			else {
				tmp = 0;
			}
		output.write(reinterpret_cast<const char*>(&tmp), 8);
	}
}
