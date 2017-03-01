/*
 *  Frame.cpp
 *
 *  Implementation for Frame objects.
 *
 *  Last update: 10/1/09 by NTO
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 *
 */

#include <fstream>

#include <Frame.h>

using namespace std;

ostream& operator<<(ostream& os, const Frame& f)
{
  for (unsigned int i = 0; i < f.pos.size(); ++i) {
    os << "\t" << f.pos[i] << "\n";
  }

  return os;
}
