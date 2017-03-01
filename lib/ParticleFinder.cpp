/*
 *  ParticleFinder.cpp
 *
 *  This is the implementation file for ParticleFinder objects.
 *
 *  Last update: 10/12/11 by NTO
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 *
 */

#include <iostream>
#include <fstream>
#include <queue>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <ParticleFinder.h>
#include <Logs.h>
#include <Position.h>

using namespace std;

ParticleFinder::ParticleFinder(int**& p, int rows, int cols, int depth, int threshold) throw(out_of_range)
: pixels(p)
{
	int colors = depth;
  
  // walk through the given array, skipping the first and last row and column
  for (int i = 1; i < (rows - 1); ++i) {
    for (int j = 1; j < (cols - 1); ++j) {
      // is this pixel a local maximum above threshold?
      if ((pixels[i][j] >= threshold) && IsLocalMax(i, j)) {
				// read in the local maximum pixel value as well as the values to its
				// left and right and top and bottom in order to calculate the 
				// particle center.  note: add 0.5 to row and column values to put
				// the pixel origin in its center
				double x1 = (j - 1) + 0.5;
				double x2 = j + 0.5;
				double x3 = (j + 1) + 0.5;
				double y1 = (i - 1) + 0.5;
				double y2 = i + 0.5;
				double y3 = (i + 1) + 0.5;
				
				// check the pixel values to make sure we have no corrupted memory
				if ((abs(pixels[i][j]) > colors) ||
						(abs(pixels[i - 1][j]) > colors) ||
						(abs(pixels[i + 1][j]) > colors) ||
						(abs(pixels[i][j - 1]) > colors) ||
						(abs(pixels[i][j + 1]) > colors)) {
					// this is a serious problem: we can't continue
					throw out_of_range("Pixel out of range!");
				}
				
				// find the column value, moving 0 intensities to 0.0001 so we can 
				// take their log
				double lnz1, lnz2, lnz3;
				if (colors == 255) {
					if (pixels[i][j - 1] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = Logs::log8bit[pixels[i][j - 1]];
					}
					if (pixels[i][j] == 0) {
            			lnz2 = log(0.0001);
					} else {
						lnz2 = Logs::log8bit[pixels[i][j]];
					}
					if (pixels[i][j + 1] == 0) {
						lnz3 = log(0.0001); 
					} else {
						lnz3 = Logs::log8bit[pixels[i][j + 1]];
					}
        } else if (colors == 65535) {
          lnz1 = Logs::log16bit[pixels[i][j - 1]];
          lnz2 = Logs::log16bit[pixels[i][j]];
					lnz3 = Logs::log16bit[pixels[i][j + 1]];
				} else {
          if (pixels[i][j - 1] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = log(static_cast<double>(pixels[i][j - 1]));
					}
					if (pixels[i][j] == 0) {
            lnz2 = log(0.0001);
					} else {
						lnz2 = log(static_cast<double>(pixels[i][j]));
					}
					if (pixels[i][j + 1] == 0) {
						lnz3 = log(0.0001); 
					} else {
						lnz3 = log(static_cast<double>(pixels[i][j + 1]));
					}					
        }
				
				double xc = -0.5 * ((lnz1 * ((x2 * x2) - (x3 * x3))) - (lnz2 * ((x1 * x1) - (x3 * x3))) + (lnz3 * ((x1 * x1) - (x2 * x2)))) / ((lnz1 * (x3 - x2)) - (lnz3 * (x1 - x2)) + (lnz2 * (x1 - x3)));
				
				// were these numbers valid?
				if (!finite(xc)) {
          // no -- we had a problem.  drop this particle
					continue;
				}
				
				// find the row value
				if (colors == 255) {
          if (pixels[i - 1][j] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = Logs::log8bit[pixels[i - 1][j]];
					}
					if (pixels[i + 1][j] == 0) {
            lnz3 = log(0.0001);
					} else {
						lnz3 = Logs::log8bit[pixels[i + 1][j]];
					}
				} else if (colors == 65535) {
					lnz1 = Logs::log16bit[pixels[i - 1][j]];
					lnz3 = Logs::log16bit[pixels[i + 1][j]];
				} else {
          if (pixels[i - 1][j] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = log(static_cast<double>(pixels[i - 1][j]));
					}
					if (pixels[i + 1][j] == 0) {
            lnz3 = log(0.0001);
					} else {
						lnz3 = log(static_cast<double>(pixels[i + 1][j]));
					}
				}
				
				double yc = -0.5 * ((lnz1 * ((y2 * y2) - (y3 * y3))) - (lnz2 * ((y1 * y1) - (y3 * y3))) +
														(lnz3 * ((y1 * y1) - (y2 * y2)))) 
														/ ((lnz1 * (y3 - y2)) - (lnz3 * (y1 - y2)) + (lnz2 * (y1 - y3)));
				
				// check these numbers too
				if (!finite(yc)) {
					// a problem occurred, so we'll drop these numbers and move on
					continue;
				}
				//cout << "xc: " << xc << " yc: " << yc << endl;
			x.push_back(xc);
			y.push_back(yc);
      }
    }
  }
}

void ParticleFinder::WriteToFile(string filename) {
  // now we have all the particle centers and can write them to a file
  ofstream outfile(filename.c_str(), ios::out);
  outfile << "# Modified Gaussian fitting (3-point method)\n";
	
  deque<double>::const_iterator xi_end = x.end();
  deque<double>::const_iterator yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    outfile << *xi << "\t" << *yi << "\n";
  }
	
  outfile.close();
}

Frame ParticleFinder::CreateFrame() {
	deque<Position> pos;
	deque<double>::const_iterator xi_end = x.end();
	deque<double>::const_iterator yi = y.begin();
	for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
		pos.push_back(Position(*xi, *yi, 0));
	}
	return Frame(pos);
}

void ParticleFinder::Squash(double rad) {
  if (rad < 1) {
    // don't do anything for small cluster radii, for efficiency
    return;
  }
  
  deque<deque<double>::const_iterator> bad;
  deque<double> newx;
  deque<double> newy;
  
  deque<double>::const_iterator xi_end = x.end();
  deque<double>::const_iterator yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    // have we looked at this entry before?
    if (find(bad.begin(), bad.end(), xi) != bad.end()) {
      // yes: skip it.
      continue;
    }
    
    double avgx = *xi;
    double avgy = *yi;
    int N = 1;
    // loop over other positions, looking for neighbors
    deque<double>::const_iterator xi2 = xi + 1;
    deque<double>::const_iterator yi2 = yi + 1;
    for (; xi2 != xi_end; ++xi2, ++yi2) {
      if (pow(*xi-*xi2,2)+pow(*yi-*yi2,2) > rad*rad) {
        continue;
      }
      // these positions are nearby
      avgx += *xi2;
      avgy += *yi2;
      ++N;
    }
    // did we find a cluster?
    if (N == 1) {
      // nope.
      continue;
    }
    
    bad.push_back(xi);
    int oldN = N;
    
    // now loop through again, looking for everything near *this average*
    // location, and re-compute the average. continue until we converge.
    while (true) {
      double ax = avgx / static_cast<double>(N);
      double ay = avgy / static_cast<double>(N);
      avgx = *xi;
      avgy = *yi;
      N = 1;
      deque<double>::const_iterator xi2 = xi + 1;
      deque<double>::const_iterator yi2 = yi + 1;
      for (; xi2 != xi_end; ++xi2, ++yi2) {
        if (pow(ax-*xi2,2)+pow(ay-*yi2,2) > rad*rad) {
          continue;
        }
        // these positions are nearby
        avgx += *xi2;
        avgy += *yi2;
        ++N;
      }
      if (oldN == N) {
        // we converged! now do bookkeeping
        newx.push_back(avgx / static_cast<double>(N));
        newy.push_back(avgy / static_cast<double>(N));
        deque<double>::const_iterator xi2 = xi + 1;
        deque<double>::const_iterator yi2 = yi + 1;
        for (; xi2 != xi_end; ++xi2, ++yi2) {
          if (pow(ax-*xi2,2)+pow(ay-*yi2,2) <= rad*rad) {
            bad.push_back(xi2);
          }
        }
        break;
      } else {
        // try again
        oldN = N;
      }
    }
  }
  
  // finally, find the good positions left in the original queue
  yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    if (find(bad.begin(), bad.end(), xi) == bad.end()) {
      newx.push_back(*xi);
      newy.push_back(*yi);
      //cout << *xi << "\t" << *yi << endl;
    }
  }
  x = newx;
  y = newy;
}

bool ParticleFinder::IsLocalMax(int r, int c)
{
  int val = pixels[r][c];

  // check only the 4 neighbors directly above, below, to the left, and to the right
  if ((pixels[r][c - 1] > val) || (pixels[r][c + 1] > val) || (pixels[r - 1][c] > val) || (pixels[r + 1][c] > val)) {
    return false;
  }
	
  return true;
}
