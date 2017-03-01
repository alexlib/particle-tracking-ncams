/*
 *  GDF.cpp
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
#include <queue>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <GDF.h>
#include <Position.h>
#include <WesleyanCPV.h>

using namespace std;

GDF::GDF(std::string filename) throw(out_of_range)
: outname(filename)
{
// Read Header:
    infile.open(outname.c_str(), ios::in | ios::binary);
	if (infile.is_open()){
		infile.read(reinterpret_cast<char*>(&magic), 4);
		// number of dimensions
		infile.read(reinterpret_cast<char*>(&tmpi), 4);
		// number of columns
		infile.read(reinterpret_cast<char*>(&cols), 4);
		// number of rows
		infile.read(reinterpret_cast<char*>(&rows), 4);
		// 4 means floating point numbers (see Matlab read_gdf function for more info)
		infile.read(reinterpret_cast<char*>(&tmpi), 4);
		// number of total points
		infile.read(reinterpret_cast<char*>(&tmpi), 4);
	}
    waiting_to_be_written = 0;
}

int GDF::seekGDF(int start) {
    while(!waiting_to_be_written) {
        filePos[first] = infile.tellg();
        infile.seekg(40, ios::cur);
        infile.read(reinterpret_cast<char*>(&fi), 8);
        currentFrameNum = fi;
        while(currentFrameNum == fi && !infile.eof()){
            filePos[tmp] = infile.tellg();
            infile.seekg(40, ios::cur);
            infile.read(reinterpret_cast<char*>(&fi), 8);
            if (infile.eof()) {
                std::cout << "\tEnd of file reached during seeking" << endl;
                exit(0);
            }
            nextFrameNum = fi;
        }
        if (prevFrameNum == currentFrameNum-1 || nextFrameNum == currentFrameNum+1) {
            if (currentFrameNum >= start) {
                cout << "\tFirst good frame number: " << currentFrameNum << endl;
                startFrameNum = currentFrameNum;
                infile.seekg(filePos[first], ios::beg);
                waiting_to_be_written = 1;
            }
            prevFrameNum = currentFrameNum;
        }
        else {
            infile.seekg(filePos[tmp], ios::beg);
        }
    }
    waiting_to_be_written = 0;
    return(startFrameNum);
}

int GDF::readGDF2D(int frame) {
	deque<double> newx;
	deque<double> newy;
	deque<double> newori;

    filePos[current] = infile.tellg();
    infile.seekg(32, ios::cur);
    infile.read(reinterpret_cast<char*>(&fi), 8);
    int particlecount = fi;
    infile.read(reinterpret_cast<char*>(&fi), 8);
    currentFrameNum = fi;
    infile.seekg(filePos[current], ios::beg);
    cout << "\tCurrent Frame Number: " << currentFrameNum << endl;
    cout << "\t" << particlecount << " particle(s) found" << endl;
    

    if (currentFrameNum==frame) {
        for(int i = 0; i < particlecount; i++){
                    infile.read(reinterpret_cast<char*>(&xi), 8);
                    newx.push_back(xi);
                    infile.read(reinterpret_cast<char*>(&yi), 8);
                    newy.push_back(yi);
                    infile.seekg(8, ios::cur); //skip brightness bytes
                    infile.read(reinterpret_cast<char*>(&orii), 8);
                    newori.push_back(orii);
                    infile.seekg(8, ios::cur); //skip number of particles bytes
                    infile.read(reinterpret_cast<char*>(&fi), 8);
                    currentFrameNum = fi;
                    if (currentFrameNum!=frame) {
                        cout << "\tIncorrect particlecount at frame" << currentFrameNum << endl;
                        break;
                    }
                    x = newx;
                    y = newy;
                    ori = newori;	
                }
        missedFrame = 0;
    }
    else {
        missedFrame = 1;
    }
    return(missedFrame);
}

void GDF::fixHeader(int nr, int cols){
	
	// now fix up the header with the proper sizes
	outfile.seekp(8, ios::beg);
	outfile.write(reinterpret_cast<const char*>(&cols), 4);
	outfile.write(reinterpret_cast<const char*>(&nr), 4);
	outfile.seekp(4, ios::cur);
	int tmpi = cols * nr;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);

}

Frame GDF::CreateFrame() {
  deque<Position> pos;
  deque<double>::const_iterator xi_end = x.end();
  deque<double>::const_iterator yi = y.begin();
  deque<double>::const_iterator orii = ori.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi, ++orii) {
  	//cout << "Create Frame: " << *xi << " " << *yi <<" "<< 0 <<" "<< *orii << endl;
    pos.push_back(Position(*xi, *yi, 0, *orii));
  }
	return Frame(pos);
}
