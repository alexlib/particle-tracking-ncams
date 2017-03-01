/*
* Based on MPIDS/Cornell tracking code, but *heavily modified*
*
* Written 11/1/11 by NTO
*
* Modified by Stefan Kramel 10/15/12
*
* Latest version: 04/09/2013
*	- fixed wrong particlecount in GDF.cpp
*	- removed cout
*
*
* Latest version: 03/28/13
*
*/

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <deque>
#include <vector>

#include <GDF.h>
#include <WesleyanCPV.h>
#include <ParticleFinder.h>
#include <Frame.h>
#include <Calibration.h>
#include <Tracker.h>

using namespace std;

// configuration parameters
struct ConfigFile {
	int ncams;
	deque<string> filenames;
	string setupfile;
	double fps;
	double threshold;
	double cluster_rad;
	int npredict;
	double max_disp;
	int memory;
	int first;
	int last;
	string stereomatched;
	string outname;
};

// globals
deque<Frame> *f;
struct ConfigFile config;

void ImportConfiguration(struct ConfigFile* config, char* name);

int main(int argc, char** argv) {
		if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <configuration file>" << endl;
		exit(1);
	}

	ImportConfiguration(&config, argv[1]);

	// are we trying to track with too much future information?
	if (config.npredict > 2) {
		cerr << "Error: too many predicted frames requested!" << endl;
		exit(1);
	}

	f = new deque<Frame>[config.ncams];

	// read the camera calibration information
	Calibration calib(config.setupfile);
		int first = config.first;
		int last = config.last;
		int threshold = config.threshold;
		double cluster_rad = config.cluster_rad;
		
	// read the data for each camera, find the particles and store them in frameobject f[camid][n]
	for (int camid = 0; camid < config.ncams; ++camid) {
	
		string files = config.filenames[camid];	
		if(files.substr(files.find_last_of(".") + 1) == "cpv") {
			std::cout << camid+1 << " .cpv file(s) detected." << std::endl;
			cout << "Processing CPV file " << files << endl;
			
			WesleyanCPV movie(files, first, last);
			
			int rows = movie.Rows();
			int cols = movie.Cols();
			//int nframes = movie.Frames();
			int** pixels = new int*[rows];
			for (int i = 0; i < rows; ++i) {
				pixels[i] = new int[cols];
			}
			//int avgnum = 0;
					
				for (int n = first; n < last; ++n) {
					cout << "\tReading frame " << n << " of " << last << " in movie " << camid+1 << endl;
					
					//Set pixel array to zero
					for (int i = 0; i < rows; ++i) {
						memset(pixels[i], 0, cols * sizeof(int));
					}
					
					int missed = movie.DecodeNextFrame(pixels, n);
					if (!missed) {
						cout << "push_back Frame: " << n << endl;
						ParticleFinder p(pixels, rows, cols, movie.Colors(), threshold);
						p.Squash(cluster_rad);
						f[camid].push_back(p.CreateFrame());
						//avgnum += (f[camid][n]).NumParticles();
						//cout << "\t(f[camid][n]).NumParticles(): " << (f[camid][n]).NumParticles() << endl;
					}
					else {
						cout << "push_back empty Frame" << endl;
						f[camid].push_back(Frame(Position())); //push_back empty frame
					}
				}
				
			//cout << "\tAveraged " << static_cast<double>(avgnum)/static_cast<double>(nframes) << " particles per frame in movie " << camid+1 << endl;
			
			for (int i = 0; i < rows; ++i) {
				delete []pixels[i];
			}
			delete []pixels;
		} 
		else if(files.substr(files.find_last_of(".") + 1) == "gdf") {
			std::cout << camid+1 << " .gdf file(s) detected." << std::endl;
			cout << "Processing GDF-file " << files << endl;
			
			//Read header information and seek to first frame
			GDF g(files);
			first = g.seekGDF(first);
            
			for (int n = first; n <= last; ++n) {
				cout << "\tReading frame " << n << " of " << last << " in GDF-file " << camid+1 << endl;	
				int missed = g.readGDF2D(n); //read in frame and find wrong and missing frame(s)
				if (!missed) {
				cout << "\tpush_back Frame: " << n << endl;
				f[camid].push_back(g.CreateFrame());
				}
				else {
				cout << "\tBad Frame" << endl;
				f[camid].push_back(Frame()); //push_back empty frame
				}
			}
		}
		else { 
			cerr << "file format unknown." << std::endl;
			exit(1);
		}
	}
	
	// do the stereomatching
	cout << "Stereomatching..." << endl;			
	calib.writeGDFHeader(config.stereomatched);
	
	vector<Frame> matched;
	int nr = 0;

	for (int i = 0; i < (last - first); ++i) {
		cout << "\tProcessing frame " << first+i << " of " << last << endl;
		deque<Frame> toMatch;
		for (int camid = 0; camid < config.ncams; ++camid) {
			//cout << "f[camid][i]: " << f[camid][i] << endl;			
			toMatch.push_back(f[camid][i]);
		}
		Frame all(calib.Stereomatch(toMatch,i));
		nr += all.end()-all.begin();
		cout << "\tCurrent Frame Number = " << i << "; nr = " << nr << endl;	
		matched.push_back(all);
	}
	
	cout << "\tTotal number of stereomatched particles: " << nr << endl;
	// first argument: total number of particles; second argument: number of columns in .gdf file;
	// framenumber, x, y, z, intersect, xy+ori, xy+ori, xy+ori, xy+ori;
	calib.fixHeader(nr,5+3*config.ncams);
    		
	// finally, do the tracking
	cout << "Tracking..." << endl;
	Tracker::TrackMode mode = Tracker::FRAME3;

	if (config.npredict == 0) {
		mode = Tracker::FRAME2;
	} 
	else if (config.npredict == 1) {
		mode = Tracker::FRAME3;
	} 
	else if (config.npredict == 2) {
		mode = Tracker::FRAME4;
	}

	Tracker t(mode, config.max_disp, config.memory, config.fps, 
	config.outname);
	t.MakeTracks(matched);

	delete []f;
	// Done!
	cout << "Done." << endl;
		
	return 0;
}

void ImportConfiguration(struct ConfigFile* config, char* name) {
		cout << "Reading configuration file..." << endl;
		ifstream file(name, ios::in);
		string line;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->ncams = atoi(line.c_str());

		for (int i = 0; i < config->ncams; ++i) {
			getline(file, line);
			line.erase(line.find_first_of(' '));
			config->filenames.push_back(line);
		}

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->setupfile = line;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->fps = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->threshold = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->cluster_rad = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->npredict = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->max_disp = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->memory = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->first = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->last = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->stereomatched = line;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->outname = line;
}
