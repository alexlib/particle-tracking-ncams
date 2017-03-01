/*
 *  Calibration.cpp
 *  
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 *
 */

#include <sstream>
#include <iostream>
#include <algorithm>
#include <list>
#include <vector>

#include <Calibration.h>
#include <GDF.h>
#include <Position.h>

using namespace std;

Calibration::Calibration(std::string& fname)
{
	// remove comments from the file
	ifstream infile(fname.c_str(), ios::in);
	string line;
	stringstream parsed;
	while (getline(infile, line)) {
		size_t commentpos = line.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) {
				line.erase(commentpos);
			}
			parsed << line << '\t';
		}
	}
	infile.close();
	
	// read the number of cameras
	parsed >> ncams;
	
	// now read in the parameters for each camera
	for (int i = 0; i < ncams; ++i) {
		cams.push_back(Camera(parsed));
	}
	
	// read in the tolerances
	parsed >> mindist_2D >> mindist_3D;
}

void Calibration::writeGDFHeader(std::string outname)
{
	// open the output file with a temporary header
	outfile.open(outname.c_str(), ios::out | ios::binary);
	int magic = 82991;
	outfile.write(reinterpret_cast<const char*>(&magic), 4);
	// number of dimensions
	int tmpi = 2;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// number of columns
	tmpi = 17;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// number of rows: we don't know this yet!
	tmpi = 0;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// a 4 means single presicion, 5 double precision
	tmpi = 5;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// number of total points: we don't know this yet!
	tmpi = 0;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);

	cout << "\nHeader information written..." << endl;

}

void Calibration::fixHeader(int nr, int cols){
	
	// now fix up the header with the proper sizes
	outfile.seekp(8, ios::beg);
	outfile.write(reinterpret_cast<const char*>(&cols), 4);
	outfile.write(reinterpret_cast<const char*>(&nr), 4);
	outfile.seekp(4, ios::cur);
	int tmpi = cols * nr;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	cout << "\nHeader information updated!" << endl;
}

Frame Calibration::Stereomatch(const deque<Frame>& iframes, int framenumber) throw(runtime_error)
{
    if (iframes.size() != cams.size()) {
        throw runtime_error("Number of cameras and number of images do not match!");
    }
    
    // step 1.
    // correct all the distortions, and move the points into a coordinate
    // system (in mm) with the origin in the middle of the image plane.
    deque<Frame> corrframes;
    cout << "\tCorrecting distortion..." << endl;
    for (int i = 0; i < ncams; ++i) {
        deque<Position> corrpos;
        Frame::const_iterator fitend = iframes[i].end();
        for (Frame::const_iterator fit = iframes[i].begin(); fit != fitend; ++fit) {
            corrpos.push_back(cams[i].UnDistort(*fit));
        }
        corrframes.push_back(Frame(corrpos));
    }
    
    // step 2.
    // for each camera, draw a line of sight through each particle on its image plane.
    // project these lines of sight onto the image planes of each other camera.
    // particles within mindist_2D of these lines are candidate matches.
    cout << "\tConstructing pair lists..." << endl;
    int avgsize = 0;
    int numlists = 0;
    mcam = -1;
    // nasty data structure; is there a better way to do this?
    list<Frame::const_iterator> ***pairlists = new list<Frame::const_iterator>**[ncams];
    for (int i = 0 ; i < ncams; ++i) {
        pairlists[i] = new list<Frame::const_iterator>*[corrframes[i].NumParticles()];
        // loop over particles in frame i
        Frame::const_iterator pAend = corrframes[i].end();
//         cout << "\tList size, corrframes[i].NumParticles(): " << corrframes[i].NumParticles() << endl;
        for (Frame::const_iterator pA = corrframes[i].begin(); pA != pAend; ++pA) {
            // this particle's coordinates in world space
//             cout << "\tChecking particle[" << pA.where() << "/" << corrframes[i].NumParticles()-1 << "] at [" << cams[i].Distort(*pA).X() << ", " << cams[i].Distort(*pA).Y() << "] on camera " << i << endl;
            Position pAworld(cams[i].ImageToWorld(*pA));
            pairlists[i][pA.where()] = new list<Frame::const_iterator>[ncams];
            for (int k = 0; k < ncams; ++k) {
                if (i == k) {
                    continue;
                }
                // position of camera i's projective center on camera k
                Position center(cams[k].WorldToImage(cams[i].Center())); /* DIFF */
                // position of this particle on camera k
                Position particle(cams[k].WorldToImage(pAworld));
                // unit vector in (projected) line of sight direction
                Position lineofsight(particle - center);
                lineofsight /= lineofsight.Magnitude();
                // unit vector normal to the line of sight
                Position perpdir(lineofsight.Y(), -lineofsight.X(), 0);
                
                // now loop over the particles in frame k
                Frame::const_iterator pBend = corrframes[k].end();
                for (Frame::const_iterator pB = corrframes[k].begin(); pB != pBend; ++pB) {
                    // if the distance from camera i's projective center along perpdir
                    // is less than mindist_2D, this is a potential match!
                    Position pBline(*pB - center);
                    if (abs(Dot(pBline, perpdir)) < mindist_2D) {
                        pairlists[i][pA.where()][k].push_back(pB);
                        Position pB_check(*pB);
//                         cout << "\t\tpairlists[" << i << "][" << pA.where() << "][" << k << "]: Particle[" << pA.where() << "] found a possible match on camera " << k << " at [" << cams[k].Distort(pB_check).X() << ", " << cams[k].Distort(pB_check).Y() << "] distance: " << abs(Dot(pBline, perpdir)) << endl;
                    }
                }
//                 cout << "\t\t" << pairlists[i][pA.where()][k].size() << " possible match(es) on camera " << k << endl;
                list<Frame::const_iterator>::const_iterator pB_check_end = pairlists[i][pA.where()][k].end();
                for (list<Frame::const_iterator>::const_iterator pB_check = pairlists[i][pA.where()][k].begin(); pB_check != pB_check_end; ++pB_check) {
//                     cout << "\t\tpairlists[" << i << "][" << pA.where() << "][" << k << "]: Particle[" << (*pB_check).where() << "] found a possible match on camera " << k << " at [" << cams[k].Distort(**pB_check).X() << ", " << cams[k].Distort(**pB_check).Y() << "]" << endl;
                }
                if (pairlists[i][pA.where()][k].size()==0) {
//                     cout << "\t\tpairlists[" << i << "][" << pA.where() << "][" << k << "]: Particle[" << pA.where() << "] found no match on camera " << k << endl;
                }
                
                avgsize += pairlists[i][pA.where()][k].size();
                numlists++;
            }
        }
    }

	cout << "\tMean pairlist size: " << static_cast<double>(avgsize)/static_cast<double>(numlists) << endl;
	
	// step 3.
	// go through the lists and search for consistency; that is, particles that show up 
	// on each other's lists. these are "real" matches, and we will compute their 3D positions.
	// note: we require the particles to be seen on all cameras -- so we only need to search along 
	// one of the lists. might as well pick cam 0!
    vector < deque<Position> > PosTouse;
	deque<Position> matchedPos;
	deque< deque<int> > frame_indices;
	deque<double> raydists;
	cout << "\tPerforming consistency checks..." << endl;
	Frame::const_iterator p0end = corrframes[0].end();
	for (Frame::const_iterator p0 = corrframes[0].begin(); p0 != p0end; ++p0) {		
		// quick check: does p0 have any other particles on all of its lists?
// 		bool toNextp0 = false;
        int ncams_missing = 0;
        for (int i = 1; i < ncams; ++i) {
            if (pairlists[0][p0.where()][i].size() == 0) {
//                 cout << "\t\tParticle[" << p0.where() << "] p0_check: no match on camera " << i << endl;
                ncams_missing += 1;
//                 toNextp0 = true;
//                 break;
            }
		}
//         cout << "\t\tParticle[" << p0.where() << "] p0_check: no match on " << ncams_missing << " cameras..." << endl;
        if (ncams_missing > 0) { // to next p0
			continue;
		}
		vector< deque<Frame::const_iterator> > toMatch;
		deque<Frame::const_iterator> tmp;
		tmp.push_back(p0);
		toMatch.push_back(tmp);
		// loop over the other cameras
		for (int i = 1; i < ncams; ++i) {
			// loop over the partially started point sets
			unsigned int s = toMatch.size();
			for (unsigned int k = 0; k < s; ++k) {
				// is this set too small already?
				if (toMatch[k].size() < static_cast<size_t>(i)) {
					// yes
					continue;
				}
				// now, can we push any points onto the end of this set?
				Frame::const_iterator pend = corrframes[i].end();
				for (Frame::const_iterator p = corrframes[i].begin(); p != pend; ++p) {
					// does this point have anything on its lists?
					bool toNextp = false;
					for (int j = i+1; j < ncams; ++j) {
						if (pairlists[i][p.where()][j].size() == 0) {
							toNextp = true;
							break;
						}
					}
					if (toNextp) {
						continue;
					}
					// loop over the points already in the point set
					bool canAdd = true;
					for (int j = 0; j < i; ++j) {
						// is this point in the point set on p's list?
						list<Frame::const_iterator>::const_iterator onlist_tm_p = find(pairlists[i][p.where()][j].begin(),pairlists[i][p.where()][j].end(),toMatch[k][j]); //(*tm)[j]);
//                         cout << "\t toMatch[k][j]: " << *toMatch[k][j] << endl;
						if (onlist_tm_p != pairlists[i][p.where()][j].end()) {
							// yes! is the converse true?
							list<Frame::const_iterator>::const_iterator onlist_p_tm = find(pairlists[j][toMatch[k][j].where()][i].begin(),pairlists[j][toMatch[k][j].where()][i].end(),p);
							if (onlist_p_tm != pairlists[j][toMatch[k][j].where()][i].end()) {
								// yes! we can potentially add it
								canAdd &= true;
							} else {
								// nope
								canAdd &= false;
							}
						} else {
							// nope
							canAdd &= false;
						}
					}
					if (canAdd) {
						// we passed all the tests. add this point!
						// note: add it as a new queue so that we can add other points too, if we want to.
						deque<Frame::const_iterator> topush = toMatch[k];//*tm;
						topush.push_back(p);
						toMatch.push_back(topush);
					}
				}
			}
		}
		// now, match the point sets that are full!
		vector< deque<Frame::const_iterator> >::const_iterator tmend = toMatch.end();
		for (vector< deque<Frame::const_iterator> >::const_iterator tm = toMatch.begin(); tm != tmend; ++tm) {
			// is this set too small?
			if (tm->size() < static_cast<size_t>(ncams)) {
				// yes
				continue;
			}
			// we can match it!
			deque<Position> PosToMatch;
			deque<int> indices;
			for (int i = 0; i < ncams; ++i) {
				PosToMatch.push_back(*((*tm)[i]));
				indices.push_back((*tm)[i].where());
			}
			pair<double,Position> wpos = WorldPosition(PosToMatch);
			if (wpos.first < mindist_3D * mindist_3D) {
				matchedPos.push_back(wpos.second);
				frame_indices.push_back(indices);
				raydists.push_back(wpos.first);
				deque<Position> ttmp;
				for (int i = 0; i < ncams; ++i) {
					ttmp.push_back(*((*tm)[i]));
				}
				PosTouse.push_back(ttmp);
			}
		}
	}
	// finally, prune the matched positions so that we only allow one per 2D point
	deque<unsigned int> bad;
	for (unsigned int i = 0; i < matchedPos.size(); ++i) {
		// have we already thrown out this point?
		if (find(bad.begin(), bad.end(), i) != bad.end()) {
			continue;
		}
		double min = raydists[i];
		for (unsigned int j = i + 1; j < matchedPos.size(); ++j) {
			bool flag = false;
			for (int k = 0; k < ncams; ++k) {
				flag |= frame_indices[i][k] == frame_indices[j][k];
			}
			if (flag) {
				// point j uses at least one of the same 2D particles
				if (min < raydists[j]) {
					// throw out j
					bad.push_back(j);
				} else {
					// throw out i
					min = raydists[j];
					bad.push_back(i);
					break;
				}
			}
		}
    }
    
	//cout << "framenumber_in_Stereomatch_bottom: " << framenumber << endl;
	deque<Position> goodPos;
    deque< deque<Position> > goodPosTouse;
    
    cout << "\tmatchedPos.size(): " << matchedPos.size() << endl;
    cout << "\tbad.size(): " << bad.size() << endl;
    
	for (unsigned int i = 0; i < matchedPos.size(); ++i) {
		if (i>=matchedPos.size())
			break;
		// have we already thrown out this point?
		if (find(bad.begin(), bad.end(), i) != bad.end()) {
			continue;
		}
//         cout << "\n\tMatched 3D pos (in mm):\t" << matchedPos[i] << endl;
// 		cout << "\n\tMatched 3D pos (in mm):\t" << matchedPos[i].X() << " " << matchedPos[i].Y() << " " << matchedPos[i].Z() << "\n\tIntersect error (in mm):\t" << raydists[i] << endl;
//         for (int kam = 0; kam < ncams; ++kam) {
//             cout << "\t\tCamera " << kam << " (used 2d coords): [" << cams[kam].Distort((PosTouse[i])[kam]).X() << ", " << cams[kam].Distort((PosTouse[i])[kam]).Y() << "]" << endl;
//         }
		double tmp = framenumber; 
		outfile.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = matchedPos[i].X();
		outfile.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = matchedPos[i].Y();
		outfile.write(reinterpret_cast<const char*>(&tmp), 8);
		tmp = matchedPos[i].Z();
		outfile.write(reinterpret_cast<const char*>(&tmp), 8);
		//cout << "\tInfo: " << wpos.second.Info() << endl;
		tmp = raydists[i];
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[0].Distort((PosTouse[i])[0]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[0].Distort((PosTouse[i])[0]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[0].Distort((PosTouse[i])[0]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[1].Distort((PosTouse[i])[1]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[1].Distort((PosTouse[i])[1]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[1].Distort((PosTouse[i])[1]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[2].Distort((PosTouse[i])[2]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[2].Distort((PosTouse[i])[2]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[2].Distort((PosTouse[i])[2]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[3].Distort((PosTouse[i])[3]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[3].Distort((PosTouse[i])[3]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[3].Distort((PosTouse[i])[3]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        
		goodPos.push_back(matchedPos[i]);
//         cout << "\tgoodPos 3D pos (in mm):\t" << goodPos[i].X() << " " << goodPos[i].Y() << " " << goodPos[i].Z() << "\n\tIntersect error (in mm):\t" << raydists[i] << endl;
	}
    
    for (int kam = 0; kam < ncams; ++kam) { 
        deque<Position> tttmp;
        for (unsigned int i = 0; i < matchedPos.size(); ++i) {
            if (i>=matchedPos.size())
                break;
            // have we already thrown out this point?
            if (find(bad.begin(), bad.end(), i) != bad.end()) {
                continue;
            }
            tttmp.push_back((PosTouse[i])[kam]);
        }
        goodPosTouse.push_back(tttmp);
    }

    cout << "\tgoodPos.size(): " << goodPos.size() << " (using all cameras)" << endl;
//     cout << "\tgoodPosTouse (size): " << " (" << goodPosTouse.size() << ") (using all cameras)" << endl;
//     
    // step 4.
	// go through the lists and search for consistency; that is, particles that show up 
	// on each other's lists. these are "real" matches, and we will compute their 3D positions.
	// note: we require the particles to be seen on all cameras -- so we only need to search along 
	// one of the lists. might as well pick cam 0!

    int pi_count = 0;
    
    vector < deque<Position> > goodPos3cams;
    vector < vector < deque<Position> > > goodPosTouse3cams;
    vector < deque<double> > raydists3cams;
    
    //vector<tempObject>().swap(tempVector); clears memory
    
    for (mcam = 0; mcam < ncams; ++mcam) {
        
        deque<Position> goodPos3;
        deque<Position> matchedPos3;
        vector < deque<Position> > PosTouse3;
        vector < deque<Position> > goodPosTouse3;
        deque< deque<int> > frame_indices3;
        deque<double> raydists3;
        
        for (int icam = 0; icam < ncams; ++icam) { // break at end of loop, so that only one icam different than mcam is picked
            
            if (icam == mcam) {
                continue;
            }
            
            int picam_count = 0;
            cout << "\tPerforming consistency checks, skipping camera " << mcam << endl;
            Frame::const_iterator picamend = corrframes[icam].end();
            for (Frame::const_iterator picam = corrframes[icam].begin(); picam != picamend; ++picam) {	
                // quick check 1: has this position already been used?
                bool toNextpicam = false;
                //cout << "\t\tParticle[" << picam.where() << "] picam_check 1: position: [" << cams[icam].Distort(*picam).X() << ", " << cams[icam].Distort(*picam).Y() << "]" << endl;
                for (unsigned int i = 0; i < goodPosTouse[icam].size(); ++i) {
                    if ((goodPosTouse[icam])[i].X() == (*picam).X() && (goodPosTouse[icam])[i].Y() == (*picam).Y()) {
                        //cout << "\t\t\t...has already been used: [" << cams[icam].Distort((PosTouse[i])[icam]).X() << ", " << cams[icam].Distort((PosTouse[i])[icam]).Y() << "]" << endl;
                        toNextpicam = true;
                        break;
                    }
                }
                if (toNextpicam) {
                    continue;
                }
                // quick check 2: does picam have any other particles on at least two of its lists?
                for (int i = 0; i < ncams; ++i) {
                    if (i == mcam || i == icam) {
                        continue;
                    }
                    if (pairlists[icam][picam.where()][i].size() == 0) {
                        toNextpicam = true;
                    }
                }
                if (toNextpicam) {
                    //cout << "\t\t\t...Particle[" << picam.where() << "] is missing on " << ncams_missing << " cameras" << endl;
                    continue;  // to next particle on camera icam
                }
                
//                 cout << "\tChecking particle [" << picam.where() << "] on camera: " << icam  << " (mcam=" << mcam <<")" << endl;
                
                // CONSIDER ONLY PARTICLES THAT ARE MISSING ON CAMERA SKIP_CAM ONLY
                
                vector< deque<Frame::const_iterator> > toMatch3;
                deque<Frame::const_iterator> tmp3;
                tmp3.push_back(picam);
                toMatch3.push_back(tmp3);

                int cam_count = 0;
                // loop over the other cameras
                for (int i = 0; i < ncams; ++i) {
                    if (i == mcam || i == icam) {
//                         cout << "\t\tskipping cam: " << i << endl;
                        continue;
                    }
                    cam_count += 1;
//                     cout << "\t\ttoMatch3.size(): " << toMatch3.size() << " cam_count: " << cam_count << endl;
                    unsigned int s = toMatch3.size();
                    for (unsigned int k = 0; k < s; ++k) {
                        // is this set too small already?
//                         cout << "\t\ttoMatch3[k=" << k << "].size(): " << toMatch3[k].size() << endl;
                        if (toMatch3[k].size() < static_cast<size_t>(cam_count)) {
//                             cout << "\t\tNot enough candidates in toMatch[k]" << endl;
                            continue;
                        }
                        // now, can we push any points onto the end of this set?
                        Frame::const_iterator piend = corrframes[i].end();
                        for (Frame::const_iterator pi = corrframes[i].begin(); pi != piend; ++pi) {
                            // does this point have anything on its lists?
                            // first check if this particle's 2d positions have already been used before!
                            bool toNextpi = false;
                            for (unsigned int ii = 0; ii < goodPosTouse[i].size(); ++ii) {
                                if ((goodPosTouse[i])[ii].X() == (*pi).X() && (goodPosTouse[i])[ii].Y() == (*pi).Y()) {
                                    toNextpi = true;
                                    break;
                                }
                            }
                            //cout << "\t\tChecking particle[" << p.where() << "] on camera " << i << endl;
                            for (int j = 0; j < ncams; ++j) {
                                if (j == i || j == mcam) {
                                    continue;
                                }
                                if (pairlists[i][pi.where()][j].size() == 0) {
                                    toNextpi = true;
                                    break;
                                }
                            }
                            if (toNextpi) {
//                                 cout << "\t\tParticle[" << pi.where() << "] on camera " << i << " is missing on more than one camera or its 2d position has already been used." << endl;
                                continue;
                            }
//                             cout << "\tPossible candidate particle [" << pi.where() << "] on camera: " << i  << " (mcam=" << mcam << ",icam=" << icam << ")" << endl;
                            // loop over the points already in the point set
                            bool canAdd = true;
                            // is this point picam in the point set on pi's list?
                            int tm_count = 0;
                            for (int j = 0; j < ncams; ++j) {
                                if (j == i || j == mcam) {
                                    continue;
                                }
                                list<Frame::const_iterator>::const_iterator onlist_tm_pi = find(pairlists[i][pi.where()][j].begin(),pairlists[i][pi.where()][j].end(),toMatch3[k][tm_count]); //(*tm)[j]);
                                if (onlist_tm_pi != pairlists[i][pi.where()][j].end()) {
                                    // yes! is the converse true?
//                                     cout << "\t\t\tParticle [" << toMatch3[k][tm_count].where() << "][" << cams[j].Distort(*(toMatch3[k][tm_count])).X() << ", " << cams[j].Distort(*(toMatch3[k][tm_count])).Y() << "] on cam " << j << " was found on pairlists["<< i <<"]["<< pi.where() <<"]["<< j <<"]" << endl;
//                                     cout << "\t\t\tChecking particle position: [" << cams[j].Distort(**(onlist_tm_pi)).X() << ", " << cams[j].Distort(**(onlist_tm_pi)).Y() << "]" << endl;
                                    list<Frame::const_iterator>::const_iterator onlist_pi_tm = find(pairlists[j][toMatch3[k][tm_count].where()][i].begin(),pairlists[j][toMatch3[k][tm_count].where()][i].end(),pi);
                                    if (onlist_pi_tm != pairlists[j][toMatch3[k][tm_count].where()][i].end()) {
                                        // yes! we can potentially add it
//                                         cout << "\t\t\tChecking the reverse: CAN ADD particle [" << pi.where() << "]" << endl;
                                        canAdd &= true;
                                    } else {
//                                         if ( s > 1) {
//                                             cout << "\t\t\tChecking the reverse failed..." << endl;
//                                         }
                                        // nope
                                        canAdd &= false;
                                    }
                                }
                                else {
                                    // nope
//                                     if ( s > 1) {
//                                         cout << "\t\t\tParticle [" << toMatch3[k][tm_count].where() << "][" << cams[j].Distort(*(toMatch3[k][tm_count])).X() << ", " << cams[j].Distort(*(toMatch3[k][tm_count])).Y() << "] on cam " << j << " was not found on pairlists["<< i <<"]["<< pi.where() <<"]["<< j <<"]" << endl;
//                                     }
                                    canAdd &= false;
                                }
                                tm_count += 1;
                                if (s == 1) {
                                    break;
                                }
                            }
                            if (canAdd) {
                                // we passed all the tests. add this point!
                                // note: add it as a new queue so that we can add other points too, if we want to.
                                deque<Frame::const_iterator> topush = toMatch3[k];//*tm;
                                topush.push_back(pi);
                                toMatch3.push_back(topush); 
//                                 cout << "\t\tpush_back pi [" << pi.where() << "]" << endl;
//                                 cout << "\t\ttoMatch3[k=" << k << "].size(): " << toMatch3[k].size() << endl;
                            }
                        }
                    }
                }
                // loop over two cameras other than icam and mcam
                // now, match the point sets that are full!
                vector< deque<Frame::const_iterator> >::const_iterator tmend = toMatch3.end();
                for (vector< deque<Frame::const_iterator> >::const_iterator tm = toMatch3.begin(); tm != tmend; ++tm) {
                    // is this set too small?
                    if (tm->size() < static_cast<size_t>(ncams-1)) {
                        // yes
                        continue;
                    }

                    // we can match it!
                    deque<Position> PosToMatch3;
                    deque<int> indices3;
                    int ic = 0;
                    for (int i = 0; i < ncams; ++i) {
                        if (i == mcam) {
    //                             cout << "*((*tm)[ic]): " << m_pos() << endl;
                            PosToMatch3.push_back(m_pos());
                            continue;
                        }
    //                         cout << "*((*tm)[ic]): " << *((*tm)[ic]) << endl;
                        PosToMatch3.push_back(*((*tm)[ic]));
                        indices3.push_back((*tm)[ic].where());
                        ic += 1;
                    }
                    ic = 0;
                    for (int i = 0; i < ncams; ++i) {
                        if (i == mcam) {
    //                             cout << "PosToMatch3 (2D): Cam " << cams[i].Distort(m_pos()).Ori() << " - [missing cam]" << endl;
                            continue;
                        }
    //                         cout << "PosToMatch3 (2D): Cam " << i << " - [" << cams[i].Distort(((*tm)[ic])[ic]).X() << ", " << cams[i].Distort(((*tm)[ic])[ic]).Y() << "][" << (*tm)[ic].where() << "]" << endl;
                        ic += 1;
                    }
                    picam_count += 1;
                    // for each camera (icam, or each iteration of a different missing cam mcam) match the pointset PosToMatch3 and if intersect error small enough add to matchedPos3
                    pair<double,Position> wpos = WorldPosition(PosToMatch3);
                    if (wpos.first < mindist_3D * mindist_3D) {
//                         cout << "\tmatchedPos3.size(): " << matchedPos3.size() << " push_back: [" << (wpos.second).X() << ", " << (wpos.second).Y() << ", " << (wpos.second).Z() << "]" << endl;
                        matchedPos3.push_back(wpos.second);
                        frame_indices3.push_back(indices3);
                        raydists3.push_back(wpos.first);
                        deque<Position> ttmp3;
                        ic = 0;
                        for (int i = 0; i < ncams; ++i) {
                            if (i == mcam) {
                                ttmp3.push_back(cams[i].UnDistort(m_pos()));
                                continue;
                            }
                            ttmp3.push_back(*((*tm)[ic]));
                            ic += 1;
                        }
                        PosTouse3.push_back(ttmp3);
                        pi_count += 1;
                    } else {
    //                         cout << "\tERROR: Particle positions were not matched! (" << wpos.first << " < " << mindist_3D * mindist_3D << ")" << endl;
                        picam_count -= 1;
                    }
                }
            }
            
//             cout << picam_count << " possible particles skipping camera: " << mcam << endl;
//             cout << pi_count << " particles have passed all tests from 3 camera configuration, skipping camera: " << mcam << endl;

            // finally, prune the matched positions again so that we only allow one per 2D point
            deque<unsigned int> bad3;
//             cout << "\tmatchedPos3.size(): " << matchedPos3.size() << " mcam: " << mcam << endl;
            for (unsigned int i = 0; i < matchedPos3.size(); ++i) {
                // have we already thrown out this point?
                if (find(bad3.begin(), bad3.end(), i) != bad3.end()) {
                    continue;
                }
                double min = raydists3[i];
                for (unsigned int j = i + 1; j < matchedPos3.size(); ++j) {
                    bool flag = false;
    //                 int kc = 0;
                    for (int k = 0; k < ncams-1; ++k) {
                        flag |= frame_indices3[i][k] == frame_indices3[j][k];
                    }
                    if (flag) {
//                         cout << "frame_indices3[i][k] appears twice, frame_indices3[j][k];" << endl;
                        // point j uses at least one of the same 2D particles
                        if (min < raydists3[j]) {
                            // throw out j
                            bad3.push_back(j);
                        } else {
                            // throw out i
                            min = raydists3[j];
                            bad3.push_back(i);
                            break;
                        }
                    }
                }
            }
        
            // check goodPos for 2d postions appearing multiple times
            for (unsigned int i = 0; i < matchedPos3.size(); ++i) {
                // don't check particles that have already been deleted
                if (find(bad3.begin(), bad3.end(), i) != bad3.end()) {
                    continue;
                }
                for (unsigned int j = 0; j < matchedPos.size(); ++j) {
                    // compare only with good Postouse from 4 cams
                    if (find(bad.begin(), bad.end(), j) != bad.end()) {
                        continue;
                    }
                    bool flagX = false;
                    bool flagY = false;
                    for (int k = 0; k < ncams; k++) {
                        flagX |= (PosTouse3[i])[k].X() == (PosTouse[j])[k].X();
                        flagY |= (PosTouse3[i])[k].Y() == (PosTouse[j])[k].Y();
                    }
                    if (flagX || flagY) {
//                         cout << "PosTouse3[i][k] appears twice, PosTouse[j][k];" << endl;
                        bad3.push_back(i);
                    }
                }
            }

//             cout << "matchedPos3.size(): " << matchedPos3.size() << endl;
//             cout << "bad3.size(): " << bad3.size() << endl;
//             cout << "raydists3.size(): " << raydists3.size() << endl;

            for (unsigned int i = 0; i < matchedPos3.size(); ++i) {
                if (i>=matchedPos3.size()) {
                    break;
                }
                // have we already thrown out this point?
                if (find(bad3.begin(), bad3.end(), i) != bad3.end()) {
                    continue;
                }
                goodPos3.push_back(matchedPos3[i]);
            }
            for (int kam = 0; kam < ncams; ++kam) { 
                deque<Position> tttmp3;
                for (unsigned int i = 0; i < matchedPos3.size(); ++i) {
                    if (i>=matchedPos3.size())
                        break;
                    // have we already thrown out this point?
                    if (find(bad3.begin(), bad3.end(), i) != bad3.end()) {
                        continue;
                    }
                    tttmp3.push_back((PosTouse3[i])[kam]);
                }
                goodPosTouse3.push_back(tttmp3);
            }
//             cout << "goodPos3 (size): " << goodPos3.size() << " (3 cameras)" << endl;
//             cout << "goodPosTouse3 (size): " << goodPosTouse3.size() << " (3 cameras)" << endl;
            break;
        }
        
        goodPos3cams.push_back(goodPos3);
        goodPosTouse3cams.push_back(goodPosTouse3);
        raydists3cams.push_back(raydists3);
        
    }
        
    deque<Position> goodPos3;
    vector < deque<Position> > goodPosTouse3;
    deque<double> raydists3;
    
    for (int mcam = 0; mcam < ncams; ++mcam) {
        for (unsigned int i = 0; i < goodPos3cams[mcam].size(); ++i) {
            goodPos3.push_back(goodPos3cams[mcam][i]);
            raydists3.push_back(raydists3cams[mcam][i]);
        }
    }
//     cout << "goodPos3.size(): " << goodPos3.size() << endl;
//     cout << "raydists3.size(): " << raydists3.size() << endl;
    
    for (int icam = 0; icam < ncams; ++icam) {
        deque<Position> temp_pos;
        for (int mcam = 0; mcam < ncams; ++mcam) {
            vector < deque<Position> > temp = goodPosTouse3cams[mcam];
//             cout << "temp.size(): " << temp[icam].size() << endl;
            for (unsigned int tt = 0; tt < temp[icam].size(); ++tt) {
                temp_pos.push_back(temp[icam][tt]);
            }
        }
        goodPosTouse3.push_back(temp_pos);
    }
    
//     cout << "goodPosTouse3[0].size(): " << goodPosTouse3[0].size() << endl;
//     cout << "goodPosTouse3[1].size(): " << goodPosTouse3[1].size() << endl;
//     cout << "goodPosTouse3[2].size(): " << goodPosTouse3[2].size() << endl;
//     cout << "goodPosTouse3[3].size(): " << goodPosTouse3[3].size() << endl;
    
    deque<unsigned int> bad3cams;
    for (unsigned int i = 0; i < goodPos3.size(); ++i) {
        double min = raydists3[i];
//         cout << "raydists3[" << i << "]:" << raydists3[i] << endl;
        for (unsigned int j = i+1; j < goodPos3.size(); ++j) {
            if (find(bad3cams.begin(), bad3cams.end(), j) != bad3cams.end()) {
                continue;
            }
            
            bool flagX = false;
            bool flagY = false;
            for (int k = 0; k < ncams; k++) {
                // BUT ONLY IF MISSING CAM IS NEGLETED!!!!!
                if ((goodPosTouse3[k])[i].Ori()==k) {
                    continue;
                }
                flagX |= (goodPosTouse3[k])[i].X() == (goodPosTouse3[k])[j].X();
                flagY |= (goodPosTouse3[k])[i].Y() == (goodPosTouse3[k])[j].Y();
            }
            if (flagX || flagY) {
//                 cout << "goodPosTouse3[i][k] appears twice, goodPosTouse3[j][k];" << endl;
                if (min < raydists3[j]) {
                    // throw out j
                    bad3cams.push_back(j);
                } else {
                    // throw out i
                    min = raydists3[j];
                    bad3cams.push_back(i);
                    break;
                }
            }
        }
    }
    
//     Check goodPosTouse from 3 cameras and eliminate positions that have been used twice
//     Check between the combinations of three cameras if a position has been used multiple times
//     cout << "goodPos3.size(): " << goodPos3.size() << endl;
//     cout << "bad3.size(): " << bad3.size() << endl;
//     cout << "raydists3.size(): " << raydists3.size() << endl;

    for (unsigned int i = 0; i < goodPos3.size(); ++i) {
        if (i>=goodPos3.size()) {
            break;
        }
        // have we already thrown out this point?
        if (find(bad3cams.begin(), bad3cams.end(), i) != bad3cams.end()) {
            continue;
        }
//         cout << "\n\tMatched 3D pos (in mm):\t" << matchedPos3[i] << endl;
//         cout << "\n\tgoodPos3 3D pos (in mm):\t" << goodPos3[i].X() << " " << goodPos3[i].Y() << " " << goodPos3[i].Z() << "\n\tIntersect error (in mm):\t" << raydists3[i] << endl;
//         for (int kam = 0; kam < ncams; ++kam) {
//             cout << "\t\tCamera " << kam << " (used 2d coords): [" << cams[kam].Distort((goodPosTouse3[kam])[i]).X() << ", " << cams[kam].Distort((goodPosTouse3[kam])[i]).Y() << "]" << endl;
//         }
        double tmp = framenumber; 
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = goodPos3[i].X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = goodPos3[i].Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = goodPos3[i].Z();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        //cout << "\tInfo: " << wpos.second.Info() << endl;
        tmp = raydists3[i];
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[0].Distort((goodPosTouse3[0])[i]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[0].Distort((goodPosTouse3[0])[i]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[0].Distort((goodPosTouse3[0])[i]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);  
        tmp = cams[1].Distort((goodPosTouse3[1])[i]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[1].Distort((goodPosTouse3[1])[i]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[1].Distort((goodPosTouse3[1])[i]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[2].Distort((goodPosTouse3[2])[i]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[2].Distort((goodPosTouse3[2])[i]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[2].Distort((goodPosTouse3[2])[i]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[3].Distort((goodPosTouse3[3])[i]).X();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[3].Distort((goodPosTouse3[3])[i]).Y();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        tmp = cams[3].Distort((goodPosTouse3[3])[i]).Ori();
        outfile.write(reinterpret_cast<const char*>(&tmp), 8);
        goodPos.push_back(goodPos3[i]);
    }
    
    // now free memory
	for (int i = 0 ; i < ncams; ++i) {	
		for (int j = 0; j < corrframes[i].NumParticles(); ++j) {
			delete []pairlists[i][j];
		}
		delete []pairlists[i];
	}
	delete pairlists;
    
    cout << "\tgoodPos.size(): " << goodPos.size() << endl;
    return Frame(goodPos);	
}


pair<double,Position> Calibration::WorldPosition(deque<Position> ipos) throw(runtime_error)
{
    int ncams_missing = 0;
//     cout << "\tMissing cam: " << mcam << endl;
    if (mcam == -1 && ipos.size() != cams.size()) {
        throw runtime_error("Number of cameras and number of images do not match!");
    }
    if (mcam != -1) {
//         cout << "\tipos.size(): " << ipos.size() << endl;
//         for (unsigned int k = 0; k < ipos.size(); ++k) {
//             cout << "ipos[i=" << k << "]: " << ipos[k] << endl;
//         }
        ncams_missing = 1;
    }
	// solve, in a least-squares sense, the intersection of the three lines of sight.
	// this will give us the distance of closest approach to the three lines.
	// embarrassingly, see 
	// http://en.wikipedia.org/wiki/Line-line_intersection

	Matrix M;
	Position P(0, 0, 0);
	
	Position sight[ncams-ncams_missing];
    
	int ic = 0;
	for (int i = 0; i < ncams; ++i) {
        if (i == mcam) {
            continue;
        }
		// construct a line of sight for this point:
		// vector from camera center to this point (in world coordinates)
		sight[ic] = cams[i].ImageToWorld(ipos[i]) - cams[i].Center();
		// normalize
		sight[ic] /= sight[ic].Magnitude();
		
		// add to the least-squares matrices
		
		// inelegant, but it works
		Matrix tmp;
		tmp.Set(0, 0, 1 - sight[ic].X() * sight[ic].X());
		tmp.Set(0, 1, -sight[ic].X() * sight[ic].Y());
		tmp.Set(0, 2, -sight[ic].X() * sight[ic].Z());
		tmp.Set(1, 0, -sight[ic].Y() * sight[ic].X());
		tmp.Set(1, 1, 1 - sight[ic].Y() * sight[ic].Y());
		tmp.Set(1, 2, -sight[ic].Y() * sight[ic].Z());
		tmp.Set(2, 0, -sight[ic].Z() * sight[ic].X());
		tmp.Set(2, 1, -sight[ic].Z() * sight[ic].Y());
		tmp.Set(2, 2, 1 - sight[ic].Z() * sight[ic].Z());
		
		P += tmp * cams[i].Center();
		M += tmp;
        ic += 1;
	}
	
	// invert the matrix and construct the 3D position
	Position tmpi(M.Invert() * P);
	Position worldpos(M.Invert() * P);
	
	// calculate the rms distance from worldpos to each ray
	double dist = 0;
    ic = 0;
	for (int i = 0; i < ncams; ++i) {
        if (i == mcam) {
            continue;
        }
		Position h = (worldpos - Dot(worldpos, sight[ic]) * sight[ic] - (cams[i].Center() - Dot(cams[i].Center(), sight[ic]) * sight[ic]));
		dist += h.Magnitude2();
        ic += 1;
	}
	dist /= static_cast<double>(ncams-ncams_missing);

	Position worldposi(tmpi.X(),tmpi.Y(),tmpi.Z(), cams[0].Distort(ipos[0]).X(), cams[0].Distort(ipos[0]).Y(), ipos[0].Ori(),cams[1].Distort(ipos[1]).X(), cams[1].Distort(ipos[1]).Y(), ipos[1].Ori(),cams[2].Distort(ipos[2]).X(), cams[2].Distort(ipos[2]).Y(),  ipos[2].Ori(), cams[3].Distort(ipos[3]).X(), cams[3].Distort(ipos[3]).Y(), ipos[3].Ori(), dist);
	return make_pair(dist,worldposi);
}
