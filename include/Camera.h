/*
 *  Camera.h
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
 */

#ifndef CAMERA_H
#define CAMERA_H

#include <fstream>

#include <Position.h>
#include <Matrix.h>
#include <string.h>
#include <limits.h>

class Camera {
public:
	Camera(std::istream& is);
	Camera(const Camera& c);
	~Camera() {};
	
	// return camera projective center (in world coordinates)
	Position Center() const;
	
	// remove distortion; return centered coordinates in physical units
	Position UnDistort(const Position& p) const;
	// add distortion back; return normal images coordinates in pixel units
	Position Distort(const Position& p) const;
	
	// project a position (distorted, in pixels!) on the image plane to 3D world coordinates (in mm)
	Position ImageToWorld(const Position& p) const;
	// project a 3D world position (in mm) to a position on the image plane (undistorted, in mm)
	Position WorldToImage(const Position& p) const;

private:
	// model parameters; names should be more or less the same as from calibTsai.m!
	int Npixw;
	int Npixh;
	double wpix;
	double hpix;
	double f_eff;
	double kr;
	double kx;
	Matrix R;
	Position T;
	Matrix Rinv;
	Position Tinv;
		
};

inline Camera::Camera(const Camera& c)
: Npixw(c.Npixw), Npixh(c.Npixh), wpix(c.wpix), hpix(c.hpix), f_eff(c.f_eff),
kr(c.kr), kx(c.kx), R(c.R), T(c.T), Rinv(c.Rinv), Tinv(c.Tinv)
{}

inline Position Camera::Center() const
{
	return Tinv;
}

inline Position Camera::UnDistort(const Position& p) const
{
	//Position centered(p);
	// shift origin to the center of the image
	//centered -= Position(Npixw/2, Npixh/2, 0);
	
	// account for left-handed coordinate system...
	Position centered(p.X() - Npixw/2, -p.Y() + Npixh/2, p.Z(), p.Ori());
	
	// scale into physical units
	centered *= Position(wpix, hpix, 1, 1);
	// remove possible cylindrical distortion
	// 	cout << corrframes[0] << endl;should this be 1.0/kx? i.e., should kx be bigger or smaller than 1?
// 	centered *= Position(kx, 1, 1);
	// compute the radial correction factor
// 	double rad = 1.0 + kr * centered.Magnitude2();
	// apply the radial correction
// 	centered *= rad;
	// finally, return the undistorted coordinates, but keep them centered and in mm
	// std::cout << centered << std::endl;
	return Position(centered.X(),centered.Y(),centered.Z(),p.Ori());
}

inline Position Camera::Distort(const Position& p) const
{
	// compute the radial distortion factor
// 	double rad = 1.0 + kr * p.Magnitude2();
// 	Position pixelcoords(p / rad);
    Position pixelcoords(p);
	// remove potential cylindrical distortion
// 	pixelcoords *= Position(1.0 / kx, 1, 0, 0);
	// scale into pixel units
	pixelcoords *= Position(1.0 / wpix, 1.0 / hpix, 1, 1);
	// shift origin
	return Position(pixelcoords.X() + Npixw/2, -1.0 * (pixelcoords.Y() - Npixh/2), p.Z(), p.Ori());
}

inline Position Camera::ImageToWorld(const Position& p) const
{
	Position pp(p.X(),p.Y(),p.Z());
	Position tmp(pp * T.Z() / f_eff);
	Position proj(tmp.X(), tmp.Y(), T.Z());
	// return coordinates of p in 3D world coordinates
	//return ((proj - T) * R);
	Position tmpi(Rinv * (proj - T));
	return Position(tmpi.X(),tmpi.Y(),tmpi.Z(),p.Ori());
}

inline Position Camera::WorldToImage(const Position& p) const
{
	//Position proj(p * Rinv + T);
	Position pp(p.X(),p.Y(),p.Z());
	Position proj(R * pp + T);
	Position tmpi(proj * (f_eff / proj.Z()));
	return Position(tmpi.X(),tmpi.Y(),tmpi.Z(),p.Ori());
}

#endif // CAMERA_H
