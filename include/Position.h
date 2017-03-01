/*
 *  Position.h
 *
 *  A Position object holds the coordinates of a single particle, and 
 *  has lots of useful operators defined.
 *
 *  Last update: 10/28/11 by NTO (made 3D)
 *
 *  In collaboration with Wesleyan Universiy.
 *  All parts of these codes have been heavily modified by Stefan Kramel.
 *  Added features: - variable number of cameras for which a particle can be missing
 *                  - data format read in. from .avi files to .cpv and .gdf files
 *                  - write out of intermediate stereomatched data
 *
 *
 */

#ifndef POSITION_H
#define POSITION_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class Position {
public:
  // constructor: does nothing
  Position() {};
  // constructor: initializes the coordinates
  Position(double NewX, double NewY, double NewZ);
  Position(double NewX, double NewY, double NewZ, double NewOri);
  Position(double NewX, double NewY, double NewZ,
			double Newx1, double Newy1, double NewOri1, 
			double Newx2, double Newy2, double NewOri2,
			double Newx3, double Newy3, double NewOri3,
			double Newx4, double Newy4, double NewOri4,
			double Newinfo);
  // copy-constructor
  Position(const Position& p);
  // destructor: nothing to do
  ~Position() {};
  
  // set this position to be fake, i.e., an estimated position
  void SetFake();
  // check if this position is fake
  bool IsFake() const;

  // get x component
  double X() const;
  // get y component
  double Y() const;
	// get z component
  double Z() const;
  
  double X1() const;
  double Y1() const;
  double Ori1() const;
  double X2() const;
  double Y2() const;
  double Ori2() const;
  double X3() const;
  double Y3() const;
  double Ori3() const;
  double X4() const;
  double Y4() const;
  double Ori4() const;
  double Ori() const;
  double Info() const;
  
  // get the squared Euclidean distance between two Positions
  // (don't take the square root for efficiency)
  friend double Distance(const Position& p1, const Position& p2);
  // get the magnitude of the vector
  double Magnitude() const;
	// get the squared magnitude
	double Magnitude2() const;

	// scalar product
  friend double Dot(const Position& left, const Position& right);
	// element-wise multiplication
	friend const Position Multiply(const Position& left, const Position& right);
	
  // member operators
  Position& operator=(const Position& p);
  // vector sum and difference
  Position& operator+=(const Position& right);
  Position& operator-=(const Position& right);
  // scalar multiplication and division
  Position& operator*=(double right);
  Position& operator/=(double right);
	// element by element multiplication
	Position& operator*=(const Position& right);

  // non-member operators

  // comparison
  friend int operator==(const Position& left, const Position& right);
  friend int operator!=(const Position& left, const Position& right);
	// compare y values (for sorting)
	friend int operator<(const Position& left, const Position& right);
	friend int operator>(const Position& left, const Position& right);
  // vector sum and difference
  friend const Position operator+(const Position& left, const Position& right);
  friend const Position operator-(const Position& left, const Position& right);
  // scalar multiplication and division
  friend const Position operator*(const Position& left, double right);
  friend const Position operator*(double left, const Position& right);
  friend const Position operator/(const Position& left, double right);
  // printing
  friend std::ostream& operator<<(std::ostream& os, const Position& p);

	
	
private:
  double x;
  double y;
  double z;
  double ori;
  double x1;
  double y1;
  double ori1;
  double x2;
  double y2;
  double ori2;
  double x3;
  double y3;
  double ori3;
  double x4;
  double y4;
  double ori4;
  double info;

  // a flag specifying whether this position is real or estimated
  bool fake;    
};

//  Inline Function Definitions
inline Position::Position(double NewX, double NewY, double NewZ) : x(NewX), y(NewY), z(NewZ), fake(false)
{}

inline Position::Position(double NewX, double NewY, double NewZ, double NewOri) 
: x(NewX), y(NewY), z(NewZ), ori(NewOri), fake(false)
{}

inline Position::Position(
	double NewX, double NewY, double NewZ,
	double Newx1, double Newy1, double NewOri1, 
	double Newx2, double Newy2, double NewOri2,
	double Newx3, double Newy3, double NewOri3,
	double Newx4, double Newy4, double NewOri4,
	double Newinfo) 
: x(NewX), y(NewY), z(NewZ),
	x1(Newx1), y1(Newy1), ori1(NewOri1),
	x2(Newx2), y2(Newy2), ori2(NewOri2),
	x3(Newx3), y3(Newy3), ori3(NewOri3),
	x4(Newx4), y4(Newy4), ori4(NewOri4),
	info(Newinfo), fake(false)
{}

inline Position::Position(const Position& p) 
: x(p.x), y(p.y), z(p.z), ori(p.ori),
	x1(p.x1), y1(p.y1), ori1(p.ori1),
	x2(p.x2), y2(p.y2), ori2(p.ori2),
	x3(p.x3), y3(p.y3), ori3(p.ori3),
	x4(p.x4), y4(p.y4), ori4(p.ori4),
	info(p.info), fake(p.fake)
{}

inline void Position::SetFake()
{
  fake = true;
}

inline bool Position::IsFake() const
{
  return fake;
}

inline double Position::Magnitude() const
{
  return sqrt(x * x + y * y + z * z);
}

inline double Position::Magnitude2() const
{
  return (x * x + y * y + z * z);
}

inline double Position::X() const
{
	return x;
}

inline double Position::Y() const
{
	return y;
}

inline double Position::Z() const
{
	return z;
}


inline double Position::X1() const
{
	return x1;
}

inline double Position::Y1() const
{
	return y1;
}


inline double Position::X2() const
{
	return x2;
}

inline double Position::Y2() const
{
	return y2;
}


inline double Position::X3() const
{
	return x3;
}

inline double Position::Y3() const
{
	return y3;
}


inline double Position::X4() const
{
	return x4;
}

inline double Position::Y4() const
{
	return y4;
}

inline double Position::Ori1() const
{
	return ori1;
}
inline double Position::Ori2() const
{
	return ori2;
}
inline double Position::Ori3() const
{
	return ori3;
}
inline double Position::Ori4() const
{
	return ori4;
}

inline double Position::Ori() const
{
	return ori;
}

inline double Position::Info() const
{
	return info;
}

inline Position& Position::operator=(const Position& p)
{
  x = p.x;
  y = p.y;
	z = p.z;
  fake = p.fake;

  return *this;
}

inline Position& Position::operator+=(const Position& right)
{
  x += right.x;
  y += right.y;
	z += right.z;

  return *this;
}

inline Position& Position::operator-=(const Position& right)
{
  x -= right.x;
  y -= right.y;
	z -= right.z;

  return *this;
}

inline Position& Position::operator*=(double right)
{
  x *= right;
  y *= right;
	z *= right;

  return *this;
}

inline Position& Position::operator/=(double right)
{
  x /= right;
  y /= right;
	z /= right;

  return *this;
}

inline Position& Position::operator*=(const Position& right)
{
	x *= right.x;
	y *= right.y;
	z *= right.z;
	
	return *this;
}


#endif // POSITION_H
