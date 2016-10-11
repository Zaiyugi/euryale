//*******************************************************************
//
//   Vector.h
//
// 3D vector class in the namespace lux
//
//
//
//*******************************************************************

#ifndef __LUX_VECTOR_H__
#define __LUX_VECTOR_H__

#include <iostream>
#include <fstream>
#include <cmath>

namespace lux
{

//! Vector is a 3D vector class
class Vector
{
  public:

   Vector(){ xyz[0] = xyz[1] = xyz[2] = 0; }

   Vector(const double v){ xyz[0] = xyz[1] = xyz[2] = v; }

   Vector(const Vector& v)
   { 
      xyz[0] = v.xyz[0];
      xyz[1] = v.xyz[1];
      xyz[2] = v.xyz[2]; 
   }
   
   Vector(const double a, const double b, const double c)
   {
      xyz[0] = a;
      xyz[1] = b;
      xyz[2] = c; 
   }

   Vector(const double v[3])
   {
      xyz[0] = v[0];
      xyz[1] = v[1];
      xyz[2] = v[2]; 
   }

   ~Vector(){}

   //!  Set all three components
   void set( const float vx, const float vy, const float vz )
   {
      xyz[0] = vx;
      xyz[1] = vy;
      xyz[2] = vz;
   }

   //! Add two vectors together
   const Vector operator+        (const Vector& v) const 
   { 
      return Vector(xyz[0]+v.xyz[0], xyz[1]+v.xyz[1], xyz[2]+v.xyz[2]); 
   }
  
   //! Subtract one vector from another
   const Vector operator-        (const Vector& v) const
   { 
      return Vector(xyz[0]-v.xyz[0], xyz[1]-v.xyz[1], xyz[2]-v.xyz[2]); 
   }

   //! Unary minus
   friend const Vector operator- (const Vector& v)
   { return Vector(-v.xyz[0],-v.xyz[1],-v.xyz[2]); }

   //! Multiplication of a constant with a vector
   friend const Vector operator* (const double w, const Vector& v)
   { return v*w; }
	  
   //! Multiplication of a vector with a constant
   const Vector operator*        (const double v) const
   { return Vector(xyz[0]*v, xyz[1]*v, xyz[2]*v); }

   const Vector operator/        (const double v) const
   { return Vector(xyz[0]/v, xyz[1]/v, xyz[2]/v); }

   // Component-wise divide
   const Vector operator/        (const Vector v) const
   { return Vector(xyz[0]/v[0], xyz[1]/v[1], xyz[2]/v[2]); }

   // Component-wise multiply
   const Vector compMult        (const Vector v) const
   { return Vector(xyz[0]*v[0], xyz[1]*v[1], xyz[2]*v[2]); }

   //! Inner product
   const double operator*        (const Vector& v) const  
   { return (xyz[0]*v.xyz[0] + xyz[1]*v.xyz[1] + xyz[2]*v.xyz[2]); }
  
   //! cross product
   const Vector operator^        (const Vector& v) const 
   { return Vector(xyz[1]*v.xyz[2] - xyz[2]*v.xyz[1], 
		   xyz[2]*v.xyz[0] - xyz[0]*v.xyz[2], 
		   xyz[0]*v.xyz[1] - xyz[1]*v.xyz[0]); }

   Vector& operator=       (const Vector& v)
   { xyz[0] = v.xyz[0]; xyz[1] = v.xyz[1]; xyz[2] = v.xyz[2]; return *this; }
  
   Vector& operator+=      (const Vector& v)
   { xyz[0] += v.xyz[0]; xyz[1] += v.xyz[1]; xyz[2] += v.xyz[2]; return *this; }
  
   Vector& operator-=      (const Vector& v)
   { xyz[0] -= v.xyz[0]; xyz[1] -= v.xyz[1]; xyz[2] -= v.xyz[2]; return *this; }
  
   Vector& operator*=      (const double v)
   { xyz[0] *= v; xyz[1] *= v; xyz[2] *= v; return *this; }
  
   Vector& operator/=      (const double v)
   { xyz[0] /= v; xyz[1] /= v; xyz[2] /= v; return *this; }
  

   const double& operator[] (const int v) const { return xyz[v]; }
         double& operator[] (const int v)       { return xyz[v]; }
   const double& operator() (const int v) const { return xyz[v]; }

   const double X() const { return xyz[0]; }
   const double Y() const { return xyz[1]; }
   const double Z() const { return xyz[2]; }

   const double x() const { return xyz[0]; }
   const double y() const { return xyz[1]; }
   const double z() const { return xyz[2]; }

   const double magnitude() const 
   { return sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] ); }
   
   const double sqr_mag() const
   { return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]; }
   
   const Vector unit() const { return *this/magnitude(); }

   const Vector unitvector() const
   {
      if(magnitude() < 5e-8)
         return *this;

      return *this/magnitude();
   }

   const Vector v_abs() const { return Vector(fabs(xyz[0]), fabs(xyz[1]), fabs(xyz[2])); }

   const Vector cross(const Vector& v) const { return *this^v; }
   const float dot(const Vector& v) const { return *this*v; }

   void normalize() 
   {
		double mag = magnitude();
		if(mag < 5e-8) return;

		xyz[0] /= mag;
		xyz[1] /= mag;
		xyz[2] /= mag;
	}

   bool isNan()
   {
      if(std::isnan(xyz[0]) || std::isnan(xyz[1]) || std::isnan(xyz[2]))
         return true;

      return false;
   }

//  Comparisons

   const bool operator==         (const Vector& v) const
       { return ( xyz[0]==v.xyz[0] && xyz[1]==v.xyz[1] && xyz[2]==v.xyz[2] ); }
  
   const bool operator!=         (const Vector& v) const
       { return ( xyz[0]!=v.xyz[0] || xyz[1]!=v.xyz[1] || xyz[2]!=v.xyz[2] ); }
  
   const bool operator<          (const Vector& v) const
       { return ( xyz[0]<v.xyz[0] && xyz[1]<v.xyz[1] && xyz[2]<v.xyz[2] ); }
  
   const bool operator<=         (const Vector& v) const
       { return ( xyz[0]<=v.xyz[0] && xyz[1]<=v.xyz[1] && xyz[2]<=v.xyz[2] ); }
  
   const bool operator>          (const Vector& v) const
       { return ( xyz[0]>v.xyz[0] && xyz[1]>v.xyz[1] && xyz[2]>v.xyz[2] ); }
  
   const bool operator>=         (const Vector& v) const
       { return ( xyz[0]>=v.xyz[0] && xyz[1]>=v.xyz[1] && xyz[2]>=v.xyz[2] ); }

   // Test for parallel
   const bool operator||         (const Vector& v) const
       { return (  fabs((*this)*v) == v.magnitude()*((*this).magnitude()) ); }

   void read(std::ifstream& ifs) { ifs >> xyz[0] >> xyz[1] >> xyz[2]; }
   void write(std::ofstream& ofs) { ofs << xyz[0] << " " << xyz[1] << " " << xyz[2]; }
 
   friend std::ostream& operator<<(std::ostream &out, const Vector& v)
   {
      out << v[0] << " " << v[1] << " " << v[2];
      return out;
   }

  private:
  double xyz[3];
};

const static lux::Vector EUCL3D_X_AXIS(1, 0, 0);
const static lux::Vector EUCL3D_Y_AXIS(0, 1, 0);
const static lux::Vector EUCL3D_Z_AXIS(0, 0, 1);

}


#endif
