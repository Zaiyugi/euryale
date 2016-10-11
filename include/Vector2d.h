//*******************************************************************
//
//   Vector2d.h
//
// 2D vector class in the namespace lux
//
//
//
//*******************************************************************

#ifndef __LUX_VECTOR2D_H__
#define __LUX_VECTOR2D_H__

#include <iostream>
#include <fstream>
#include <cmath>

namespace lux
{

//! Vector is a 2D vector class
class Vector2d
{
  public:

   Vector2d(){ xyz[0] = xyz[1] = 0; }

   Vector2d(const double v){ xyz[0] = xyz[1] = v; }

   Vector2d(const Vector2d& v)
   { 
      xyz[0] = v.xyz[0];
      xyz[1] = v.xyz[1];
   }
   
   Vector2d(const double a, const double b)
   {
      xyz[0] = a;
      xyz[1] = b;
   }

   Vector2d(const double v[2])
   {
      xyz[0] = v[0];
      xyz[1] = v[1];
   }

   ~Vector2d(){}

   //!  Set all three components
   void set( const float vx, const float vy )
   {
      xyz[0] = vx;
      xyz[1] = vy;
   }

   //! Add two vectors together
   const Vector2d operator+        (const Vector2d& v) const 
   { 
      return Vector2d(xyz[0]+v.xyz[0], xyz[1]+v.xyz[1]); 
   }
  
   //! Subtract one vector from another
   const Vector2d operator-        (const Vector2d& v) const
   { 
      return Vector2d(xyz[0]-v.xyz[0], xyz[1]-v.xyz[1]); 
   }

   //! Unary minus
   friend const Vector2d operator- (const Vector2d& v)
   { return Vector2d(-v.xyz[0],-v.xyz[1]); }

   //! Multiplication of a constant with a vector
   friend const Vector2d operator* (const double w, const Vector2d& v)
   { return v*w; }
	  
   //! Multiplication of a vector with a constant
   const Vector2d operator*        (const double v) const
   { return Vector2d(xyz[0]*v, xyz[1]*v); }

   const Vector2d operator/        (const double v) const
   { return Vector2d(xyz[0]/v, xyz[1]/v); }

   // Component-wise divide
   const Vector2d operator/        (const Vector2d v) const
   { return Vector2d(xyz[0]/v[0], xyz[1]/v[1]); }

   // Component-wise multiply
   const Vector2d compMult        (const Vector2d v) const
   { return Vector2d(xyz[0]*v[0], xyz[1]*v[1]); }

   //! Inner product
   const double operator*        (const Vector2d& v) const  
   { return (xyz[0]*v.xyz[0] + xyz[1]*v.xyz[1]); }
  
   Vector2d& operator=       (const Vector2d& v)
   { xyz[0] = v.xyz[0]; xyz[1] = v.xyz[1]; return *this; }
  
   Vector2d& operator+=      (const Vector2d& v)
   { xyz[0] += v.xyz[0]; xyz[1] += v.xyz[1]; return *this; }
  
   Vector2d& operator-=      (const Vector2d& v)
   { xyz[0] -= v.xyz[0]; xyz[1] -= v.xyz[1]; return *this; }
  
   Vector2d& operator*=      (const double v)
   { xyz[0] *= v; xyz[1] *= v; return *this; }
  
   Vector2d& operator/=      (const double v)
   { xyz[0] /= v; xyz[1] /= v; return *this; }
  

   const double& operator[] (const int v) const { return xyz[v]; }
         double& operator[] (const int v)       { return xyz[v]; }
   const double& operator() (const int v) const { return xyz[v]; }

   const double X() const { return xyz[0]; }
   const double Y() const { return xyz[1]; }

   const double x() const { return xyz[0]; }
   const double y() const { return xyz[1]; }

   const double magnitude() const 
   { return sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] ); }
   
   const double sqr_mag() const
   { return xyz[0]*xyz[0] + xyz[1]*xyz[1]; }
   
   const Vector2d unit() const { return *this/magnitude(); }

   const Vector2d unitvector() const
   {
      if(magnitude() < 5e-8)
         return Vector2d();

      return *this/magnitude();
   }

   const Vector2d v_abs() const { return Vector2d(fabs(xyz[0]), fabs(xyz[1])); }

   const float dot(const Vector2d& v) const { return *this*v; }

   void normalize() 
   {
		double mag = magnitude();
		if(mag < 5e-8) return;

		xyz[0] /= mag;
		xyz[1] /= mag;
	}

   bool isNan()
   {
      if(std::isnan(xyz[0]) || std::isnan(xyz[1]))
         return true;

      return false;
   }

   bool isInf()
   {
      if(std::isinf(xyz[0]) || std::isinf(xyz[1]))
         return true;

      return false;
   }

//  Comparisons

   const bool operator==         (const Vector2d& v) const
       { return ( xyz[0]==v.xyz[0] && xyz[1]==v.xyz[1] ); }
  
   const bool operator!=         (const Vector2d& v) const
       { return ( xyz[0]!=v.xyz[0] || xyz[1]!=v.xyz[1] ); }
  
   const bool operator<          (const Vector2d& v) const
       { return ( xyz[0]<v.xyz[0] && xyz[1]<v.xyz[1] ); }
  
   const bool operator<=         (const Vector2d& v) const
       { return ( xyz[0]<=v.xyz[0] && xyz[1]<=v.xyz[1] ); }
  
   const bool operator>          (const Vector2d& v) const
       { return ( xyz[0]>v.xyz[0] && xyz[1]>v.xyz[1] ); }
  
   const bool operator>=         (const Vector2d& v) const
       { return ( xyz[0]>=v.xyz[0] && xyz[1]>=v.xyz[1] ); }

   // Test for parallel
   const bool operator||         (const Vector2d& v) const
       { return (  fabs((*this)*v) == v.magnitude()*((*this).magnitude()) ); }

   void read(std::ifstream& ifs) { ifs >> xyz[0] >> xyz[1]; }
   void write(std::ofstream& ofs) { ofs << xyz[0] << " " << xyz[1]; }
 
   friend std::ostream& operator<<(std::ostream &out, const Vector2d& v)
   {
      out << v[0] << " " << v[1];
      return out;
   }

  private:
  double xyz[2];
};

const static lux::Vector2d EUCL2D_X_AXIS(1, 0);
const static lux::Vector2d EUCL2D_Y_AXIS(0, 1);

}


#endif
