/* Name: Zachary Shore
 * Date: 2014-07-14
 * Class: Research
 * Spec: Array (3D)
 */

#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <ostream>

#include "Vector.h"
#include "Utility.h"
#include "LinearAlgebra.h"

/* Also contains Light Grid */
namespace lux
{

class Coord
{
   public:
      Coord() : i(0), j(0), k(0), valid(false) {};

      Coord(int a, int b, int c) : i(a), j(b), k(c), valid(true) {};
      Coord(lux::Vector v) : i(v[0]), j(v[1]), k(v[2]), valid(true) {};
      Coord(const Coord& C) : i(C.i), j(C.j), k(C.k), valid(C.valid) {};

      bool isValid(int nx, int ny, int nz) const
      {
         bool bx = (i >= 0 && i <= nx-1);
         bool by = (j >= 0 && j <= ny-1);
         bool bz = (k >= 0 && k <= nz-1);

         return bx && by && bz;
      }

      Coord Clamp(int nx, int ny, int nz) const
      {
         Coord ndx;
         ndx.i = (i < 0) ? 0 : (i >= nx) ? nx-1 : i;
         ndx.j = (j < 0) ? 0 : (j >= ny) ? ny-1 : j;
         ndx.k = (k < 0) ? 0 : (k >= nz) ? nz-1 : k;

         return ndx;
      }

      Coord Wrap(int nx, int ny, int nz) const
      {
         Coord ndx(*this);
         Coord size(nx, ny, nz);

         for(size_t a = 0; a < 3; ++a)
         {
            if( ndx[a] < 0 )
               ndx[a] = size[a] + ndx[a];

            else if( ndx[a] >= size[a] )
               ndx[a] = ndx[a] - size[a];
         }

         return ndx;
      }

      bool operator==(const Coord C)
      { return (i == C.i && j == C.j && k == C.k); }

      bool operator!=(const Coord C)
      { return (i != C.i || j != C.j || k != C.k); }
      
      // Add Coord
      const Coord operator+(const Coord& v) const
      { return Coord(i+v.i, j+v.j, k+v.k); }

      // Subtract Coord
      const Coord operator-(const Coord& v) const
      { return Coord(i-v.i, j-v.j, k-v.k); }

      // Indexing
      const int& operator[](const size_t ndx) const
      {
         size_t mod = ndx % 3;
         return (mod == 0) ? i : (mod == 1) ? j : k;
      }

      int& operator[](const size_t ndx)
      {
         size_t mod = ndx % 3;
         return (mod == 0) ? i : (mod == 1) ? j : k;
      }

      // Assignment
      Coord& operator=(const Coord& v)
      { i = v.i; j = v.j; k = v.k; return *this; }
      
      // Output
      friend std::ostream& operator<<(std::ostream &out, const Coord& v)
      {
         out << v[0] << " " << v[1] << " " << v[2];
         return out;
      }

      int i, j, k;
      bool valid;
};

template <typename T>
class Array
{
   public:
      typedef T ArrayDataType;

      Array() {};

      Array(lux::Vector urc, lux::Vector llc, lux::Vector res, ArrayDataType def)
      {
         URC = urc;
         LLC = llc;
         _res = res;
         _default = def;

         _Nx = 1 + ceil((URC.x() - LLC.x()) / _res[0]);
         _Ny = 1 + ceil((URC.y() - LLC.y()) / _res[1]);
         _Nz = 1 + ceil((URC.z() - LLC.z()) / _res[2]);

         data = new ArrayDataType[_Nx * _Ny * _Nz];
      }

      ~Array() { delete [] data; }

      const lux::Vector mapTo(lux::Vector i) const
      {
         i = i - LLC;
         i = i / _res;

         return i;
      }

      const lux::Vector evalP(Coord ndx) const
      {
         return LLC + lux::Vector(
            static_cast<float>(ndx.i) * _res.x(), 
            static_cast<float>(ndx.j) * _res.y(), 
            static_cast<float>(ndx.k) * _res.z()
            );
      }

      const ArrayDataType eval(const Coord ndx) const
      {
         if(ndx.isValid(_Nx-1, _Ny-1, _Nz-1))
            return data[ndx.i + _Nx * (ndx.j + _Ny * ndx.k)];
         
         return _default;
      }

      ArrayDataType& eval(const Coord ndx)
      {
         if(ndx.isValid(_Nx-1, _Ny-1, _Nz-1))
            return data[ndx.i + _Nx * (ndx.j + _Ny * ndx.k)];

         return _default;
      }

      const void set(const Coord ndx, const ArrayDataType value) const
      {
         //util::gettype(value);
         if(ndx.isValid(_Nx, _Ny, _Nz))
            data[ndx.i + _Nx * (ndx.j + _Ny * ndx.k)] = value;
      }

      const bool testBoundingBox(const lux::Vector& P) const
      {
         if(P <= URC && P >= LLC)
            return true;
         return false;
      }

      void clear()
      {
         for(int i = 0; i < _Nx * _Ny * _Nz; i++)
            data[i] = _default;
      }

      const ArrayDataType get_default() const { return _default; }
      void set_default(ArrayDataType n_def) { _default = n_def; }

      ArrayDataType& operator[](const Coord ndx)
      {
         if(ndx.isValid(_Nx-1, _Ny-1, _Nz-1))
            return data[ndx.i + _Nx * (ndx.j + _Ny * ndx.k)];

         return _default;
      }

      const ArrayDataType& operator[](const Coord ndx) const
      {
         if(ndx.isValid(_Nx-1, _Ny-1, _Nz-1))
            return data[ndx.i + _Nx * (ndx.j + _Ny * ndx.k)];
         
         return _default;
      }

      const ArrayDataType& operator[](size_t ndx) const
      {
         if( ndx < (_Nx * _Ny * _Nz) )
            return data[ndx];
         
         return _default;
      }

      const int toRawIndex(lux::Coord ndx) const
      {
         return ndx.i + _Nx * (ndx.j + _Ny * ndx.k);
      }

      const lux::Vector getURC() const { return URC; }
      const lux::Vector getLLC() const { return LLC; }
      const lux::Vector getRes() const { return _res; }

      const int getNx() const { return _Nx-1; }
      const int getNy() const { return _Ny-1; }
      const int getNz() const { return _Nz-1; }
      const lux::Coord getN() const { return lux::Coord(_Nx-1, _Ny-1, _Nz-1); }

   private:
      ArrayDataType* data;
      int _Nx, _Ny, _Nz;
      lux::Vector LLC, URC, _res;
      ArrayDataType _default;

};

}

#endif
