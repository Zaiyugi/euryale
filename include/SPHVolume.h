/* Name: Zachary Shore
 * Date: 2014-07-14
 * Spec: SPH Occupancy Volume
 */

#ifndef __SPHVOLUME_H__
#define __SPHVOLUME_H__

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <ostream>

#include "Vector2d.h"

namespace sim
{

class SPHVolume
{
   public:
      SPHVolume() {};

      SPHVolume(lux::Vector2d urc, lux::Vector2d llc, lux::Vector2d res)
      {
         resize(urc, llc, res);
      }

      ~SPHVolume() {}

      void resize(lux::Vector2d urc, lux::Vector2d llc, lux::Vector2d res)
      {
         _urc = urc;
         _llc = llc;
         _res = res;

         _Nx = 1 + (_urc.x() - _llc.x()) / _res[0];
         _Ny = 1 + (_urc.y() - _llc.y()) / _res[1];

         //_res[0] = (_urc[0] - _llc[0]) / static_cast<float>(_Nx-1);
         //_res[1] = (_urc[1] - _llc[1]) / static_cast<float>(_Ny-1);

         _data.resize(_Nx * _Ny);
      }

      const lux::Vector2d mapTo(lux::Vector2d i) const
      {
         i = i - _llc;
         i = i / _res;

         return i;
      }

      const lux::Vector2d evalP(const int i, const int j) const
      {
         return _llc + lux::Vector2d(i * _res.x(), j * _res.y());
      }

      const std::vector<size_t>& eval(const int i, const int j) const
      {
         if(i < 0 || i >= _Nx || j < 0 || j >= _Ny)
            return _default;

         return _data[j * _Nx + i];
      }

      std::vector<size_t>& eval(const int i, const int j)
      {
         if(i < 0 || i >= _Nx || j < 0 || j >= _Ny)
            return _default;

         return _data[j * _Nx + i];
      }

      size_t ndx(const int i, const int j) const
      { return j * _Nx + i; }

      void reset()
      {
         _data.clear();
      }

      const lux::Vector2d getURC() const { return _urc; }
      const lux::Vector2d getLLC() const { return _llc; }
      const lux::Vector2d getRes() const { return _res; }

      const int getNx() const { return _Nx; }
      const int getNy() const { return _Ny; }

   private:
      std::vector< std::vector<size_t> > _data;
      std::vector<size_t> _default;
      int _Nx, _Ny;
      lux::Vector2d _urc, _llc, _res;

};

}

#endif
