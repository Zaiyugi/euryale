#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdint.h>

#include "Vector2d.h"
#include "Array.h"
#include "Material.h"

namespace sim
{

class Particle
{
   public:
      lux::Vector2d _p;
      lux::Vector2d _v;
      lux::Vector2d _a;

      float _radius;

      float _density;
      float _pressure;

      size_t _id;
      bool _dead;

      Material* _mat;

      Particle() {}

      Particle(const Particle& P) :
         _p(P._p), _v(P._v), _a(P._a),
         _radius(P._radius),
         _density(P._density),
         _pressure(P._pressure),
         _id(P._id),
         _dead(P._dead),
         _mat(P._mat)
      {}

      float weight(lux::Vector2d b)
      {
         float q = (b - _p).magnitude() / _radius;

         float ret = 0;
         if( q < 1 )
         {
            ret = 10.0 / (M_PI * _radius * _radius);
            ret *= std::pow(1.0 - q, 3.0);
         }

         return ret;
      }

      lux::Vector2d weight_gradient(lux::Vector2d b)
      {
         float q = (b - _p).magnitude() / _radius;

         lux::Vector2d ret = (b - _p).unitvector();

         float s = 0;
         if( q < 1 )
         {
            s = -30.0 / (M_PI * std::pow(_radius, 3.0));
            s *= std::pow(1.0 - q, 2.0);
         }

         return ret * s;
      }

      float viscosity(const Particle& b, const float dt)
      {
         lux::Vector2d Vab = _v - b._v;
         lux::Vector2d Xab = _p - b._p;
         float xmag = Xab.magnitude();

         float ret = 0;
         float c = Vab * Xab;

         if( c < 0 )
         {
            float u = c / (xmag * xmag + _mat->_viscosity_epsilon * _radius * _radius);
            float v = _mat->_viscosity_coeff * _radius / (_density + b._density);

            ret = -v * u;
         }

         return ret;
      }

      float getPressure() { return _mat->pressure(_density) / (_density * _density); }
};

}

#endif
