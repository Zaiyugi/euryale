#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <cstdlib>
#include <cmath>
#include <cstring>

namespace sim
{

class Material
{
   public:
      std::string _name;
      lux::Vector _color;

      float _mass;
      float _reference_density;

      // Pressure
      float _pressure_scale;
      float _pressure_gamma;

      // Viscosity
      float _viscosity_coeff;
      float _viscosity_epsilon;

      Material() :
         _name("Default"),
         _color(lux::Vector(1.0, 1.0, 1.0)),
         _mass(1.0),
         _reference_density(1.0),
         _pressure_scale(1.0),
         _pressure_gamma(3.0),
         _viscosity_coeff(1.0),
         _viscosity_epsilon(0.1)
      {}

      float pressure(float density)
      {
         return _pressure_scale * (std::pow(density / _reference_density, _pressure_gamma) - 1.0);
      }

};

}

#endif
