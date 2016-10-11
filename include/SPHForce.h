#ifndef __SPHFORCE_H__
#define __SPHFORCE_H__

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdint.h>
#include <vector>

#include "Vector2d.h"
#include "Particle.h"
#include "SPHVolume.h"

namespace sim
{

class SPHForce
{
   public:
      lux::Vector2d _gravity;

      // Constructor(s)
      SPHForce() {}

      // Methods
      lux::Vector2d eval(const std::vector<Particle*>& sphParticles, const SPHVolume* bins, Particle& A, const float dt);

      void calculateDensity(const std::vector<Particle*>& sphParticles, const SPHVolume* bins);

};

}

#endif
