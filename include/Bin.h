/* Creation Date: 2014-07-14
 * Last edit: 2016-02-23
 */
#ifndef __BIN_H__
#define __BIN_H__

#include <iostream>
#include <cstdlib>
#include <vector>

#include "Particle.h"

const static double DIST_TOL = 1e-8;

namespace sim
{

class Bin
{

   public:

      Bin() {};

      Bin(size_t idx, size_t idy, size_t idz) {};

      size_t addTo(Particle* P)
      {
         _particles.push_back(P);
         return _particles.size() - 1;
      }

      void removeFrom(size_t ndx)
      {
         if(ndx >= _particles.size() || _particles.size() == 0)
         {
            std::cerr << "Invalid index for removal from bin " << ndx << "\n";
            return;
         }

         if(ndx != _particles.size()-1)
         {
            size_t end = _particles.size() - 1;

            std::swap(_particles[ndx]->id_in_bin, _particles[end]->id_in_bin);
            std::swap(_particles[ndx], _particles[end]);
         }

         _particles.pop_back();

      }

      Particle* get(size_t ndx) { return _particles[ndx]; }
      size_t size() { return _particles.size(); }

   private:

      std::vector<Particle*> _particles;

};

}

#endif
