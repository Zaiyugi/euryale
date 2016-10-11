#include "SPHForce.h"

namespace sim
{

lux::Vector2d SPHForce::eval(const std::vector<Particle*>& sphParticles, const SPHVolume* bins, Particle& A, const float dt)
{
   // Iterate over bins adjacent to p's bin
   lux::Vector2d F = _gravity;
   float p_a = A._pressure;

   lux::Vector2d ndx = bins->mapTo(A._p);
   int nx = static_cast<int>(ndx[0]);
   int ny = static_cast<int>(ndx[1]);

   for(int i = nx-1; i < nx+2; ++i)
   {
      for(int j = ny-1; j < ny+2; ++j)
      {
         const std::vector<size_t>& cell = bins->eval(i, j);

         for( size_t n = 0; n < cell.size(); ++n)
         {
            Particle& B = *sphParticles[cell[n]];

            // Calculate pressure
            float p_b = B._pressure;

            float vis = B.viscosity(A, dt);

            // Calculate gradient of weight function
            F -= B._mat->_mass * (p_b + p_a + vis) * B.weight_gradient(A._p);
         }
      }
   }

   return (1.0 / A._mat->_mass) * F;
}

void SPHForce::calculateDensity(const std::vector<Particle*>& sphParticles, const SPHVolume* bins)
{
   #pragma omp parallel for
   for(size_t x = 0; x < sphParticles.size(); ++x)
   {
      Particle& B = *sphParticles[x];

      B._density = 0;

      // Using adjacent bins, calculate density for particle B
      lux::Vector2d ndx = bins->mapTo(B._p);
      int nx = static_cast<int>(ndx[0]);
      int ny = static_cast<int>(ndx[1]);
      for(int i = nx-1; i < nx+2; ++i)
      {
         for(int j = ny-1; j < ny+2; ++j)
         {
            const std::vector<size_t>& cell = bins->eval(i, j);

            for( size_t n = 0; n < cell.size(); ++n)
            {
               Particle& A = *sphParticles[cell[n]];

               B._density += A._mat->_mass * A.weight(B._p);
            }
         }
      }

      B._pressure = B.getPressure();
   }
}





// end namespace
}
