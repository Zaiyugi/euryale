/* Date: 2016-02-23
 */
#include "SPHEngine.h"

namespace sim
{

void SPHEngine::initializeRandomDevice(const float min, const float max)
{
   global_gen = std::mt19937(global_rd());
   global_dis = std::uniform_real_distribution<>(min, max);

   std::cout << "Entropy of random_device: " << global_rd.entropy() << "\n";
}

void SPHEngine::genMaterial(lux::Vector C, float m, float bd, float ps, float pg, float v, float ve)
{
   _materials.push_back( std::unique_ptr<Material>(new Material()) );

   // Material Properties
   _materials.back()->_color = C;
   _materials.back()->_mass = m;
   _materials.back()->_reference_density = bd;
   _materials.back()->_pressure_scale = ps;
   _materials.back()->_pressure_gamma = pg;
   _materials.back()->_viscosity_coeff = v;
   _materials.back()->_viscosity_epsilon = ve;
}

void SPHEngine::addMaterial(std::unique_ptr<Material> mat)
{
   _materials.push_back( std::move(mat) );
}

void SPHEngine::printCurrentMaterialAttributes()
{
   size_t cm = _current_material;
   std::cout << "Material " << cm << ": " << _materials[cm]->_name << std::endl;
   std::cout << "  Color: " << _materials[cm]->_color << std::endl;
   std::cout << "  Mass: " << _materials[cm]->_mass << std::endl;
   std::cout << "  Reference Density: " << _materials[cm]->_reference_density << std::endl;
   std::cout << "  Pressure Scale: " << _materials[cm]->_pressure_scale << std::endl;
   std::cout << "  Pressure Gamma: " << _materials[cm]->_pressure_gamma << std::endl;
   std::cout << "  Viscosity: " << _materials[cm]->_viscosity_coeff << std::endl;
   std::cout << "  Viscosity Epsilon: " << _materials[cm]->_viscosity_epsilon << std::endl;
}

void SPHEngine::update()
{
   // Update
   if( _flags.use_leapfrog )
      leapfrog(_dt);
   else if( _flags.use_sixth )
      sixth(_dt);

   // Compile stats
   _stats.max_velocity = -1.0f;
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      float v = _particles[i]->_v.magnitude();
      _stats.max_velocity = std::max(_stats.max_velocity, v);
   }

   buildGrids();
}

// Integration Schemes

// Leapfrog Integration -- 2nd order numerical integration
void SPHEngine::leapfrog(const float dt)
{
   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      _particles[i]->_p = _particles[i]->_p + _particles[i]->_v * dt * 0.5;
   }

   applyBoundaryConditions();

   // Update binning structure using last known positions
   updateOccupancyVolume();

   // Calculate densities for all particles
   _force->calculateDensity(_particles, _occupancy_volume.get());

   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      _particles[i]->_a = _force->eval(_particles, _occupancy_volume.get(), *(_particles[i]), dt);
   }

   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      _particles[i]->_v = _particles[i]->_v + _particles[i]->_a * dt;
   }

   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      _particles[i]->_p = _particles[i]->_p + _particles[i]->_v * dt * 0.5;
   }

   applyBoundaryConditions();
   testBoundaries();

}

void SPHEngine::sixth(const float dt)
{
   float a = 1.0 / (4.0 - std::pow(4.0, 1.0/3.0));
   float b = 1.0 - 4.0 * a;
   a *= dt;
   b *= dt;

   leapfrog(a);
   leapfrog(a);
   leapfrog(b);
   leapfrog(a);
   leapfrog(a);
}

// Boundaries

void SPHEngine::applyBoundaryConditions()
{
   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      Particle &P = *(_particles[i]);

      for(int c = 0; c < 2; c++)
      {
         int d = (c + 1) % 2;
         if(P._p[c] < _llc[c])
         {
            P._p[c] = 2*_llc[c] - P._p[c];

            P._v[c] = -_cor * P._v[c];
            P._v[d] = _cof * P._v[d];
         }

         if(P._p[c] > _urc[c])
         {
            P._p[c] = 2*_urc[c] - P._p[c];

            P._v[c] = -_cor * P._v[c];
            P._v[d] = _cof * P._v[d];
         }
      }

   }
}

void SPHEngine::injectParticles(const size_t amount)
{
   for(size_t i = 0; i < amount; i++)
   {
      double x = global_dis(global_gen);
      double y = global_dis(global_gen);

      Particle *p_i = new Particle();
      p_i->_p = lux::Vector2d(x, y);
      p_i->_v = lux::Vector2d(0, 0);
      p_i->_radius = _radius;

      p_i->_mat = _materials[_current_material].get();

      p_i->_dead = false;

      p_i->_id = global_id++;
      _particles.push_back(p_i);
   }
}

void SPHEngine::injectParticle(const float x_s, const float y_s)
{
   float x = x_s * (_urc[0] - _llc[0]) + _llc[0];
   float y = y_s * (_urc[1] - _llc[1]) + _llc[1];

   Particle *p_i = new Particle();
   p_i->_p = lux::Vector2d(x, y);
   p_i->_v = lux::Vector2d(0, 0);
   p_i->_radius = _radius;

   p_i->_mat = _materials[_current_material].get();

   p_i->_dead = false;

   p_i->_id = global_id++;
   _particles.push_back(p_i);
}

void SPHEngine::reset()
{
   for(auto particle : _particles)
   {
      delete particle;
   }
   _particles.clear();
}

void SPHEngine::resetParticle(Particle& P)
{
   double x = global_dis(global_gen);
   double y = global_dis(global_gen);

   P._p = lux::Vector2d(x, y);
   P._v = lux::Vector2d(0, 0);
   P._radius = _radius;

   P._mat = _materials[_current_material].get();

   P._dead = false;
}

void SPHEngine::testBoundaries()
{
   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      // If a particle has escaped or has an invalid position, reset it
      lux::Vector2d P = _particles[i]->_p;
      if( !(P < _urc && P > _llc) || P.isInf() || P.isNan() )
         resetParticle(*_particles[i]);
   }
}

void SPHEngine::updateOccupancyVolume()
{
   _occupancy_volume->reset();

   lux::Vector2d hi(std::numeric_limits<int>::min());
   lux::Vector2d lo(std::numeric_limits<int>::max());
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      // Bin: min/max bounds
      for(int j = 0; j < 2; ++j)
      {
         hi[j] = std::max(_particles[i]->_p[j], hi[j]);
         lo[j] = std::min(_particles[i]->_p[j], lo[j]);
      }
   }

   // Resize volume
   _occupancy_volume->resize(hi, lo, lux::Vector2d(_radius));

   // Sort particles into bins
   for(size_t i = 0; i < _particles.size(); ++i)
   {
      lux::Vector2d ndx = _occupancy_volume->mapTo(_particles[i]->_p);

      _occupancy_volume->eval(ndx.x(), ndx.y()).push_back(i);
   }

}

// Invert the tank
// Make some NOISE!
void SPHEngine::flip()
{
   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); i++)
   {
      _particles[i]->_p[1] *= -1.0;
   }

}

// Rotate the tank
// Make some NOISE!
void SPHEngine::swapCoords()
{
   #pragma omp parallel for
   for(size_t i = 0; i < _particles.size(); i++)
   {
      std::swap(_particles[i]->_p[0], _particles[i]->_p[1]);
   }

}

void SPHEngine::buildGrids()
{
   std::memset(_color_grid, 0, sizeof(float) * _Nx * _Ny * 3);

   lux::Vector2d grid_urc(_Nx, _Ny), grid_llc;
   for(int i = 0; i < _Nx; i++)
   {
      #pragma omp parallel for
      for(int j = 0; j < _Ny; j++)
      {
         lux::Vector2d P(i, j);
         P = castToSpace(P, grid_urc, grid_llc, _urc, _llc);

         lux::Vector2d mapped = _occupancy_volume->mapTo(P);
         int px = mapped[0];
         int py = mapped[1];

         lux::Vector color;
         int nP = 0;
         for(int s = px-1; s < px+2; s++)
         {
            for(int t = py-1; t < py+2; t++)
            {
               const std::vector<size_t>& cell = _occupancy_volume->eval(s, t);

               for( size_t k = 0; k < cell.size(); ++k)
               {
                  lux::Vector2d dist = (P - _particles[cell[k]]->_p);
                  float q = (dist[0]*dist[0] + dist[1]*dist[1]) / std::pow(_radius/3.0f, 2.0f);
                  if( q < 1.0 )
                  {
                     q = 1.0 - q;
                     color += _particles[cell[k]]->_mat->_color * q;
                     nP++;
                  }
               }
            }
         }

         int ndx = coord(i, j);
         color = color / static_cast<float>(nP);
         _color_grid[3*ndx] = color[0];
         _color_grid[3*ndx+1] = color[1];
         _color_grid[3*ndx+2] = color[2];
      }
   }

}

lux::Vector2d SPHEngine::castToSpace(
   lux::Vector2d x,
   lux::Vector2d from_urc, lux::Vector2d from_llc,
   lux::Vector2d to_urc, lux::Vector2d to_llc)
{
   lux::Vector2d dx = (x - from_llc) / (from_urc - from_llc);

   lux::Vector2d xT;
   xT[0] = dx[0] * (to_urc - to_llc)[0];
   xT[1] = dx[1] * (to_urc - to_llc)[1];

   return xT + to_llc;
}

// Interaction
void SPHEngine::useLeapfrog()
{
   _flags.use_leapfrog = 1;
   _flags.use_sixth = 0;
}

void SPHEngine::useSixth()
{
   _flags.use_leapfrog = 0;
   _flags.use_sixth = 1;
}

}

// End definitions
