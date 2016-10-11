/* Date: 2016-02-23
 */
#ifndef __SPHENGINE_H__
#define __SPHENGINE_H__

#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <limits>
#include <vector>
#include <omp.h>

#include "Utility.h"
#include "Vector.h"
#include "Vector2d.h"

#include "Particle.h"
#include "RigidBody.h"

#include "SPHForce.h"
#include "SPHVolume.h"

namespace sim
{

const static double TOL = 0.00005;
const static double DIST_TOL = 10e-7;

class EngineStats
{
   public:
      float max_velocity;

};

struct EngineFlags
{
   int use_leapfrog = 1;
   int use_sixth = 0;

};

class SPHEngine
{
   public:

      lux::Vector2d _urc, _llc;
      int _Nx, _Ny;
      std::vector<Particle*> _particles;

      std::vector< std::unique_ptr<Material> > _materials;
      size_t _current_material;

      SPHEngine() {};

      SPHEngine(lux::Vector2d hi, lux::Vector2d lo, float dt, float r, int nx, int ny) :
         _urc(hi), _llc(lo),
         _Nx(nx), _Ny(ny),
         _current_material(0),
         _dt(dt), _radius(r),
         _cor(1.0), _cof(1.0),
         _force(new SPHForce()),
         _occupancy_volume(new SPHVolume())
      {
         global_id = 0;

         global_gen = std::mt19937(global_rd());
         global_gen.seed(1);

         global_dis = std::uniform_real_distribution<>(0, 1);
         std::cout << "Entropy of random_device: " << global_rd.entropy() << "\n";

         _materials.push_back( std::unique_ptr<Material>(new Material()) );

         // Limit threads
         omp_set_num_threads(_thread_cnt);

         _color_grid = new float[_Nx * _Ny * 3];
         std::memset(_color_grid, 0, sizeof(float) * _Nx * _Ny * 3);
      }

      ~SPHEngine()
      {
         for(auto particle : _particles)
         {
            delete particle;
         }

         if( _color_grid )
            delete [] _color_grid;
      }

      void initializeRandomDevice(const float min, const float max);

      void genMaterial(lux::Vector C, float m, float bd, float ps, float pg, float v, float ve);
      void addMaterial(std::unique_ptr<Material> mat);
      void printCurrentMaterialAttributes();

      void update();

      // Integration schemes
      void leapfrog(const float dt);
      void sixth(const float dt);

      void applyBoundaryConditions();
      void bounceParticle(Particle& P, lux::Vector2d& B, lux::Vector2d& N);

      void injectParticles(const size_t amount);
      void injectParticle(const float x_s, const float y_s);
      void updateOccupancyVolume();

      void reset();
      void resetParticle(Particle& P);
      void testBoundaries();

      void flip();
      void swapCoords();

      void buildGrids();

      lux::Vector2d castToSpace(
         lux::Vector2d x,
         lux::Vector2d from_urc, lux::Vector2d from_llc,
         lux::Vector2d to_urc, lux::Vector2d to_llc);

      // Interaction
      void useLeapfrog();
      void useSixth();

      // Get and Set
      void setCollisionCoeffs(const float cor, const float cof) { _cor = cor; _cof = cof; }

      float getRadius() const { return _radius; }
      SPHForce* getForce() const { return _force.get(); }
      SPHVolume* getVolume() const { return _occupancy_volume.get(); }
      float* getColorField() const { return _color_grid; }
      const EngineStats& getStats() const { return _stats; }

   private:
      // Variables
      uint64_t global_id;

      float _dt, _radius;
      float _cor, _cof;
      std::unique_ptr<SPHForce> _force;
      std::unique_ptr<SPHVolume> _occupancy_volume;

      // Color Grid
      float* _color_grid;

      std::random_device global_rd;
      std::mt19937 global_gen;
      std::uniform_real_distribution<> global_dis;

      EngineFlags _flags;
      EngineStats _stats;

      const static int _thread_cnt = 3;

      // Utilities
      inline bool isValid(const int x, const int y) const
      { return ((x >= 0 && x < _Nx) && (y >= 0 && y < _Ny)); }

      int coord(const int x, const int y) const
      {
         if(!isValid(x, y))
            return -1;

         return x + _Nx * y;
      }
};

}

#endif
