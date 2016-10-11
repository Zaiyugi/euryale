/* Date: 2014-06-16
 */
#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <cstdlib>
#include <cstdio>
#include <stdint.h>

#include "Vector.h"
#include "Matrix.h"
#include "LinearAlgebra.h"
using namespace std;

namespace lux
{

class Tensor
{
   public:
      Vector pos;

      Vector e_vector[3];
      double e_value[3];
      
      uint64_t id;
      bool is_deleted;

      Tensor() : pos(), id(0), is_deleted(false)
      {
         e_vector[0] = Vector(1,0,0);
         e_vector[1] = Vector(0,1,0);
         e_vector[2] = Vector(0,0,1);

         e_value[0] = e_value[1] = e_value[2] = 1;
      }

      Tensor(double e1, double e2, double e3, Vector ev1, Vector ev2, Vector ev3) : pos(), id(0), is_deleted(false)
      {
         e_vector[0] = ev1;
         e_vector[1] = ev2;
         e_vector[2] = ev3;

         e_value[0] = e1;
         e_value[1] = e2;
         e_value[2] = e3;
      }

      lux::Matrix getMatrix()
      {
         Matrix R(e_vector[0], e_vector[1], e_vector[2]);

         Matrix S(
            e_value[0],0,0,
            0,e_value[1],0,
            0,0,e_value[2]
         );

         Matrix M = R * S * R.transpose();

         return M;
      }

      lux::Matrix transform()
      {
         Matrix R(e_vector[0], e_vector[1], e_vector[2]);

         Matrix S(
            std::sqrt(e_value[0]),0,0,
            0,std::sqrt(e_value[1]),0,
            0,0,std::sqrt(e_value[2])
         );

         R = R * S.inverse();

         return R;
      }

      void setValues(const double e1, const double e2, const double e3)
      {
         e_value[0] = e1;
         e_value[1] = e2;
         e_value[2] = e3;
      }

      void setVectors(const Vector ev1, const Vector ev2, const Vector ev3)
      {
         e_vector[0] = ev1;
         e_vector[1] = ev2;
         e_vector[2] = ev3;
      }

      const double eigen_value(const size_t id) const
      { return e_value[(id-1) % 3]; }

      const Vector eigen_vector(const size_t id) const
      { return e_vector[(id-1) % 3]; }

      bool is_isotropic() const
      {
         bool xy = std::fabs(e_value[0] - e_value[1]) < TOL;
         bool xz = std::fabs(e_value[0] - e_value[2]) < TOL;
         bool yz = std::fabs(e_value[1] - e_value[2]) < TOL;

         return xy && xz && yz;
      }

   private:
      double TOL = 5.0e-6;
};

}

#endif
