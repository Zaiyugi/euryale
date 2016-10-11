/* Date: 2014-06-16
 */
#ifndef __TENSORFIELD_H__
#define __TENSORFIELD_H__

#include <iostream>
#include <cstdlib>
#include <vector>

#include "Vector.h"
#include "Matrix.h"
#include "LinearAlgebra.h"
#include "Tensor.h"

#include "Array.h"

namespace lux
{

class TensorField
{

   typedef std::vector<Tensor*> TensorVoxel;

   public:

      Array<TensorVoxel> tensors;

      TensorField() {}

      TensorField(Vector urc, Vector llc, Vector res) : tensors( urc, llc, res, TensorVoxel() ) {}

      ~TensorField()
      {
         for(int k = 0; k < tensors.getNz(); ++k)
         {
            for(int j = 0; j < tensors.getNy(); ++j)
            {
               for(int i = 0; i < tensors.getNx(); ++i)
               {
                  lux::Coord ndx(i, j, k);

                  TensorVoxel tv = tensors[ndx];
                  for(auto tensor : tv)
                  {
                     delete tensor;
                  }
               }
            }
         }

      }

      void addTensor(Tensor* t)
      {
         if(!tensors.testBoundingBox(t->pos))
         {
            std::cerr << "ERROR: (" << t->pos[0] << ", " << t->pos[1] << ", " << t->pos[2] << ") is not within bounds of tensor field\n";
            return;
         }

         Vector q = tensors.mapTo(t->pos);

         Vector ndx;
         q[0] = std::modf(q.x(), &ndx[0]);
         q[1] = std::modf(q.y(), &ndx[1]);
         q[2] = std::modf(q.z(), &ndx[2]);

         tensors[Coord(ndx)].push_back(t);
      }

      Tensor eval(Vector p)
      {
         Vector ev1, ev2, ev3;
         double e1(0), e2(0), e3(0);

         interp(p, e1, e2, e3, ev1, ev2, ev3);

         return Tensor(e1, e2, e3, ev1, ev2, ev3);
      }

      Matrix transform(Vector p)
      {
         Vector ev1, ev2, ev3;
         double e1(0), e2(0), e3(0);

         interp(p, e1, e2, e3, ev1, ev2, ev3);

         Matrix S(
            std::sqrt(e1),0,0,
            0,std::sqrt(e2),0,
            0,0,std::sqrt(e3)
         );

         Matrix R(ev1, ev2, ev3);

         R = R * S.inverse();

         return R;
      }

      int size() { return tensors.getNx() * tensors.getNy() * tensors.getNz(); }

   private:

      Vector mapTo(Vector p)
      {
         if(!tensors.testBoundingBox(p))
         {
            std::cerr << "ERROR: (" << p[0] << ", " << p[1] << ", " << p[2] << ") is not within bounds of tensor field\n";
            return Vector();
         }

         Vector q = tensors.mapTo(p);
         Vector ndx;
         q[0] = std::modf(q.x(), &ndx[0]);
         q[1] = std::modf(q.y(), &ndx[1]);
         q[2] = std::modf(q.z(), &ndx[2]);

         return ndx;
      }

      void interp(Vector p, double &e1, double &e2, double &e3, Vector &ev1, Vector &ev2, Vector &ev3)
      {
         double ws = 0.0;
         Vector q = mapTo(p);
         Coord ndx(q);

         for(int a = ndx.i-1; a < ndx.i+2; ++a)
            for(int b = ndx.j-1; b < ndx.j+2; ++b)
               for(int c = ndx.k-1; c < ndx.k+2; ++c)
               {
                  Coord aug_ndx(a, b, c);
                  if( !aug_ndx.isValid(tensors.getNx(), tensors.getNy(), tensors.getNz()) )
                     continue;

                  // aug_ndx = aug_ndx.Wrap(tensors.getNx(), tensors.getNy(), tensors.getNz());

                  int n = tensors[aug_ndx].size();
                  for(int it = 0; it < n; ++it)
                  {
                     Tensor* T = tensors[aug_ndx][it];

                     Vector r = T->pos - p;
                     double a = std::pow(1.0 / (1 + r.magnitude()), 2.0);

                     ws += a;
                     ev1 += a * T->eigen_vector(1); e1 += a * T->eigen_value(1);
                     ev2 += a * T->eigen_vector(2); e2 += a * T->eigen_value(2);
                     ev3 += a * T->eigen_vector(3); e3 += a * T->eigen_value(3);
                  }
               }

         ev1 = ev1 / ws;
         ev2 = ev2 / ws;
         ev3 = ev3 / ws;
         e1 /= (ws);
         e2 /= (ws);
         e3 /= (ws); 
      }

};

}

#endif
