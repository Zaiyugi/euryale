#ifndef __RIGIDBODY_H__
#define __RIGIDBODY_H__

#include <cstdlib>
#include <cmath>

#include "Vector.h"
#include "Vector2d.h"

namespace sim
{

class RigidBody
{
   public:
      lux::Vector2d _p;
      lux::Vector2d _v;
      lux::Vector2d _a;

      float _radius;

      lux::Vector _color;

      RigidBody() : _id(-1) {}
      RigidBody(lux::Vector2d p, lux::Vector2d v, lux::Vector2d a, float r) :
         _p(p), _v(v), _a(a),
         _radius(r), _color(lux::Vector(1.0)),
         _id(-1)
      {}

      RigidBody(const RigidBody& rb) :
         _p(rb._p), _v(rb._v), _a(rb._a),
         _radius(rb._radius), _color(rb._color),
         _id(rb._id)
      {}

      const int ID() const { return _id; }
      void setID(int id) { _id = id; }

   private:
      int _id;

};

}

#endif
