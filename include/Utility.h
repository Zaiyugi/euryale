/* Name: Zachary Shore
 * Created: 2013-12-30
 * Edited: 2014-08-04
 * Spec: Utility Functions
 */

#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <cstdio>
#include <string>
#include <cxxabi.h>
#include <vector>

#include <Vector.h>

namespace util
{

template <typename T>
void gettype(const T obj);
void writeKarl(std::string path, std::vector<lux::Vector> &data);
void writeObj(std::string path, std::vector<lux::Vector> &data);
float round_to_nearest(float number, float base);

// Template Functions (Strange Magic...)
template <typename T> int sgn(T val);
template <typename T> void swap(T &a, T &b);

}

#endif

