/* Name: Zachary Shore
 * Created: 2016-02-24
 * Edited: 2016-02-24
 * Spec: Utility Functions
 */

#include "Utility.h" 

namespace util
{

template <typename T>
void gettype(const T obj)
{
   int* status = (int*)malloc(sizeof(int));
   char * name = abi::__cxa_demangle(typeid(obj).name(), 0, 0, status);
   if(*status == 0 && name != NULL) {
      std::cout << "type is: " << name << "\n";
      free(name);
   } else if(*status == -1)
      std::cerr << "Memory Allocation failed\n";
   else if(*status == -2)
      std::cerr << typeid(obj).name() << " is not a valid name under C++ ABI mangling rules\n";
   else
      std::cerr << "Invalid Parameter\n";

   free(status);
}

void writeKarl(std::string path, std::vector<lux::Vector> &data)
{
   std::cout << "\tWriting: " << path << "\n";
   std::ofstream ofs(path.c_str(), std::ofstream::binary | std::ofstream::out);

   size_t size = data.size();
   ofs << size;

   for(size_t i = 0; i < data.size(); ++i)
      data[i].write(ofs);

   ofs.close();
}

void writeObj(std::string path, std::vector<lux::Vector> &data)
{
   std::cout << "\tWriting: " << path.c_str() << "\n";
   std::ofstream ofs(path.c_str(), std::ofstream::out);

   for(auto vec : data)
   {
      ofs << "v ";
      vec.write(ofs);
      ofs << "\n";
   }

   ofs.close();
}

float round_to_nearest(float number, float base)
{
   if(base != 0 && number != 0)
   {
      float sign = number > 0 ? 1 : -1;
      number *= sign;
      number /= base;
      int fixedPoint = static_cast<int>(ceil(number));
      number = fixedPoint*base;
      number *= sign;
   }
   return number;
}

// Template Functions (Strange Magic...)

template <typename T> int sgn(T val)
{
   return (T(0) < val) - (val < T(0));
}

template <typename T> void swap(T &a, T &b)
{
   T temp = a;
   a = b;
   b = temp;
}

}

