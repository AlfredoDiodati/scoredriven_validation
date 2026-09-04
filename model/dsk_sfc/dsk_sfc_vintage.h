// The per-vintage arrays as one block of memory each.
//
// They were a vector of vectors of vectors, so every row of firms was its own
// heap allocation and consecutive rows landed wherever the allocator put them.
// The loops that matter walk vintages with the firm held fixed, which touches
// one element of each row in turn: 400 loads from 400 addresses the processor
// cannot guess, once per firm per period, in the four functions that are most
// of the run.
//
// Laid out as one block the same walk reads addresses a fixed distance apart,
// which is a stride the hardware prefetcher recognises and can run ahead of.
// The index order is unchanged - period, supplier, firm - so every expression
// of the form X[tt-1][i-1][j-1] still means what it meant and still reads the
// same value in the same order. Only where the memory sits changes.
//
// The subscripts return small proxies rather than references, which is what
// lets the 185 element accesses in the model stay exactly as they were written.
// One asymmetry is deliberate and worth knowing about: copying a row proxy
// copies the pointer, so binding one to a name is free, while assigning to a
// row copies the firms' values, which is what g_c[tt-1][i-1]=g[tt-1][i-1] has
// always meant.

#ifndef DSK_SFC_VINTAGE_H
#define DSK_SFC_VINTAGE_H

#include <vector>
#include <cstring>
#include <cstddef>

template <class T> struct VintageRow
{
   T* p;
   int n;

   T& operator[](int j) const { return p[j]; }

   // Assigning a row copies the values, as the vector of vectors did.
   const VintageRow& operator=(const VintageRow& other) const
   {
      std::memcpy(p, other.p, (std::size_t)n * sizeof(T));
      return *this;
   }
};

template <class T> struct VintagePlane
{
   T* p;
   int n2;

   VintageRow<T> operator[](int i) const { return VintageRow<T>{p + (std::size_t)i * n2, n2}; }
};

template <class T> struct Vintage3D
{
   std::vector<T> data;
   int d1, d2;

   Vintage3D() : d1(0), d2(0) {}

   void resize3(int periods, int suppliers, int firms)
   {
      d1 = suppliers;
      d2 = firms;
      data.assign((std::size_t)periods * suppliers * firms, T());
   }

   VintagePlane<T> operator[](int tt) const
   {
      return VintagePlane<T>{const_cast<T*>(data.data()) + (std::size_t)tt * d1 * d2, d2};
   }
};

#endif
