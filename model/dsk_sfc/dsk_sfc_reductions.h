// Column reductions that do not build the column first.
//
// newmat evaluates M.Column(i) into a temporary before reducing it, and since
// a Matrix is stored row-major that copies the column out one row at a time.
// In the two ranking loops that rank each bank's borrowers by debt service it
// happens once per firm per bank per period, which made it the single largest
// cost in a run: 264 million row copies, around 30% of the time. These walk
// the column where it lies instead, by its stride through the matrix's own
// storage.
//
// The scan order, the comparisons and the tie-breaking are newmat's own, from
// GeneralMatrix::Maximum and GeneralMatrix::Minimum1 in newmat10/newmat8.cpp,
// so the value returned and the index reported are the ones the library
// returned. Maximum keeps the first of equal maxima, Minimum1 the last of
// equal minima, and neither lets a NaN win.

#ifndef DSK_SFC_REDUCTIONS_H
#define DSK_SFC_REDUCTIONS_H

#include "newmat10/newmat.h"

static inline Real ColumnMaximum(const Matrix& M, int col)
{
   const int n = M.Nrows(), stride = M.Ncols();
   const Real* s = M.Store() + (col - 1);
   Real maxval = *s;
   for (int k = 1; k < n; k++)
   {
      const Real a = s[(size_t)k * stride];
      if (maxval < a) maxval = a;
   }
   return maxval;
}

// index comes back 1-based, as Minimum1's own out-parameter does.
static inline Real ColumnMinimum1(const Matrix& M, int col, int& index)
{
   const int n = M.Nrows(), stride = M.Ncols();
   const Real* s = M.Store() + (col - 1);
   Real minval = *s;
   int li = n - 1;
   for (int k = 1; k < n; k++)
   {
      const Real a = s[(size_t)k * stride];
      const int l = n - 1 - k;
      if (minval >= a) { minval = a; li = l; }
   }
   index = n - li;
   return minval;
}

// Sum of one row, where it lies. M.Row(i).Sum() copies the row into a fresh
// vector first, and a Matrix is stored row-major, so the row is already
// contiguous and the copy is pure cost. The accumulation starts at zero and
// runs left to right, which is GeneralMatrix::Sum's own order, so the value is
// the one that call returned.
static inline Real RowSum(const Matrix& M, int row)
{
   const int ncols = M.Ncols();
   const Real* s = M.Store() + (size_t)(row - 1) * ncols;
   Real sum = 0.0;
   for (int k = 0; k < ncols; k++) sum += s[k];
   return sum;
}

// Minimum over a block of whole rows, with GeneralMatrix::Minimum's own scan
// and comparison, for the case where the rest of the matrix is known to hold
// values that cannot win.
static inline Real RowRangeMinimum(const Matrix& M, int firstRow, int lastRow)
{
   const int ncols = M.Ncols();
   const Real* s = M.Store() + (size_t)(firstRow - 1) * ncols;
   const size_t n = (size_t)(lastRow - firstRow + 1) * ncols;
   Real minval = *s++;
   for (size_t k = 1; k < n; k++)
   {
      const Real a = *s++;
      if (minval > a) minval = a;
   }
   return minval;
}

#endif
