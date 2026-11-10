// Writing numbers too large for their output column.
//
// Every output file writes each value right-aligned in a fixed width, 60 or
// 20 characters, in fixed point with 10 or 4 decimals, and relies on the
// padding to separate one value from the next. A value whose fixed-point form
// is as long as the width gets no padding and is written against the previous
// value with no space between them. In the results file that first happens at
// period 13,273 of the model's own parameter file, seed 1, when nominal GDP
// (column 39) passes 1e49; read back, the two values become one malformed
// number or two wrong ones. docs/DSK_LONG_HORIZON.md has the measurement.
//
// The fix is a number formatter installed as the global locale at the start
// of main, so every file stream the model opens uses it without any write
// being touched. It changes only the values that would have run into their
// neighbour: those are written in scientific notation with 17 significant
// digits, which is enough to give back the same double, after at least one
// space. Every value that fitted before is written by the standard formatter
// exactly as before, so output that was readable is byte-identical.
//
// Only fixed-point output with a width set is affected. Values written without
// a width, as in the rescaling logs, are separated by explicit spaces already.

#ifndef DSK_SFC_WIDE_VALUES_H
#define DSK_SFC_WIDE_VALUES_H

#include <cmath>
#include <cstdio>
#include <ios>
#include <locale>

// Enough significant digits for any double to be read back exactly.
const int wide_value_significant_digits=17;

class WideValueNumPut : public std::num_put<char>
{
protected:
   iter_type do_put(iter_type out, std::ios_base &stream, char fill, double value) const
   {
      std::streamsize width=stream.width();
      bool fixed=(stream.flags() & std::ios_base::floatfield)==std::ios_base::fixed;
      if (width>0 && fixed && std::isfinite(value))
      {
         int fixed_length=std::snprintf(NULL,0,"%.*f",(int)stream.precision(),value);
         if (fixed_length>=width)
         {
            char text[64];
            int length=std::snprintf(text,sizeof text," %.*e",wide_value_significant_digits-1,value);
            stream.width(0);
            for (std::streamsize pad=width-length; pad>0; pad--)
            {
               *out++=fill;
            }
            for (int i=0; i<length; i++)
            {
               *out++=text[i];
            }
            return out;
         }
      }
      return std::num_put<char>::do_put(out,stream,fill,value);
   }
};

// Call once at the top of main, before any output file is opened. A locale
// holding a facet of its own has no name, so this leaves the C locale, and
// with it printf, as it was.
inline void INSTALL_WIDE_VALUE_FORMAT(void)
{
   std::locale::global(std::locale(std::locale(),new WideValueNumPut));
}

#endif
