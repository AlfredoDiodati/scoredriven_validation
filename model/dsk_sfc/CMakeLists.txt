cmake_minimum_required(VERSION 3.6)

# Optimised rather than Debug. The shipped Debug build carries no optimisation
# flag at all and runs 3.8 times slower for output that is identical byte for
# byte; equivalence_test.c checks that claim against an unmodified build.
# -O2 rather than -O3: -O3 measured no faster here. -march=native is
# deliberately absent, since it lets the compiler contract multiply and add
# into one instruction, which changes the last digits of every result.
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_C_FLAGS_RELEASE   "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fno-math-errno -msse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -msse")

# Link-time optimisation. newmat defines its element accessors out of line, so
# without it every one of the hundreds of millions of A(t,i) reads per run is a
# real function call.
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

project(dsk_SFC)

set(SOURCE_FILES
        include/CLI11/CLI11.hpp
        newmat10/bandmat.cpp
        newmat10/cholesky.cpp
        newmat10/controlw.h
        newmat10/evalue.cpp
        newmat10/fft.cpp
        newmat10/hholder.cpp
        newmat10/jacobi.cpp
        newmat10/myexcept.cpp
        newmat10/newfft.cpp
        newmat10/newmat1.cpp
        newmat10/newmat2.cpp
        newmat10/newmat3.cpp
        newmat10/newmat4.cpp
        newmat10/newmat5.cpp
        newmat10/newmat6.cpp
        newmat10/newmat7.cpp
        newmat10/newmat8.cpp
        newmat10/newmat9.cpp
        newmat10/newmatap.h
        newmat10/newmatex.cpp
        newmat10/newmatnl.cpp
        newmat10/newmatnl.h
        newmat10/newmatrc.h
        newmat10/newmatrm.cpp
        newmat10/newmatrm.h
        newmat10/precisio.h
        newmat10/solution.cpp
        newmat10/solution.h
        newmat10/sort.cpp
        newmat10/submat.cpp
        newmat10/svd.cpp
        newmat10/newmatio.h
        newmat10/newmat.h
        newmat10/include.h
		
        auxiliary/betadev.cpp
        auxiliary/betadev.h
        auxiliary/bnldev.cpp
        auxiliary/bnldev.h
        auxiliary/gammln.cpp
        auxiliary/gammln.h
        auxiliary/gasdev.cpp
        auxiliary/gasdev.h
        auxiliary/boolean.h
        auxiliary/myexcept.h
        auxiliary/ran1.cpp
        auxiliary/ran1.h

        dsk_sfc_main.cpp
        dsk_sfc_functions.h
        dsk_sfc_globalvars.h
	dsk_sfc_flags.h
        dsk_sfc_include.h
	dsk_sfc_inits.h
	dsk_sfc_parameters.h
        modules/module_energy_sfc.cpp
        modules/module_energy_sfc.h
        modules/module_climate_sfc.cpp
        modules/module_climate_sfc.h
        modules/module_macro_sfc.cpp
        modules/module_macro_sfc.h
	modules/module_finance_sfc.cpp
	modules/module_finance_sfc.h
        modules/WITCH_input/WITCH_input.cpp
        modules/WITCH_input/WITCH_input.h
    )

add_executable(dsk_SFC ${SOURCE_FILES})
if (NOT APPLE)
	target_link_libraries(dsk_SFC -static)
endif(NOT APPLE)