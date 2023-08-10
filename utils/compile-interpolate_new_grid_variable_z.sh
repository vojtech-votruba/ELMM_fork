gfortran -O3 -g -fbacktrace -cpp -Wall -Wno-maybe-uninitialized -Wno-unused-dummy-argument -fsanitize=address,undefined \
../../bspline-fortran/src/bspline_sub_module.f90 \
../../bspline-fortran/src/bspline_oo_module.f90 \
../src/lists.f90  \
../src/strings.f90 \
../src/stop.f90 \
../src/parameters.f90 \
../src/arrayutilities.f90 \
../src/boundaries.f90 \
../src/scalarboundaries.f90 \
interpolate_new_grid_variable_z.f90 \
-o ~/bin/interpolate_new_grid_variable_z \
"$@"
