# ELMM - Extended Large Eddy Microscale Model #

### About ELMM: ###

ELMM is an LES model for atmospheric flows. Especially used for atmospheric dispersion problems in complex geometry.


### Prerequisities ###

ELMM depends on libraries:

* PoisFFT [https://github.com/LadaF/PoisFFT](https://github.com/LadaF/PoisFFT)
* PFFT [https://github.com/mpip/pfft/](https://github.com/mpip/pfft/) (required for the MPI and hybrid MPI/OpenMP mode)
* BLAS and LAPACK (often optimized versions provided by your system administrators or Linux distributions)

The build scripts require:

* Scons [http://scons.org/](http://scons.org/)

### Installation ###

The scons build scripts are invoked using the `make_...` wrappers in the `src/` subdirectory.

For example:

```
#!bash

git clone https://github.com/LadaF/PoisFFT.git
git clone https://LadaF@bitbucket.org/LadaF/elmm.git
cd elmm/src
./make_release

```
or

```
#!bash

./make_mpi_release
```


### Examples ###

Examples of configuration are located in the `examples/` subdirectory.


### Contacts: ###

Repository owner: Vladimir Fuka

Issues should be reported at https://bitbucket.org/LadaF/elmm/issues?status=new&status=open

Any comments and questions are welcome by e-mail: vladimir.fuka@gmail.com. Prefer the above link for more efficient issue solving.