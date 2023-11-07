# ELMM - Extended Large Eddy Microscale Model #

### Poděkování SFG a vedoucímu práce Vladimíru Fukovi, Ph.D. /  Thanks to the SFG and to the supervisor of this project Vladimír Fuka, Ph.D. ###
Děkuji SFG za možnost vyzkoušení si vědecké práce v rámci studia na MFF UK, děkuji také obecně za finanční i jinou podporu řešení projektu.
Dále děkuji svému vedoucímu projektu panu Vladimíru Fukovi, Ph.D., který mi trpělivě vysvětloval všechno, čemu jsem nerozuměl.

Zbytek tohoto README souboru byl převzat z repozitáře ELMM pana doktora Fuky.

I thank the SFG (Student Faculty Grant) for an opportunity to try out real scientific work while studying at CUNI MFF, I also thank the SFG for all the support (financial and other) it gave me while I was working on this project.
In addition a huge thank you belongs to my supervisor Vladimír Fuka, Ph.D., who explained to me all the things I didn't understand with patience. 

The rest of this README file was copied from ELMM repository of Doctor Fuka.

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
