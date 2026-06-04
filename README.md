# slmpi_emcee

**The MPI-based Parallelization of the S-Lang MCMC Hammer**

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/mcfit/slmpi_emcee/blob/master/LICENSE)
[![Zenodo](https://img.shields.io/badge/DOI-10.5281/zenodo.4495925-blue.svg)](https://doi.org/10.5281/zenodo.4495925)

## Description

[slmpi_emcee](https://www.sternwarte.uni-erlangen.de/wiki/index.php/Emcee)
is an MPI-based Parallelization of the S-Lang *Markov chain Monte Carlo
(MCMC) Hammer* ([isis_emcee](https://github.com/mcfit/isis_emcee)) based
on algorithm proposed by [Goodman & Weare
(2010)](https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract),
implemented in Python ([emcee](https://github.com/dfm/emcee)) by
[Foreman-Mackey et al.
(2013)](https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract),
which has then been implemented in the [Interactive Spectral
Interpretation System (ISIS)](http://space.mit.edu/cxc/isis/) in S-Lang
by [M. A. Nowak](http://space.mit.edu/home/mnowak/isis_vs_xspec/),
modified by T. Beuchert, L. Corrales, and M. Kuehnel, and extended for
MPI-based Parallelization by A. Danehkar. Both the S-Lang MCMC hammer
([emcee](https://www.sternwarte.uni-erlangen.de/wiki/index.php/Emcee))
and the MPI-parallelized S-Lang MCMC hammer (`mpi_emcee`) are included
in the [Remeis
ISISscripts](http://www.sternwarte.uni-erlangen.de/isis/). It utilizes
the S-Lang MPI Interface Package
([slmpi](http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:mpi))
implemented by Thomas Dauser and Fritz Schwarm.

## Installation

1- Install HEASoft: download the [HEASOFT
Software](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/):

``` 
gunzip -c heasoft6.19.tar.gz | tar xf -
cd heasoft-6.19/BUILD_DIR/
./configure
make
make install
```

2- Install ISIS:

``` 
wget http://space.mit.edu/cxc/isis/install-isis.sh
setenv HEADAS /usr/local/headas/x86_64-unknown-linux-gnu-libc2.7/
sh install-isis.sh DIR
```

3- Install ISISscripts: download the [Remeis
ISISscripts](http://www.sternwarte.uni-erlangen.de/isis/)

or obtain the developing version from the ISISscripts git-repository:

``` 
git clone http://www.sternwarte.uni-erlangen.de/git.public/isisscripts 
tar xfz isisscripts.tgz
```

Add it to to your personal `~/.isisrc` startup file:

``` prolog
add_to_isis_load_path("/path/to/isisscripts/");
```

Load the ISISscripts in your S-Lang script:

``` prolog
require("isisscripts");
```

4- Install the S-Lang Interface Package
([slirp](http://space.mit.edu/cxc/slirp/)):

download [slirp-pre2.0.0-31.tar.gz](http://www.jedsoft.org/snapshots/)

Then:

``` 
gunzip -c slirp-pre2.0.0-31.tar.gz | tar xf -
./configure
make
make install
```

Add `slirp-module.so` to the `isis_load_path`.

5- Install the S-lang MPI Interface Package
([slmpi](http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:mpi)):

download
[slmpi](http://www.sternwarte.uni-erlangen.de/git.public/?p=slmpi.git)

Then:

``` 
gunzip -c slmpi.tar.gz | tar xf -
make
```

Add `slmpi-module.so` to the `isis_load_path`.

6- Install the S-lang MCMC Hammer:

Obtain the developing version from the github:

``` 
git clone https://github.com/mcfit/slmpi_emcee.git
```

Load the S-lang MCMC Hammer in your S-Lang script:

``` prolog
require("mpi_isis_emcee");
```

7- Run the example:

``` 
mpirun -np 8 isis ./example_mpi_emcee.sl
```

## Acknowledgement

If you employ **slmpi_emcee** in your work, please acknowledge the usage
by citing the following reference:

- Danehkar, A. (2025). slmpi_emcee: MPI Parallelization of the S-Lang
  MCMC Hammer. *ASP Conf. Ser.*, **538**, 405. doi:
  [10.26624/LJPV6003](https://doi.org/10.26624/LJPV6003) ads:
  [2025ASPC..538..405D](https://ui.adsabs.harvard.edu/abs/2025ASPC..538..405D/abstract).

``` bibtex
@article{Danehkar2025,
  author = {{Danehkar}, Ashkbiz,
  title = {slmpi\_emcee: MPI Parallelization of the S-Lang MCMC Hammer},
  journal = {ASP Conf. Ser.},
  volume = {538},
  pages = {405},
  year = {2025},
  doi = {10.26624/LJPV6003}
}
```

## Learn More

|                    |                                                                  |
|--------------------|------------------------------------------------------------------|
| **Repository**     | <https://github.com/mcfit/slmpi_emcee>                           |
| **Issues & Ideas** | <https://github.com/mcfit/slmpi_emcee/issues>                    |
| **DOI**            | [10.26624/LJPV6003](https://doi.org/10.26624/LJPV6003)           |
| **Archive**        | [10.5281/zenodo.4495925](https://doi.org/10.5281/zenodo.4495925) |
