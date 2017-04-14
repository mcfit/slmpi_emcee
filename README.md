# slmpi_emcee
**The MPI-based Parallelization of the S-Lang MCMC Hammer**

[slmpi_emcee](http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:emcee) is an MPI-based Parallelization of the S-Lang *Markov chain Monte Carlo (MCMC) Hammer* based on algorithm proposed by [Goodman & Weare (2010)](http://dx.doi.org/10.2140/camcos.2010.5.65), implemented in Python ([emcee](https://github.com/dfm/emcee)) by [Foreman-Mackey et al. (2013)](http://adsabs.harvard.edu/abs/2013PASP..125..306F), which has then been implemented in the [Interactive Spectral Interpretation System (ISIS)](http://space.mit.edu/cxc/isis/) in S-Lang by [Michael A. Nowak](http://space.mit.edu/home/mnowak/isis_vs_xspec/), modified by Tobias Beuchert, Lia Corrales, and Matthias Kuehnel, and extended for MPI-based Parallelization by Ashkbiz Danehkar. It is now included in the [Remeis ISISscripts](http://www.sternwarte.uni-erlangen.de/isis/). It utilizes S-Lang MPI ([slmpi](http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:mpi)) implemented by Thomas Dauser and Fritz Schwarm. 

### Installation

1- Install HEASoft: download [the HEASOFT Software](https://heasarc.nasa.gov/lheasoft/download.html)

    gunzip -c heasoft6.19.tar.gz | tar xf -
    cd heasoft-6.19/BUILD_DIR/
    ./configure
    make
    make install
    
2- Install ISIS:

    wget http://space.mit.edu/cxc/isis/install-isis.sh
    setenv HEADAS /usr/local/headas/x86_64-unknown-linux-gnu-libc2.7/
    sh install-isis.sh DIR

3- Install ISISscripts: download [the Remeis ISISscripts](http://www.sternwarte.uni-erlangen.de/isis/)

or obtain the developing version from the ISISscripts git-repository:

    git clone http://www.sternwarte.uni-erlangen.de/git.public/isisscripts 

    tar xfz isisscripts.tgz
    
Add it to to your personal `~/.isisrc` startup file

    add_to_isis_load_path("/path/to/isisscripts/");
    
Load the ISISscripts in your S-Lang script:

    require("isisscripts");

4- Install the S-Lang Interface Package ([slirp](http://space.mit.edu/cxc/slirp/)):

download [slirp-pre2.0.0-31.tar.gz](http://www.jedsoft.org/snapshots/)

Then: 

    gunzip -c slirp-pre2.0.0-31.tar.gz | tar xf -
    ./configure
    make
    make install

Add slirp-module.so to the isis_load_path.

5- Install the S-lang MPI Interface ([slmpi](http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:mpi)):

download [slmpi](http://www.sternwarte.uni-erlangen.de/git.public/?p=slmpi.git)

Then: 

    gunzip -c slmpi.tar.gz | tar xf -
    make

Add slirp-module.so to the isis_load_path.

6- Install the S-lang MCMC Hammer:

Obtain the developing version from the github:

    git clone https://github.com/mcfit/sl_emcee.git
    
Load the S-lang MCMC Hammer in your S-Lang script:

    require("mpi_isis_emcee");
