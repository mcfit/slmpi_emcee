===========
slmpi_emcee
===========

**The MPI-based Parallelization of the S-Lang MCMC Hammer**

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://github.com/mcfit/slmpi_emcee/blob/master/LICENSE
    :alt: GitHub license
      
.. image:: https://img.shields.io/badge/DOI-10.5281/zenodo.4495925-blue.svg
    :target: https://doi.org/10.5281/zenodo.4495925
    :alt: Zenodo

Description
===========

`slmpi_emcee <https://www.sternwarte.uni-erlangen.de/wiki/index.php/Emcee>`_ is an MPI-based Parallelization of the S-Lang *Markov chain Monte Carlo (MCMC) Hammer* based on algorithm proposed by `Goodman & Weare (2010) <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_, implemented in Python (`emcee <https://github.com/dfm/emcee>`_) by `Foreman-Mackey et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_, which has then been implemented in the `Interactive Spectral Interpretation System (ISIS) <http://space.mit.edu/cxc/isis/>`_ in S-Lang by `M. A. Nowak <http://space.mit.edu/home/mnowak/isis_vs_xspec/>`_, modified by T. Beuchert, L. Corrales, and M. Kuehnel, and extended for MPI-based Parallelization by A. Danehkar. It is now included in the `Remeis ISISscripts <http://www.sternwarte.uni-erlangen.de/isis/>`_. It utilizes the S-Lang MPI Interface Package (`slmpi <http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:mpi>`_) implemented by Thomas Dauser and Fritz Schwarm. 

Installation
============

1- Install HEASoft: download the `HEASOFT Software <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_:

.. code-block::

    gunzip -c heasoft6.19.tar.gz | tar xf -
    cd heasoft-6.19/BUILD_DIR/
    ./configure
    make
    make install
    
2- Install ISIS:

.. code-block::

    wget http://space.mit.edu/cxc/isis/install-isis.sh
    setenv HEADAS /usr/local/headas/x86_64-unknown-linux-gnu-libc2.7/
    sh install-isis.sh DIR

3- Install ISISscripts: download the `Remeis ISISscripts <http://www.sternwarte.uni-erlangen.de/isis/>`_

or obtain the developing version from the ISISscripts git-repository:

.. code-block::

    git clone http://www.sternwarte.uni-erlangen.de/git.public/isisscripts 
    tar xfz isisscripts.tgz
    
Add it to to your personal ``~/.isisrc`` startup file:

.. code-block:: prolog

    add_to_isis_load_path("/path/to/isisscripts/");
    
Load the ISISscripts in your S-Lang script:

.. code-block:: prolog

    require("isisscripts");

4- Install the S-Lang Interface Package (`slirp <http://space.mit.edu/cxc/slirp/>`_):

download `slirp-pre2.0.0-31.tar.gz <http://www.jedsoft.org/snapshots/>`_

Then:

.. code-block::

    gunzip -c slirp-pre2.0.0-31.tar.gz | tar xf -
    ./configure
    make
    make install

Add ``slirp-module.so`` to the ``isis_load_path``.

5- Install the S-lang MPI Interface Package (`slmpi <http://www.sternwarte.uni-erlangen.de/wiki/doku.php?id=isis:mpi>`_):

download `slmpi <http://www.sternwarte.uni-erlangen.de/git.public/?p=slmpi.git>`_

Then:

.. code-block::

    gunzip -c slmpi.tar.gz | tar xf -
    make

Add ``slmpi-module.so`` to the ``isis_load_path``.

6- Install the S-lang MCMC Hammer:

Obtain the developing version from the github:

.. code-block::

    git clone https://github.com/mcfit/slmpi_emcee.git
    
Load the S-lang MCMC Hammer in your S-Lang script:

.. code-block:: prolog

    require("mpi_isis_emcee");

7- Run the example:

.. code-block::

    mpirun -np 8 isis ./example_mpi_emcee.sl

Learn More
==========

==================  =============================================
**Repository**      https://github.com/mcfit/slmpi_emcee
**Issues & Ideas**  https://github.com/mcfit/slmpi_emcee/issues
**Archive**         `10.5281/zenodo.4495925 <https://doi.org/10.5281/zenodo.4495925>`_
==================  =============================================
