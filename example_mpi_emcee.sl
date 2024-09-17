%%! /usr/bin/env isis
% -*- slang -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_import_module_path( strcat(".", char(path_get_delimiter()), get_import_module_path()) );
require("mpi_isis_emcee");
require("isisscripts");

EM_Hash_Table_Size_Hint=500000;

% power-law parameters
variable Gamma = 2.;
% energy grid
variable Emin = 1.; % keV
variable Emax = 10.;
variable nbins = 100;
variable lo, hi;
% fake
(lo,hi) = log_grid(Emin, Emax, nbins);
fit_fun("powerlaw");
set_par("powerlaw(1).norm", 1000);
set_par("powerlaw(1).PhoIndex", Gamma);
variable flux = eval_fun_keV(lo, hi);
variable err = sqrt(flux);
()=define_counts(
  _A(hi), _A(lo), reverse(flux + grand(nbins)*err), reverse(err)
);

mpi_emcee(
  30, % number of walkers, nw, per free parameter
  1000; % number of iterations, nsim, i.e., the number of "walker"-steps
  % init_chain= "init_chain.fits", 
  output = "emcee-chain.fits", % output FITS-filename for the chain
  sim_count = 50
);

save_par("model-parameters.par");

variable nw, nfreepar, freepar, i;
variable chain_statistic, chain;

(nw, nfreepar, freepar, chain_statistic, chain) = read_chain("emcee-chain.fits");

% New Windows
variable id1 = open_plot("/xwin");
close_plot; % start with a fresh plot window
_for i (0, nw-1, 1) { % loop over all walkers
                      % and plot their path
  chain_plot(
    1, % number of the free parameter (here: Gamma)
    i, % number of the walker to plot
    nw, % totalnumber of walkers
    chain; % the chain
    oplot  % overplot the window
  ); 
};

% New Windows
variable id1 = open_plot("/xwin");
variable frac_update;
frac_update = read_chain("emcee-chain.fits"; frac_update);
xlabel("Iteration Step");
ylabel("Acceptance Rate");
plot([0:length(frac_update)-1], frac_update);

% New Windows
variable id3 = open_plot("/xwin");
chain_hist(
  0,          % number of the free parameter (here: power-law norm)
  chain;      % the chain
  fcut = .3,  % ignore the first fcut fraction of the data (see below)
  pmin = 0, pmax = 2000 % parameter range for the histogram (see below)
);
chain_hist(1, chain; fcut = .3, % parameter plotted here: photon index
  nbin = 100  % number of bins of the histogram
);

%exit();
