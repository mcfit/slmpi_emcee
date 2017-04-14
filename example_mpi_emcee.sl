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

exit();
