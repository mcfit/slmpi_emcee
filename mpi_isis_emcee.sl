
% v. 0.9.3 -- Mar. 08, 2017 (MPI-based Parallelization Version)

%%%  NOTE: Breaks backward compatability with v0.7.x !!!
%%%  v0.8.+ uses a new file format for output.  Conversion
%%%  routine included to change old fits files into new format

% A simple MCMC code for ISIS, based upon the parallel "simple
% stretch" method from "emcee: the MCMC Hammer" by Foreman-Mackey et
% al., (2012); arXiv:1202.3665

% Written by M. Nowak, mnowak@space.mit.edu, with modifications by
% Tobias Beuchert, tobias.beuchert@sternwarte.uni-erlangen.de
% Lia Corrales, lia@space.mit.edu
% Matthias Kuehnel, matthias.kuehnel@sternwarte.uni-erlangen.de, 12/09/2016 (v. 0.9.2)
% Extended for MPI-based Parallelization by Ashkbiz Danehkar, 08/03/2017 (v. 0.9.3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import("slmpi");

% Some definitions for people with older versions of S-lang

#ifnexists wherefirst
define wherefirst (a)
     {
        variable i = where (a);
        if (length(i))
          return i[0];
        else
          return NULL;
     }
#endif

#ifnexists wherefirstmax
public define wherefirstmax(a)
{
   return wherefirst (a==max(a));
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the free model parameters from an array of parameter values.
% Return a 1 if any of the steps violates a parameter boundary

#ifnexists set_pars
private define set_pars(fp,x)
{
   variable i;
   _for i (0,length(x)-1,1)
   {
      if( x[i] < get_par_info(fp[i]).min ||
          x[i] > get_par_info(fp[i]).max    )
      {
         return 1;
      }
      else
      {  
         set_par(fp[i],x[i]);
      }
   }
   return 0;
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate the log probability difference between the current walker
% and the proposal.  This presumes that P(D|M) is given by
% P(M)*exp(-chi^2/2), where P(M) is uniform within the min/max bounds
% of the parameter file, and 0 outside of that range. Note that this
% presumes that for the already known walkers, one has the fit
% statistic already calculated (sent via xstat).

#ifnexists gauss_prior
private define gauss_prior(priors)
{
    % Expected INPUTS:
    % priors    = struct{parlist, priorlist}
    % where
    %    parlist   = list of parameter ids that have priors
    %    priorlist = list of prior structures with field names "center", "sigma"
    % The priors are assumed to be gaussian and are not normalized with respect to the hard limits
    %
    % OUTPUT:
    % Sum of log-probabilities to be added to the log-likelihood of the raw fit statistic (logpn)
    % Note that chi^2 stat used in logpn is positive (negative is propogated through dlogpn), 
    %   so this needs to return a positive value.  Also do not divide by two because of treatment in logpn
    
    variable i, par, pprob, result = 0.0;
    _for i (0,length(priors.parlist)-1, 1)
    {
	par   = get_par(priors.parlist[i]);
	pprob = priors.priorlist[i];
	result += (par - pprob.mean)^2 / (pprob.sigma^2);
    }

    %message(sprintf("Adding %.2f to the chi-squared stat", result));
    return result;
}
#endif

#ifnexists logpn
private define logpn(fp,xstat,y,z)
{
    % First, deal with priors
    variable log_prior = 0.0;
    if (qualifier_exists("priors")) log_prior = gauss_prior(qualifier("priors"));

    % Now do the rest
    variable infoy, ystat=1.e32, dlogpn=-1.e32, fv_hold = Fit_Verbose, 
            i, nfp = length(fp);
    
    ifnot(set_pars(fp,y))
    {
	Fit_Verbose = -1;
	
	() = eval_counts(&infoy);  
	ystat  = infoy.statistic + log_prior;  % Log-likelihood = log-fitstat + log-priors
	%message(sprintf("Total ystat = %.2f", ystat));
	dlogpn = (xstat-ystat)/2. +  (nfp-1)*log(z);
	
	Fit_Verbose = fv_hold;
    }

    return dlogpn,ystat;
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Revamped output routine so that it appends the file, rather than
% write it from scratch each time

#ifnexists write_chain
public define write_chain( fname, nw, freepar, cdf_scl_lo, cdf_scl_hi,
                           chain, chainstat, chng, pfile )
{
   variable i, fp;
   if(qualifier_exists("create"))
   {
      fp = fits_open_file(fname,"c");

      %%% First extension: list of free parameters.
      %%%
      variable extname = "PARAMETERS";
      variable colname = ["FREE_PAR"];
      variable nrows   = length(freepar);
      variable ttype   = ["J"];
      variable tunit   = [" parameter indices"];

      fits_create_binary_table(fp, extname, nrows, colname, ttype, tunit);

      fits_update_key(fp, "PFILE", pfile, "");
      fits_update_key(fp, "MODEL", get_fit_fun, "");

      variable dlist;
      list_data(&dlist);
      if(typeof(dlist)==Undefined_Type) dlist="   Data not defined at time of file write";
      variable hist_comment = 
         ["This file contains the results of a Markov chain analysis",
          "generated by the emcee hammer subroutine",
          "This extension contains information about the model",
          "PFILE are the parameter files that started the chain",
          "MODEL is the fit function that was applied to the data",
          "FREE_PAR are the parameter indices of the free parameters",
          "The following is the list of data fit:","",
          strchop(dlist,'\n',0)
         ];
      array_map(Void_Type, &fits_write_comment, fp, hist_comment);

      if(_fits_write_col(fp,fits_get_colnum(fp,"FREE_PAR"),1,1,freepar))
         throw IOError;

      %%% Second extension: chain results for each walker/parameter
      %%%
      extname = "MCMCCHAIN";
      colname = ["FITSTAT","UPDATE"], ttype = ["D","J"], 
      tunit =  [" fit statistics"," update indicator"];
      nrows = length(chainstat);
      variable nfrpar = length(freepar);
      variable nsteps = length(chainstat[*,0]);

      _for i (0,length(freepar)-1,1)
      {
         colname = [colname,"CHAINS"+string(freepar[i])];
         ttype = [ttype,"D"];
         tunit = [tunit, " parameter values"];
      }

      fits_create_binary_table(fp, extname, nrows, colname, ttype, tunit);

      fits_update_key(fp, "NWALKERS", nw, "");
      fits_update_key(fp, "NFREEPAR", nfrpar, "");
      fits_update_key(fp, "NSTEPS", nsteps, "");
      fits_update_key(fp, "CDF_SCL_LO", cdf_scl_lo, "");
      fits_update_key(fp, "CDF_SCL_HI", cdf_scl_hi, "");

      hist_comment = 
         ["A Markov chain generated by the emcee hammer subroutine",
          "NWALKERS is the number of walkers *per* free parameter",
          "NFREEPAR is the number of free parameters",
          "NSTEPS is the number of iterations for each walker",
          "CDF_SCL_LO: The step amplitude distribution goes as 1/sqrt(z),",
          "CDF_SCL_HI: bounded by z = 1/CDF_SCL_LO -> CDF_SCL_HI","",
          "The following file columns are the results unpacked",
          "as 1D vectors of length NWALKERS*NFREEPAR*NSTEPS:","",
          "FITSTAT contains the vector of walker fit statistics",
          "UPDATE is a yes/no answer as to whether a walker updated",
          "CHAINS# are the values for the free parameter given by",
          "parameter index #" 
         ];
      array_map(Void_Type, &fits_write_comment, fp, hist_comment);

      _for i (0,nfrpar-1,1)
      {
         if(_fits_write_col(fp,
             fits_get_colnum(fp,"CHAINS"+string(freepar[i])),1,1,chain[i,*]))
            throw IOError;
      }

      _for i (0,nsteps-1,1)
      {
         if(_fits_write_col(fp,fits_get_colnum(fp,"FITSTAT"),1+i*nw*nfrpar,1,chainstat[i,*]))
            throw IOError;

         if(_fits_write_col(fp,fits_get_colnum(fp,"UPDATE"),1+i*nw*nfrpar,1,chng[i,*]))
            throw IOError;
      }

      %%% Third extension: chain statistics at each iteration
      %%%
      extname = "CHAINSTATS";
      colname = ["FRAC_UPDATE","MIN_STAT","MED_STAT","MAX_STAT"];
      nrows   = length(chainstat[*,0]);
      ttype   = ["D","D","D","D"];
      tunit   = [" fraction"," chi2"," chi2"," chi2"];

      fits_create_binary_table(fp, extname, nrows, colname, ttype, tunit);

      hist_comment = 
         ["This extension contains some useful summary information",
          "for the individual chain steps", 
          "FRAC_UPDATE is the fraction of walkers that updated",
          "MIN_STAT: minimum chi2 for a given step",
          "MED_STAT: median chi2 for a given step",
          "MAX_STAT: maximum chi2 for a given step"
         ];
      array_map(Void_Type, &fits_write_comment, fp, hist_comment);

      variable update=Double_Type[nrows], mnstat=@update, mdstat=@update, mxstat=@update;
      _for i (0,nrows-1,1)
      {
         update[i] = sum(chng[i,*])/length(chng[i,*]);
         mnstat[i] = min(chainstat[i,*]);
         mdstat[i] = median(chainstat[i,*]);
         mxstat[i] = max(chainstat[i,*]);
      }

      if(_fits_write_col(fp,fits_get_colnum(fp,"FRAC_UPDATE"),1,1,update))
         throw IOError;

      if(_fits_write_col(fp,fits_get_colnum(fp,"MIN_STAT"),1,1,mnstat))
         throw IOError;

      if(_fits_write_col(fp,fits_get_colnum(fp,"MED_STAT"),1,1,mdstat))
         throw IOError;

      if(_fits_write_col(fp,fits_get_colnum(fp,"MAX_STAT"),1,1,mxstat))
         throw IOError;
   }
   else
   {
      %%% Second extension (first unchanged): chain results
      %%%
      fp = fits_open_file(fname+"[MCMCCHAIN]","w");

      variable nsteps_prev = fits_read_key(fname+"[MCMCCHAIN]","NSTEPS");

      nfrpar = length(freepar);
      nsteps = length(chainstat[*,0]);

      nrows = nsteps_prev*nfrpar*nw;

      _for i (0,nfrpar-1,1)
      {
         if(_fits_write_col(fp,fits_get_colnum(fp,"CHAINS"+string(freepar[i])),
                                         nrows+1,1,chain[i,*]))
            throw IOError;
      }

      _for i (0,nsteps-1,1)
      {
         if(_fits_write_col(fp,fits_get_colnum(fp,"FITSTAT"),nrows+1+i*nw*nfrpar,1,chainstat[i,*]))
            throw IOError;

         if(_fits_write_col(fp,fits_get_colnum(fp,"UPDATE"),nrows+1+i*nw*nfrpar,1,chng[i,*]))
            throw IOError;
      }

      fits_update_key(fp, "NSTEPS", nsteps+nsteps_prev, "");

      %%% Third extension: chain statistics
      %%%
      fp = fits_open_file(fname+"[CHAINSTATS]","w");

      nrows   = length(chainstat[*,0]);
      update=Double_Type[nrows], mnstat=@update, mdstat=@update, mxstat=@update;
      _for i (0,nrows-1,1)
      {
         update[i] = sum(chng[i,*])/length(chng[i,*]);
         mnstat[i] = min(chainstat[i,*]);
         mdstat[i] = median(chainstat[i,*]);
         mxstat[i] = max(chainstat[i,*]);
      }

      if(_fits_write_col(fp,fits_get_colnum(fp,"FRAC_UPDATE"),nsteps_prev+1,1,update))
         throw IOError;

      if(_fits_write_col(fp,fits_get_colnum(fp,"MIN_STAT"),nsteps_prev+1,1,mnstat))
         throw IOError;

      if(_fits_write_col(fp,fits_get_colnum(fp,"MED_STAT"),nsteps_prev+1,1,mdstat))
         throw IOError;

      if(_fits_write_col(fp,fits_get_colnum(fp,"MAX_STAT"),nsteps_prev+1,1,mxstat))
         throw IOError;
  }
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private variable fs = stderr;

% A subroutine to read out a chain from a FITS file, of the format
% from the subroutine above.  If the qualifier mpi_init_walker is set, it
% will spit out an array suitable for restarting the chain.

public define mpi_old_read_chain()
{
   if(_NARGS==1)
   {
      variable fname = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  (nw, nfreepar, freepar, chain) = mpi_old_read_chain(fname);
  (nw, nfreepar, freepar, 
   cdf_scl_lo, cdf_scl_hi,
   chain, chain_statistic, 
   update_indicator, pfile, 
   walkers_init, statistic_init) = mpi_old_read_chain(fname; mpi_init_walker);
   frac_update                   = mpi_old_read_chain(fname; frac_update);

   Read the MCMC results stored in a fits file of the format produced
   by versions <0.7 of this code.  Normally, it will return just the
   minimum results required for plotting.  Qualifiers, however, can be
   used to have it return enough information to initialize a chain, or
   just return an array that has the fraction of walkers updated at
   each step.
`); %}}}
      return;
   }

   variable nw, nfreepar, freepar, cdf_scl_lo, cdf_scl_hi, pfile,
            chain, chainstat, chng, frac_update, statinit;

   (nw, nfreepar, pfile, cdf_scl_lo, cdf_scl_hi ) = 
      fits_read_key(fname+"[MCMCCHAIN]","NWALKERS","NFREEPAR","PFILE",
                                        "CDF_SCL_LO","CDF_SCL_HI");
   (freepar, chain) = fits_read_col(fname+"[MCMCCHAIN]",
                                    "FREE_PAR","CHAIN"  );
   chainstat        = fits_read_col(fname+"[CHAINSTEPS]","FITSTAT");
   chng             = fits_read_col(fname+"[CHAINSTEPS]","UPDATE");
   frac_update      = fits_read_col(fname+"[CHAINSTEPS]","FRAC_UPDATE");

   if(qualifier_exists("mpi_init_walker"))
   {
      variable iwalk   = Double_Type[nfreepar,nw*nfreepar];
      variable clength = length(chain[0,*]);
      iwalk            = chain[*,[clength-nw*nfreepar:clength-1]];
      statinit         = chainstat[-1,*];
      return nw, nfreepar, freepar, cdf_scl_lo, cdf_scl_hi, 
             chain, chainstat, chng, pfile, iwalk, statinit;
   }
   else if(qualifier_exists("frac_update"))
   {
      return frac_update;
   }
   else
   {
      return nw, nfreepar, freepar, chain;
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define mpi_old_to_new_chain()
{
   if(_NARGS==2)
   {
      variable fname_old, fname_new;
      (fname_old, fname_new) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
   mpi_old_to_new_chain(fname_old,fname_new);

   Takes the chain results stored in a version <= 0.7.x file,
   fname_old, and writes them to a version >= 0.8.0 file, fname_new.
`); %}}}
      return;
   }

   variable nw, nfreepar, freepar, cdf_scl_lo, cdf_scl_hi, chain,
    chain_statistic, update_indicator, pfile;

   (nw, nfreepar, freepar, cdf_scl_lo, cdf_scl_hi, 
    chain, chain_statistic, update_indicator, pfile, , ) = 
                      mpi_old_read_chain(fname_old; mpi_init_walker);

   write_chain(fname_new, nw, freepar, cdf_scl_lo, cdf_scl_hi, chain,
    chain_statistic, update_indicator, pfile; create); 
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private variable fs = stderr;

% A subroutine to read out a chain from a FITS file, of the format
% from the subroutine above.  If the qualifier mpi_init_walker is set, it
% will spit out an array suitable for restarting the chain.
  
public define mpi_read_chain()
{
   if(_NARGS==1)
   {
      variable fname = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  (nw, nfreepar, freepar, 
   chain_statistic, chain)       = mpi_read_chain(fname);

  (nw, nfreepar, freepar, 
   cdf_scl_lo, cdf_scl_hi,
   pfile, walkers_init, 
   statistic_init, change)       = mpi_read_chain(fname; mpi_init_walker);

   frac_update                   = mpi_read_chain(fname; frac_update);

  (frac_update, min_statistic,
   med_statistic, max_statistic) = mpi_read_chain(fname; chain_stats);

   Read the MCMC results stored in a fits file.  Normally, it will
   return just the minimum results required for plotting.  Qualifiers,
   however, can be used to have it return enough information to
   initialize a chain, return an array that has the fraction of
   walkers updated at each step, or return the fractional update array
   along with the min, median, and maximum value of the statistic at
   each step.
`); %}}}
      return;
   }
   variable i, nw, nfreepar, nsteps, cdf_scl_lo, cdf_scl_hi, pfile, 
            freepar, frac_update, chain, chain_stat, chain_chng, 
            tmp_chain, tmp_stat, tmp_chng, indces; 

   variable statinit, walkinit;

   (pfile) = fits_read_key(fname+"[PARAMETERS]","PFILE");

   (nw, nfreepar, nsteps, cdf_scl_lo, cdf_scl_hi ) = 
      fits_read_key(fname+"[MCMCCHAIN]","NWALKERS","NFREEPAR","NSTEPS",
                                        "CDF_SCL_LO","CDF_SCL_HI");

   freepar = fits_read_col(fname+"[PARAMETERS]","FREE_PAR");

   if( qualifier_exists("mpi_init_walker") )
   {
      chain = Double_Type[nfreepar,nw*nfreepar]; 
      indces = [nsteps*nw*nfreepar-nw*nfreepar:nsteps*nw*nfreepar-1];
      _for i (0,nfreepar-1,1)
      {
         tmp_chain = fits_read_col(fname+"[MCMCCHAIN]","CHAINS"+string(freepar[i]));
         chain[i,*] = tmp_chain[indces];
         () = __tmp(tmp_chain);
      }
      
      tmp_stat = fits_read_col(fname+"[MCMCCHAIN]","FITSTAT");
      chain_stat = tmp_stat[indces];
      () = __tmp(tmp_stat);

      tmp_chng = fits_read_col(fname+"[MCMCCHAIN]","UPDATE");
      chain_chng = tmp_chng[indces];
      () = __tmp(tmp_chng);

      return nw, nfreepar, freepar, cdf_scl_lo, cdf_scl_hi, pfile, 
             chain, chain_stat, chain_chng;
   }
   else if(qualifier_exists("frac_update"))
   {
      frac_update = fits_read_col(fname+"[CHAINSTATS]","FRAC_UPDATE");
      return frac_update;
   }
   else if(qualifier_exists("chain_stats"))
   {
      variable min_statistic, med_statistic, max_statistic;
      (frac_update, min_statistic,
       med_statistic, max_statistic)=fits_read_col(fname+"[CHAINSTATS]",
                                           "FRAC_UPDATE", "MIN_STAT",
                                           "MED_STAT", "MAX_STAT"     );
      return frac_update, min_statistic, med_statistic, max_statistic;
   }
   else
   {
      chain = Double_Type[nfreepar,nsteps*nw*nfreepar]; 

      _for i (0,nfreepar-1,1)
      {
         chain[i,*] = fits_read_col(fname+"[MCMCCHAIN]","CHAINS"+string(freepar[i]));
      }
      chain_stat = fits_read_col(fname+"[MCMCCHAIN]","FITSTAT");
      reshape(chain_stat,[nsteps,nw*nfreepar]);

      return nw, nfreepar, freepar, chain_stat, chain;
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inverse of the Cumulative Distribution Function for 1/sqrt(z),
% if 1/a <= z <= b, to use as a random number generator for a
% 1/sqrt(z) probability distribution.  Parameter steps will have an
% amplitude governed by random draws from this distribution.

#ifnexists icdf
private define icdf(x,a,b)
{
   variable c1 = 1/(sqrt(a*b)-1);
   variable c2 = 1/c1^2/a;
   return c2*(x+c1)^2;
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return an array with the indices of the free parameters

#ifnexists free_par
private define free_par()
{
   variable i, free_par = Integer_Type[0];
   _for i (1,get_num_pars,1)
   {
      if(get_par_info(i).freeze == 0 &&
         get_par_info(i).tie == NULL)
      {
         free_par = [free_par,i];
      }
   }
   return free_par;
}
#endif
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize nw walkers *per free parameter*, whose indices are in fp.
% For now, we will initialize uniformly, or as a gaussian, over the
% range given for the min/max values of the free parameters, or over
% the min/max values for an input parameter file, both modulo a scale
% factor (==1 for full min/max range, uniform draw)

private define mpi_init_walker(fp,nw,proc_rank,proc_total)
{
   variable iret=0; 
   variable npfile=1, fname = "./tmp_mcmc_pfile_"+string(getpid)+".par";
   %variable npfile=1, fname = "./tmp_mcmc_pfile_"+string(proc_rank)+".par";
   variable fname_len=int(1);
   if(proc_rank==0)
   {
      try{
         save_par(fname);
      }
      catch AnyError: {
         () = printf("\n\n Failed saving initial parameters \n\n");
          plot_pause;
      }

      if(qualifier_exists("pfile"))
      {
         variable pfile_init, pfiles = strchop(qualifier("pfile"),';',0);
         npfile = length(pfiles);

         try{
            foreach pfile_init (pfiles)
            {
               load_par(pfile_init);
            }
         }
         catch AnyError: {
            () = printf("\n\n Cannot read parameter file: %s \n",pfile_init);
            () = printf(" Retaining original parameter file for initialization. \n\n");
            npfile = 1;
            pfiles = [fname];
            load_par(fname);
         }
      }
      else
      {
         pfiles = [fname];
      }
      fname_len=strlen(fname);
      fname_len=int(fname_len);
   }
   iret=rcl_mpi_barrier();
   iret=rcl_mpi_bcast_int([npfile],1);
   %iret=rcl_mpi_bcast_int([fname_len],1);
   %iret=rcl_mpi_bcast_char(fname,26);
      
   variable sc, sig;
   if(qualifier_exists("gaussian"))
   {
      sc = qualifier("scale",1./3.);
   }
   else
   {
      sc = qualifier("scale",1.);
   }

   variable fv_hold = Fit_Verbose;
   Fit_Verbose = -1;

   % The walkers are 2D arrays, with the first index covering the free
   % parameters, and the second index being the specific instances
   % (i.e., the length of the 2nd index is the number of walkers).

   variable i, j, k, info, nfp = length(fp),
            xp    = Double_Type[nfp, npfile], xl=@xp, xh=@xp, 
            xstrt = Double_Type[nfp,nw*nfp],
            xstat = Double_Type[nw*nfp];
   
   iret=rcl_mpi_barrier();
   if(proc_rank==0)
   {
      _for k (0,npfile-1,1)
      {   
         load_par(pfiles[k]);
         () = eval_counts;

         _for i (0,nfp-1,1)
         {
            xp[i,k]      = get_par_info(fp[i]).value;
            xl[i,k]      = get_par_info(fp[i]).min;
            xh[i,k]      = get_par_info(fp[i]).max;
         }
      }
   }
   iret=rcl_mpi_barrier();
   iret=rcl_mpi_bcast_double(xp,nfp*npfile);
   iret=rcl_mpi_bcast_double(xl,nfp*npfile);
   iret=rcl_mpi_bcast_double(xh,nfp*npfile);
   iret=rcl_mpi_barrier();
   % Reload the original parameter file, which should have the widest
   % allowable min/max values (enough to encompass any ranges present
   % in the input parameter files)

   if(qualifier_exists("pfile"))
   {
      try{ 
         load_par(fname);
      }
      catch AnyError: {
         () = printf("\n\n Failed reloading original  parameters \n\n");
         plot_pause;
      }
      if(eval_counts)
      {
         () = printf("\n\n Failed evaluating original  parameters \n\n");
          plot_pause;
      }
   }
   variable nw_nfp_section= round(nw*nfp/proc_total);
   variable nw_nfp_min=0;
   variable nw_nfp_max=0;
   if (proc_rank == proc_total-1) 
   {
      nw_nfp_min=int(proc_rank*nw_nfp_section);
      nw_nfp_max=int(nw*nfp-1);
   } 
   else
   {
      nw_nfp_min=int(proc_rank*nw_nfp_section);
      nw_nfp_max=int(((proc_rank+1)*nw_nfp_section)-1);
   }
   %_for j (0,nw*nfp-1,1)
   _for j (nw_nfp_min,nw_nfp_max,1)
   {
      k = j*npfile/(nw*nfp);

      _for i (0,nfp-1,1)
      {
         if(qualifier_exists("gaussian"))
         {
            sig = (@(qualifier("grnd", &grand)))();
            if(sig < 0)
            {
               xstrt[i,j] = xp[i,k] + sig*sc*(xp[i,k]-xl[i,k]);
            }
            else
            {
               xstrt[i,j] = xp[i,k] + sig*sc*(xh[i,k]-xp[i,k]);
            }
         }
         else
         {
            xstrt[i,j] = (1-sc)*xp[i,k] + 
              sc*(xl[i,k]+(xh[i,k]-xl[i,k])*(@(qualifier("urnd", &urand)))());
         }

         if(xstrt[i,j]<xl[i,k]) xstrt[i,j] = xl[i,k];
         if(xstrt[i,j]>xh[i,k]) xstrt[i,j] = xh[i,k];
         set_par(fp[i],xstrt[i,j]);
      }

      () = eval_counts(&info);
      variable stat = info.statistic;
      if (qualifier_exists("priors"))
      {
         stat += gauss_prior(qualifier("priors"));
      }
      xstat[j] = stat;
   }
   iret=rcl_mpi_barrier();
   variable lentag=0;
   if(proc_rank==0)
   {
      _for i (1,proc_total-1,1)
      {
         if (i == proc_total-1) 
         {
            nw_nfp_min=int(i*nw_nfp_section);
            nw_nfp_max=int(nw*nfp-1);
         } 
         else
         {
            nw_nfp_min=int(i*nw_nfp_section);
            nw_nfp_max=int(((i+1)*nw_nfp_section)-1);
         }
         variable xstat_recv=xstat[[nw_nfp_min:nw_nfp_max]];
         iret=rcl_mpi_org_recv_double(xstat_recv, nw_nfp_max-nw_nfp_min+1, i, lentag);
         xstat[[nw_nfp_min:nw_nfp_max]]=xstat_recv[[0:nw_nfp_max-nw_nfp_min]];
      }
   }
   else
   {
      variable xstat_send=xstat[[nw_nfp_min:nw_nfp_max]];
      iret=rcl_mpi_org_send_double(xstat_send, nw_nfp_max-nw_nfp_min+1, 0, lentag);
   }
   iret=rcl_mpi_barrier();
   iret=rcl_mpi_bcast_double(xstat[*],nw*nfp);
   iret=rcl_mpi_barrier();
   _for k (0,nfp-1,1)
   {
      if(proc_rank==0)
      {
         _for i (1,proc_total-1,1)
         {
            if (i == proc_total-1) 
            {
               nw_nfp_min=int(i*nw_nfp_section);
               nw_nfp_max=int(nw*nfp-1);
            } 
            else
            {
               nw_nfp_min=int(i*nw_nfp_section);
               nw_nfp_max=int(((i+1)*nw_nfp_section)-1);
            }
            variable xstrt_recv=xstrt[k,[nw_nfp_min:nw_nfp_max]];
            iret=rcl_mpi_org_recv_double(xstrt_recv, nw_nfp_max-nw_nfp_min+1, i, lentag);
            xstrt[k,[nw_nfp_min:nw_nfp_max]]=xstrt_recv[[0:nw_nfp_max-nw_nfp_min]];
         }
      }
      else
      {
         variable xstrt_send=xstrt[k,[nw_nfp_min:nw_nfp_max]];
         iret=rcl_mpi_org_send_double(xstrt_send, nw_nfp_max-nw_nfp_min+1, 0, lentag);
      }
      iret=rcl_mpi_barrier();
   }
   iret=rcl_mpi_barrier();
   iret=rcl_mpi_bcast_double(xstrt[*,*],nw*nfp*nfp);
   iret=rcl_mpi_barrier();
   Fit_Verbose = fv_hold;
   return xstrt, xstat;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The function to create a proposition (i.e., the new trial
% parameters), and test whether or not that proposition is accepted.

private define mpi_update_walker(cdf_scl_lo,cdf_scl_hi,rnd,fp,xa,xstat,xb)
{
   % xa is a vector of parameters for a specific walker to be updated,
   % while xstat is the fit statistic of that walker. xb is an array
   % of walker parameters (a second, independent group) from which an
   % update is drawn.  The first index of xb spans the parameter,
   % while the second spans specific walker instances.

   variable lb = length(xb[0,*]);                   % No. of walkers
   variable xuse = xb[*,int(rnd[0]*lb)];            % Choose one
   variable z = icdf(rnd[1],cdf_scl_lo,cdf_scl_hi); % Random step size

   % New trial parameters
   xuse += z*(xa-xuse); 
   
   %variable priors = qualifier("priors",0);
   variable dlogpn, ystat;
   (dlogpn,ystat) = logpn(fp,xstat,xuse,z;; __qualifiers);

   variable ptest = log(rnd[2]);
   if(ptest <= dlogpn)
   {
      return xuse, ystat, 1;             % 1 = changed value
   }
   else
   {
      return xa, xstat, 0;               % 0 = same value
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This variable will be used to receive their results.

private variable MPI_Task_List = 
   struct{ index, length, rands, free_par, walkers, wstats, change, 
           proposals, output };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The function to do the emcee hammer.  It takes as input the number
% of walkers per free parameter, and the number of simulations per set
% of walkers.  Optionally, it can take a qualifier to read the
% parameter min/max values from other file(s).  Or it can be
% initialized from an existing chain.

public define mpi_emcee()
{
   variable proc_rank=rcl_mpi_init();
   variable proc_total=rcl_mpi_numtasks(); 
   if(_NARGS==2)
   {
      variable nw, nsim;
      (nw, nsim) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
   mpi_emcee( nw, nsim; pfile="file1;file2;...;filen", 
                    scale=#, gaussian, 
                    cdf_scl_lo=#, cdf_scl_hi=#, 
                    init_chain=file, 
                    output=file, sim_count=#, 
                    nocopy, priors, urnd=&urand, grnd=&grand );

  Create a set of MCMC chains for an ensemble of walkers (i.e.,
  propositions) consisting of a number of walkers equal to nw X the
  number of free parameters.  This ensemble will be evolved for nsim
  iterations.  All results will be stored in an output file.

  The initial distribution of walkers is taken from the current
  parameter values, or from the parameter values in the files listed
  in pfile (file names, including full paths if not in the local
  directory, separated by semi-colons).  If the gaussian qualifier is
  not set, each walker has its parameters uniformly distributed from
  value-scale*(value-min) to value+scale*(max-value). (Default
  scale=1.)

  Otherwise, the walkers have initial parameters gaussianly
  distributed about value with sigma = scale*(max-value) or sigma =
  scale*(value-min) (i.e., normalized to the possibly asymmetric
  min/max values).  (Default scale=1/3.)

  Alternatively, the walkers can be initialized from the last
  iteration of an existing chain stored in a file, init_chain.

  The walkers are split into groups (0) and (1), with the update
  proposition (i.e., parameter vector) for the group (0) j^th walker
  being:
            X(0)'_j = X(1)_k + z*(X(0)_j-X(1)_k)
  X(1)_k is randomly drawn from group (1) and z is randomly
  distributed as 1/sqrt(z) between values of 1/cdf_scl_lo and
  cdf_scl_hi (default for both=2). Group (1) is then subsequently
  updated in a similar manner.

  The output chain is placed in the file specified by the output
  qualifier (default name is mcmc_results_processid.fits).
  Intermediate results are written to this file every sim_count
  (default=50) iterations.  If the chains were started from a previous
  chain file, the final output file will contain the initial chains as
  well unless the nocopy qualifier has been set.

  The optional qualifiers grnd and urnd can be set to user defined
  functions to, e.g., allow random numbers drawn from a server or file
  in case of a manual multi-job run (among several machines).

  PRIORS
  One can apply gaussian priors to particular parameters by
  constructing a structure with the following fields
  priors.parlist -- an array of parameter ids
  priors.priorlist -- an array of structures describing each prior. A
     prior structure needs two fields: 'mean' and 'sigma'

  For example:
  myprior = struct{mean=0.1, sigma=0.2}
  mpi_emcee(nw, nsim; priors=struct{parlist=[1], priorlist=[myprior]})
`); %}}}
      return;
   }
   % Store initial parameters so that we can restore them at the end
   variable pars = get_params;

   variable isc=0, sim_count_array, sim_count = int(qualifier("sim_count",50));
   if(sim_count <= 0) sim_count=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   sim_count_array=[0:nsim+sim_count:sim_count];
   sim_count_array=[sim_count-1:nsim-1+sim_count:sim_count];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   variable ochain, ostat, ochng;

   variable screwup    = qualifier("infile","NONE");
            screwup    = qualifier("in_file",screwup);
            screwup    = qualifier("start",screwup);
            screwup    = qualifier("start_file",screwup);
            screwup    = qualifier("initfile",screwup);
            screwup    = qualifier("init_file",screwup);
            screwup    = qualifier("input",screwup);
            screwup    = qualifier("initchain",screwup);
   variable init_chain = qualifier("init_chain",screwup);
   
   if ( screwup != "NONE" ) 
      if (proc_rank==0) ()=printf("\n *** Using %s for chain initialization  *** \n\n",init_chain);
      
            screwup    = qualifier("outfile","mcmc_results_");
            screwup    = qualifier("out_file",screwup);
            screwup    = qualifier("stop",screwup);
            screwup    = qualifier("stop_file",screwup);
            screwup    = qualifier("outchain",screwup);
            screwup    = qualifier("out_chain",screwup);
   variable outfile    = qualifier("output",screwup);

   if ( outfile == "mcmc_results_" ) outfile += "_"+string(getpid)+".fits";
   if (proc_rank==0) ()=printf("\n *** Using %s for output file *** \n\n",outfile);

   variable pfile      = qualifier("pfile","NONE");
   variable cdf_scl_lo = qualifier("cdf_scl_lo",2.);
   variable cdf_scl_hi = qualifier("cdf_scl_hi",2.);

   nw   = int(nw);
   nsim = int(nsim);

   variable xwalk, xstat, xwalk_a, xwalk_b, freepar = free_par();
   
   if(init_chain=="NONE")
   {
      (xwalk, xstat) = mpi_init_walker(freepar,nw,proc_rank,proc_total;;__qualifiers);
   }
   else
   {
      variable nwc=0, nfpc=0, fpc;
      try
      {
         (nwc, nfpc, fpc, cdf_scl_lo, cdf_scl_hi, pfile,
          xwalk, xstat, ochng)=mpi_read_chain(init_chain;mpi_init_walker);
         if (proc_rank==0)
         {
            if(qualifier_exists("nocopy"))
            {
               variable wstat=@xstat, wchng=@ochng;
               reshape(wstat,[1,length(wstat)]);
               reshape(wchng,[1,length(wchng)]);
               write_chain(outfile, nw, freepar, cdf_scl_lo, cdf_scl_hi,
                           xwalk, wstat, wchng, pfile;create);
            }
            else
            {
               () = system("cp "+init_chain+" "+outfile);
            }
         }
      }
      catch AnyError:
      {
         if (proc_rank==0) ()=printf("\n Failed to read or create chain file \n");
      }
      % check, if present fit corresponds to the one, one wants to append to (init_chain)
      if(nwc!=nw || nfpc!=length(freepar) || max(fpc-freepar)!=0)
      {
         if (proc_rank==0) ()=printf("\n Chain file inconsistent, initializing from data \n");
         (xwalk, xstat) = mpi_init_walker(freepar,nw,proc_rank,proc_total;;__qualifiers);
         init_chain="NONE";
         pfile="NONE";
      }
   }

   if(1/cdf_scl_lo >= cdf_scl_hi || cdf_scl_lo <= 0 || cdf_scl_hi <= 0)
   {
      if (proc_rank==0) ()=printf("\n1/CDF_SCL_LO must be < CDF_SCL_HI and > 0.\n");
      return;
   }

   variable iwlk_a, iwlk_b, la, lb, im, i, j, k, s, 
            nwf       = nw*num_free_params,
            xchng     = @xstat,
            xout      = Double_Type[num_free_params,sim_count*nwf],
            xout_stat = Double_Type[sim_count,nwf],
            xout_chng = Double_Type[sim_count,nwf]; 
   
   iwlk_a = [0:nwf-1:2];
   iwlk_b = [1:nwf-1:2];

   la = length(iwlk_a);
   lb = length(iwlk_b);

   MPI_Task_List.free_par  = freepar;
   MPI_Task_List.output    = Array_Type[max([la,lb])];
   MPI_Task_List.wstats    = Double_Type[max([la,lb])];
   MPI_Task_List.change    = Long_Type[max([la,lb])];

   variable array_xwalk      = Array_Type[max([la,lb])];
   variable array_xwalk_stat = Double_Type[max([la,lb])];

   % Seed all the required random numbers at the start, so we don't
   % have issues with stuff when it comes time to parallelize
   variable rands_a = Array_Type[1];%[nsim*la];
   variable rands_b = Array_Type[1];%[nsim*lb];
   variable rnd_var = Double_Type[3];

   variable rands_a_new = Double_Type[nsim*la,3];
   variable rands_b_new = Double_Type[nsim*lb,3];
   
   %
   if (proc_rank==0)
   {
      _for i (0,nsim*la-1,1) 
      { 
         rands_a[0] = (@(qualifier("urnd", &urand)))(3); 
         rands_a_new[i,0]=rands_a[0][0];
         rands_a_new[i,1]=rands_a[0][1];
         rands_a_new[i,2]=rands_a[0][2];
      }
      _for i (0,nsim*lb-1,1) 
      { 
         rands_b[0] = (@(qualifier("urnd", &urand)))(3); 
         rands_b_new[i,0]=rands_b[0][0];
         rands_b_new[i,1]=rands_b[0][1];
         rands_b_new[i,2]=rands_b[0][2];
      }
   }
   variable iret=0;  
   iret=rcl_mpi_barrier();
   iret=rcl_mpi_bcast_double(rands_a_new,3*nsim*la);
   iret=rcl_mpi_bcast_double(rands_b_new,3*nsim*lb);
   iret=rcl_mpi_barrier();
   
   variable irnd;
   
   
   variable lentag=0;
    
   variable la_section= int(ceil(la/proc_total));
   variable la_min=0;
   variable la_max=0;
   variable iwlk_a_list1= Integer_Type[la_section];
   variable iwlk_a_list2= Integer_Type[la_section];
   
   variable xwalk_recv_a=Double_Type[num_free_params,la_section];
   variable xwalk_recv_a_len=0;
   variable xstat_recv_a= Double_Type[la_section];
   variable xstat_recv_a_len=0;
   variable xchng_recv_a= Double_Type[la_section];
   variable xchng_recv_a_len=0;

   variable xwalk_send_a=Double_Type[num_free_params,la_section];
   variable xwalk_send_a_len=0;
   variable xstat_send_a= Double_Type[la_section];
   variable xstat_send_a_len=0;
   variable xchng_send_a= Double_Type[la_section];
   variable xchng_send_a_len=0;  
     
   variable lb_section= int(ceil(lb/proc_total));
   variable lb_min=0;
   variable lb_max=0;
   variable iwlk_b_list1= Integer_Type[lb_section];
   variable iwlk_b_list2= Integer_Type[lb_section];
   
   variable xwalk_recv_b=Double_Type[num_free_params,lb_section];
   variable xwalk_recv_b_len=0;
   variable xstat_recv_b= Double_Type[lb_section];
   variable xstat_recv_b_len=0;
   variable xchng_recv_b= Double_Type[lb_section];
   variable xchng_recv_b_len=0;
   
   variable xwalk_send_b=Double_Type[num_free_params,lb_section];
   variable xwalk_send_b_len=0;
   variable xstat_send_b= Double_Type[lb_section];
   variable xstat_send_b_len=0;
   variable xchng_send_b= Double_Type[lb_section];
   variable xchng_send_b_len=0;  
   variable xdata_send_b=Double_Type[length(xwalk)+length(xstat)+length(xchng)];
   variable xdata_recv_b=Double_Type[length(xwalk)+length(xstat)+length(xchng)];
   
   variable xwalk_len=0;
   variable xstat_len=0;
   variable xchng_len=0;
       
   iret=rcl_mpi_barrier();
   iret=rcl_mpi_bcast_double(xwalk,length(xwalk));
   iret=rcl_mpi_bcast_double(xstat,length(xstat));
   iret=rcl_mpi_bcast_double(xchng,length(xchng));
   la_section= ceil(la/proc_total);
   lb_section= ceil(lb/proc_total);
   iret=rcl_mpi_barrier();
   _for i (0,nsim-1,1)
   {
      % Update first half of walkers ...
      irnd = i*la+[0:la-1];
      
      if (proc_rank == proc_total-1) 
      {
         la_min=int(proc_rank*la_section);
         la_max=int(la-1);
      } 
      else
      {
         la_min=int(proc_rank*la_section);
         la_max=int(((proc_rank+1)*la_section)-1);
      }
      %_for j (0,la-1,1)
      _for j (la_min,la_max,1)
      {
         array_xwalk[j]      = xwalk[*,iwlk_a[j]];
         array_xwalk_stat[j] = xstat[iwlk_a[j]];
         if(i==0 && j==0) () = printf("\n Start Run \n");
         rnd_var=[rands_a_new[irnd[j],0],rands_a_new[irnd[j],1],rands_a_new[irnd[j],2]];
         (MPI_Task_List.output[j],MPI_Task_List.wstats[j],
          MPI_Task_List.change[j]) = mpi_update_walker(cdf_scl_lo,cdf_scl_hi,
                                               rnd_var, 
                                               freepar,
                                               array_xwalk[j],
                                               array_xwalk_stat[j],
                                               xwalk[*,iwlk_b]);
      }

      %_for j (0,la-1,1)
      _for j (la_min,la_max,1)
      {
         xwalk[*,iwlk_a[j]] = MPI_Task_List.output[j];
         xstat[iwlk_a[j]]   = MPI_Task_List.wstats[j];
         xchng[iwlk_a[j]]   = MPI_Task_List.change[j];
      }
      
      % Update second half of walkers ...
      irnd = i*lb+[0:lb-1];
      
      if (proc_rank == proc_total-1) 
      {
         lb_min=int(proc_rank*lb_section);
         lb_max=int(lb-1);
      } 
      else
      {
         lb_min=int(proc_rank*lb_section);
         lb_max=int(((proc_rank+1)*lb_section)-1);
      }
      
      %_for j (0,lb-1,1)
      _for j (lb_min,lb_max,1)
      {
         array_xwalk[j]      = xwalk[*,iwlk_b[j]];
         array_xwalk_stat[j] = xstat[iwlk_b[j]];
         rnd_var=[rands_b_new[irnd[j],0],rands_b_new[irnd[j],1],rands_b_new[irnd[j],2]];
         (MPI_Task_List.output[j],MPI_Task_List.wstats[j],
          MPI_Task_List.change[j]) = mpi_update_walker(cdf_scl_lo,cdf_scl_hi,
                                               rnd_var, 
                                               freepar,
                                               array_xwalk[j],
                                               array_xwalk_stat[j],
                                               xwalk[*,iwlk_a]);
      }

      %_for j (0,lb-1,1)
      _for j (lb_min,lb_max,1)
      {
         xwalk[*,iwlk_b[j]] = MPI_Task_List.output[j];
         xstat[iwlk_b[j]]   = MPI_Task_List.wstats[j];
         xchng[iwlk_b[j]]   = MPI_Task_List.change[j];
      }
      
      lentag=0;
      iret=rcl_mpi_barrier();
      
      if(proc_rank==0)
      {
         _for k (1,proc_total-1,1)
         {
            if (k == proc_total-1) 
            {
               la_min=int(k*la_section);
               la_max=int(la-1);
            } 
            else
            {
               la_min=int(k*la_section);
               la_max=int(((k+1)*la_section)-1);
            }
            iwlk_a_list1=int(iwlk_a[[la_min:la_max]]);
            xwalk_recv_a_len=int(length(xwalk[*,iwlk_a_list1]));
            xstat_recv_a_len=int(length(xstat[iwlk_a_list1]));
            xchng_recv_a_len=int(length(xchng[iwlk_a_list1]));
            
            if (k == proc_total-1) 
            {
               lb_min=int(k*lb_section);
               lb_max=int(lb-1);
            } 
            else
            {
               lb_min=int(k*lb_section);
               lb_max=int(((k+1)*lb_section)-1);
            }
            iwlk_b_list1=int(iwlk_b[[lb_min:lb_max]]);
            xwalk_recv_b_len=int(length(xwalk[*,iwlk_b_list1]));
            xstat_recv_b_len=int(length(xstat[iwlk_b_list1]));
            xchng_recv_b_len=int(length(xchng[iwlk_b_list1]));
            
            variable xdata_recv_a= Double_Type[xwalk_recv_a_len+3*xstat_recv_a_len+xwalk_recv_b_len+3*xstat_recv_b_len];
            
            iret=rcl_mpi_org_recv_double(xdata_recv_a, xwalk_recv_a_len+3*xstat_recv_a_len + xwalk_recv_b_len+3*xstat_recv_b_len, k, 0);
            
            iwlk_a_list1=int(xdata_recv_a[[0:xstat_recv_a_len-1]]);
            xwalk[*,iwlk_a_list1]=xdata_recv_a[[xstat_recv_a_len+0:xwalk_recv_a_len+xstat_recv_a_len-1]];
            xstat[iwlk_a_list1]=xdata_recv_a[[xwalk_recv_a_len+xstat_recv_a_len:xwalk_recv_a_len+2*xstat_recv_a_len-1]];
            xchng[iwlk_a_list1]=xdata_recv_a[[xwalk_recv_a_len+2*xstat_recv_a_len:xwalk_recv_a_len+3*xstat_recv_a_len-1]];

            iwlk_b_list1=int(xdata_recv_a[[xwalk_recv_a_len+3*xstat_recv_a_len + 0:xwalk_recv_a_len+3*xstat_recv_a_len + xstat_recv_b_len-1]]);
            xwalk[*,iwlk_b_list1]=xdata_recv_a[[xwalk_recv_a_len+3*xstat_recv_a_len + xstat_recv_b_len+0:xwalk_recv_a_len+3*xstat_recv_a_len + xwalk_recv_b_len+xstat_recv_b_len-1]];
            xstat[iwlk_b_list1]=xdata_recv_a[[xwalk_recv_a_len+3*xstat_recv_a_len + xwalk_recv_b_len+xstat_recv_b_len:xwalk_recv_a_len+3*xstat_recv_a_len + xwalk_recv_b_len+2*xstat_recv_b_len-1]];
            xchng[iwlk_b_list1]=xdata_recv_a[[xwalk_recv_a_len+3*xstat_recv_a_len + xwalk_recv_b_len+2*xstat_recv_b_len:xwalk_recv_a_len+3*xstat_recv_a_len + xwalk_recv_b_len+3*xstat_recv_b_len-1]];
         }
         _for k (1,proc_total-1,1)
         {
            xdata_send_b=Double_Type[length(xwalk)+length(xstat)+length(xchng)];
            xdata_send_b[[0:length(xwalk)-1]]=xwalk[*,*];
            xdata_send_b[[length(xwalk):length(xwalk)+length(xstat)-1]]=xstat[*];
            xdata_send_b[[length(xwalk)+length(xstat):length(xwalk)+length(xstat)+length(xchng)-1]]=xchng[*];
            iret=rcl_mpi_org_send_double(xdata_send_b,length(xwalk)+length(xstat)+length(xchng), k, 0); 
         }
      }
      else
      {
         iwlk_a_list2=int(iwlk_a[[la_min:la_max]]);
         xwalk_send_a_len=int(length(xwalk[*,iwlk_a_list2]));
         xstat_send_a_len=int(length(xstat[iwlk_a_list2]));
         xchng_send_a_len=int(length(xchng[iwlk_a_list2]));
         
         iwlk_b_list2=int(iwlk_b[[la_min:la_max]]);
         xwalk_send_b_len=int(length(xwalk[*,iwlk_b_list2]));
         xstat_send_b_len=int(length(xstat[iwlk_b_list2]));
         xchng_send_b_len=int(length(xchng[iwlk_b_list2]));
         
         variable xdata_send_a= Double_Type[xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+3*xstat_send_b_len];

         xdata_send_a[[0:xstat_send_a_len-1]]=double(iwlk_a_list2);
         xdata_send_a[[xstat_send_a_len+0:xwalk_send_a_len+xstat_send_a_len-1]]=xwalk[*,iwlk_a_list2];
         xdata_send_a[[xwalk_send_a_len+xstat_send_a_len:xwalk_send_a_len+2*xstat_send_a_len-1]]=xstat[iwlk_a_list2];
         xdata_send_a[[xwalk_send_a_len+2*xstat_send_a_len:xwalk_send_a_len+3*xstat_send_a_len-1]]=xchng[iwlk_a_list2];
         
         xdata_send_a[[xwalk_send_a_len+3*xstat_send_a_len + 0:xwalk_send_a_len+3*xstat_send_a_len + xstat_send_b_len-1]]=double(iwlk_b_list2);
         xdata_send_a[[xwalk_send_a_len+3*xstat_send_a_len + xstat_send_b_len+0:xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+xstat_send_b_len-1]]=xwalk[*,iwlk_b_list2];
         xdata_send_a[[xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+xstat_send_b_len:xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+2*xstat_send_b_len-1]]=xstat[iwlk_b_list2];
         xdata_send_a[[xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+2*xstat_send_b_len:xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+3*xstat_send_b_len-1]]=xchng[iwlk_b_list2];
         
         iret=rcl_mpi_org_send_double(xdata_send_a, xwalk_send_a_len+3*xstat_send_a_len + xwalk_send_b_len+3*xstat_send_b_len, 0, 0); 
         
         xdata_recv_b=Double_Type[length(xwalk)+length(xstat)+length(xchng)];
         iret=rcl_mpi_org_recv_double(xdata_recv_b, length(xwalk)+length(xstat)+length(xchng), 0, 0);
         xwalk[*,*]=xdata_recv_b[[0:length(xwalk)-1]];
         xstat[*]=xdata_recv_b[[length(xwalk):length(xwalk)+length(xstat)-1]];
         xchng[*]=xdata_recv_b[[length(xwalk)+length(xstat):length(xwalk)+length(xstat)+length(xchng)-1]];
      }
      
      iret=rcl_mpi_barrier();
      
      %xwalk_len=length(xwalk[*,*]);
      %xstat_len=length(xstat[*]);
      %xchng_len=length(xchng[*]);
      %iret=rcl_mpi_bcast_double(xwalk[*,*],xwalk_len);
      %iret=rcl_mpi_bcast_double(xstat[*],xstat_len);
      %iret=rcl_mpi_bcast_double(xchng[*],xchng_len);
      %iret=rcl_mpi_barrier();
      
      if (proc_rank==0)
      {
         im = i mod sim_count;
         
         xout[*,im*nwf+[0:nwf-1:2]] = xwalk[*,iwlk_a];
         xout[*,im*nwf+[1:nwf-1:2]] = xwalk[*,iwlk_b];

         xout_stat[im,[0:nwf-1:2]] = xstat[iwlk_a];
         xout_stat[im,[1:nwf-1:2]] = xstat[iwlk_b];

         xout_chng[im,[0:nwf-1:2]] = xchng[iwlk_a];
         xout_chng[im,[1:nwf-1:2]] = xchng[iwlk_b];
         
         if(i==0){ () = printf("\n Sim loop: \n"); }
         if( i==sim_count_array[isc] || i==nsim-1 )
         {
            if(isc==0)
            {
               if ((init_chain=="NONE") or (qualifier_exists("nocopy")))
               {
                  write_chain(outfile, nw, freepar, cdf_scl_lo, cdf_scl_hi,
                              xout, xout_stat, xout_chng, pfile; create);
               }
               else
               {
                  write_chain(outfile, nw, freepar, cdf_scl_lo, cdf_scl_hi,
                              xout, xout_stat, xout_chng, pfile);
               }
            }
            else
            {
               write_chain(outfile, nw, freepar, cdf_scl_lo, cdf_scl_hi,
                                xout, xout_stat, xout_chng, pfile);
            }
         
            () = printf("%d ... \n",i+1);
            isc++;
         }
      }
   }
   iret=rcl_mpi_barrier();
   if (proc_rank==0)
   {
      if(sim_count > 0) () = printf("\n");
   }
   set_params(pars);
   
   rcl_mpi_finalize();
   
}

%%%%%%%%%%%%%%%%%%%%%%%

#ifnexists unity
private define unity(x)
{
   return x;
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifnexists chain_hist
public define chain_hist()
{
   if(_NARGS==2)
   {
      variable i, y;
      (i, y) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  chain_hist(i, chain; fcut=#, lcut=#, nbin=#, color=#, linestyle=#,
                       pmin=#, pmax=#, 
                       xrange={# or NULL,# or NULL}, 
                       yrange={# or NULL,# or NULL},
                       transform=&func, 
                       function=&func, func_vars=[Array],
                       clevel=#, print, return_map, return_median,
                       xlabel=String, ylabel=String, oplot); 

  Plot a normalized histogram of the MCMC chain for the i^th free
  parameter, between values pmin and pmax (defaults to the min/max
  values in the chain; normalization, however, is always for the full
  range of the chain subsample). The plot axes have the ranges given
  in the lists, xrange and yrange (NULL==autoscale). Throw out the
  first fcut fraction of the data (default=0.2) and the last lcut
  fraction of the data (default=0), and use nbin bins (default=50).

  The x-axis variable can be transformed using any monotonic function
  passed as a reference via the transform qualifier (e.g.,
  transform=&log10).  pmin/pmax will then be in reference to this
  transform variable.

  Alternatively, if a multivariate function is passed via the function
  qualifier, with func_vars being an array indicating which of the i
  free parameters serve as variables.  These choices supercede both
  the parameter index input, i, and any transform input.  The function
  should expect input of the form x[i,k], where the i index is for the
  free parameters of the fit, and the k index is a specific instance
  from the chain.

  Plot using pgplot color given by qualifier (default=1, black) and
  linestyle given by qualifier (default=1, solid line).  The plot will
  be an overplot if the oplot qualifier is set. The plot function
  presumes that the appropriate fit model has been loaded.  Custom X/Y
  labels will be applied if the x/ylabel qualifiers have been set.
  The print qualifier has the lower and upper clevel (default=0.9)
  confidence levels printed out. The qualifier return_map, returns the
  MAP (maximum a posteriori) value, i.e., the peak of the
  distribution, return_median the median of the posterior
  distribution.

`); %}}}
      return;
   }
   variable xr = qualifier("xrange",{NULL,NULL});
   variable yr = qualifier("yrange",{NULL,NULL});
   xrange(xr[0],xr[1]); yrange(yr[0],yr[1]);
   connect_points(-1);
   variable fp = free_par;
   variable xlab = qualifier("xlabel",NULL);
   variable ylab = qualifier("ylabel",NULL);
   variable clevel = qualifier("clevel",0.9);
   variable transf = qualifier("transform",&unity);
   if(clevel >= 1) clevel=0.9;
   if(xlab==NULL)
   {
      xlabel("\\fr "+get_par_info(free_par[i]).name);
   }
   else
   {
      xlabel(xlab);
   }
   if(ylab==NULL)
   {
      ylabel("\\fr Probability");
   }
   else
   {
      ylabel(ylab);
   };
   variable lo,hi,lo_fine_wide,hi_fine_wide, real_xmin, real_xmax;
   variable nbin = qualifier("nbin",50);
   variable fcut = qualifier("fcut",0.2);
   variable lcut = qualifier("lcut",0);
   variable col = [qualifier("color",1)][0];
   variable lstyle = [qualifier("linestyle",1)][0];
   variable l = length(y[i,*]);
   variable ikeep=[int(fcut*l):l-1-int(lcut*l)];
   variable ntot = 1.*length(ikeep);
   variable ykeep;
   if( qualifier_exists("function") && qualifier_exists("func_vars") ) 
   {
      variable func = qualifier("function",&unity);
      variable func_vars = qualifier("func_vars",i);
      ykeep = @func(y[func_vars,ikeep]);
   }
   else
   {
      ykeep = @transf(y[i,ikeep]);
   }

   real_xmin = min(ykeep);
   real_xmax = max(ykeep);
   variable xmin = qualifier("pmin",real_xmin);
   variable xmax = qualifier("pmax",real_xmax);
   (lo,hi) = linear_grid(xmin,xmax,nbin);
   (lo_fine_wide,hi_fine_wide) = linear_grid(real_xmin,real_xmax,4*nbin);
   variable n = histogram(ykeep,lo,hi);
   variable n_fine_wide = histogram(ykeep,lo_fine_wide,hi_fine_wide);

   variable csum = cumsum(n_fine_wide)/ntot;

   % determine MAP or median values of posterior distribution
   if(qualifier_exists("return_map") || qualifier_exists("return_median")) {    
     variable grid_mean = (lo_fine_wide+hi_fine_wide)/2.;
     % P. Gregory, "Bayesian Logical Data Analysis for the Physical Sciences": median robust measure!
     variable p_map = grid_mean[where(n_fine_wide==max(n_fine_wide))][0];
     if (length(wherelast(csum <= 0.5))>=1 && wherelast(csum <= 0.5) != NULL) {
       variable p_median = grid_mean[wherelast(csum <= 0.5)][0];
     }
     else {
       p_median = p_map;
       message("Distribution unresolved with current binning. Adapt parameter ranges or change binning. Using MAP value.");
     }
     vmessage("MAP  %.2f   MEDIAN  %.2f",p_map,p_median);
   }
   ifnot(qualifier_exists("oplot") || qualifier_exists("return_map") || qualifier_exists("return_median"))
   { 
      linestyle(lstyle);
      hplot(lo,hi,n/(hi-lo)/ntot,col);
   }
   else ifnot(qualifier_exists("return_map") || qualifier_exists("return_median"))
   {
      linestyle(lstyle);
      ohplot(lo,hi,n/(hi-lo)/ntot,col);
   }
   % If posterior distribution is not clearly centrally peaked but shows secondary maxima,
   % the cummulative sum over the whole distribution may not be the best way.
   % 
   if(qualifier_exists("print"))
   {
      variable ifirst = wherelast(csum <= (1-clevel)/2.);
      variable ilast = wherefirst(csum >= (1+clevel)/2.);
      if(length(ifirst) > 0 && length(ilast) > 0)
      {
         () = printf("\n %.3e confidence levels are: %.4e %.4e \n\n",
                         clevel,hi_fine_wide[ifirst],lo_fine_wide[ilast]);
      }
      else
      {
         () = printf("\n Failed finding confidence levels.");
      }
   }
   if (qualifier_exists("return_map")) return p_map;
   else if (qualifier_exists("return_median")) return p_median;
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%

#ifnexists chain_set_params
public define chain_set_params()
{
  if(_NARGS==1)
  {
    variable y;
    (y) = ();
  }
  else
  {
    () = fprintf(fs,  "%s\n", %{{{
`
		 chain_set_params(chain; );

  Set parameters based on the posterior probability distributions gained
  in the emcee chain. This function loops over all free parameters and
  uses the function chain_hist to extract either the MAP or the median
  values by setting the appropriate qualifiers that are passed to the
  chain_hist function. 		 
`		 
		); %}}}
    return;
  }
  variable fp = free_par;
  variable i;
  _for i (0,length(fp)-1,1) {
    set_par(fp[i],chain_hist(i,y;;__qualifiers));
  }  
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%

#ifnexists vec_hist
public define vec_hist()
{
   if(_NARGS==1)
   {
      variable i, y;
      y = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  vec_hist(chain_vector; fcut=#, lcut=#, nbin=#, color=#, linestyle=#,
                       pmin=#, pmax=#, 
                       xrange={# or NULL,# or NULL}, 
                       yrange={# or NULL,# or NULL},
                       clevel=#, print, 
                       xlabel=String, ylabel=String, oplot); 

  Plot a normalized histogram for a vector from the MCMC chain,
  between values pmin and pmax (defaults to the min/max values in the
  chain; normalization, however, is always for the full range of the
  chain subsample). The plot axes have the ranges given in the lists,
  xrange and yrange (NULL==autoscale). Throw out the first fcut
  fraction of the data (default=0.2) and the last lcut fraction of the
  data (default=0), and use nbin bins (default=50).

  Plot using pgplot color given by qualifier (default=1, black) and
  linestyle given by qualifier (default=1, solid line).  The plot will
  be an overplot if the oplot qualifier is set.  The print qualifier
  has the lower and upper clevel (default=0.9) confidence levels
  printed out.
`); %}}}
      return;
   }
   variable xr = qualifier("xrange",{NULL,NULL});
   variable yr = qualifier("yrange",{NULL,NULL});
   xrange(xr[0],xr[1]); yrange(yr[0],yr[1]);
   connect_points(-1);
   variable fp = free_par;
   variable xlab = qualifier("xlabel","\\frX");
   variable ylab = qualifier("ylabel","\\frY");
   variable clevel = qualifier("clevel",0.9);
   if(clevel >= 1) clevel=0.9;

   xlabel(xlab);
   ylabel(ylab);

   variable lo,hi,lo_fine_wide,hi_fine_wide, real_xmin, real_xmax;
   variable nbin = qualifier("nbin",50);
   variable fcut = qualifier("fcut",0.2);
   variable lcut = qualifier("lcut",0);
   variable col = [qualifier("color",1)][0];
   variable lstyle = [qualifier("linestyle",1)][0];
   variable l = length(y);
   variable ikeep=[int(fcut*l):l-1-int(lcut*l)];
   variable ntot = 1.*length(ikeep);
   variable ykeep = y[ikeep];
   real_xmin = min(ykeep);
   real_xmax = max(ykeep);
   variable xmin = qualifier("pmin",real_xmin);
   variable xmax = qualifier("pmax",real_xmax);
   (lo,hi) = linear_grid(xmin,xmax,nbin);
   (lo_fine_wide,hi_fine_wide) = linear_grid(real_xmin,real_xmax,4*nbin);
   variable n = histogram(ykeep,lo,hi);
   variable n_fine_wide = histogram(ykeep,lo_fine_wide,hi_fine_wide);
   ifnot(qualifier_exists("oplot"))
   { 
      linestyle(lstyle);
      hplot(lo,hi,n/(hi-lo)/ntot,col);
   }
   else
   {
      linestyle(lstyle);
      ohplot(lo,hi,n/(hi-lo)/ntot,col);
   }
   if(qualifier_exists("print"))
   {
      variable csum = cumsum(n_fine_wide)/ntot;
      variable ifirst = wherefirst(csum >= (1-clevel)/2.);
      variable ilast = wherelast(csum >= (1+clevel)/2.);
      if(length(ifirst) > 0 && length(ilast) > 0)
      {
         () = printf("\n %.3e confidence levels are: %.4e %.4e \n\n",
                         clevel,hi_fine_wide[ifirst],lo_fine_wide[ilast]);
      }
      else
      {
         () = printf("\n Failed finding confidence levels.");
      }
   }
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifnexists chain_plot
public define chain_plot()
{
   if(_NARGS==4)
   {
      variable i, nc, nw, ans;
      (i, nc, nw, ans) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  chain_plot(i, nc, nw, chain; xrange={# or NULL,# or NULL}, 
                               yrange={# or NULL,# or NULL},
                               xlabel=String, ylabel=String, 
                               color=[#], transform=&func, oplot,
                               ); 

  Plot the nc^th parameter chain for the i^th free parameter from an
  MCMC chain ensemble that was made up of nw walkers per free
  parameter.  Alternatively, plot the chain with a function applied,
  passed as a reference via the transform qualifier (e.g.,
  transform=&log10). Plot using pgplot color given by qualifier
  (default=[1], black).
  The plot ranges can be set via the
  xrange/yrange qualifiers, and the plot labels can be overridden via
  the xlabel/ylabel qualifiers.  The plot will be an overplot if the
  oplot qualifier is set. The plot function presumes that the
  appropriate fit model has been loaded.
`); %}}}
      return;
   }
   connect_points(0);
   variable xr = qualifier("xrange",{NULL,NULL});
   variable yr = qualifier("yrange",{NULL,NULL});
   xrange(xr[0],xr[1]); yrange(yr[0],yr[1]);

   variable xlab = qualifier("xlabel",NULL);
   variable ylab = qualifier("ylabel",NULL);
   if(ylab==NULL)
   {
      ylabel("\\fr "+get_par_info(free_par[i]).name);
   }
   else
   {
      ylabel(ylab);
   }
   if(xlab==NULL)
   {
      xlabel("\\fr Iteration");
   }
   else
   {
      xlabel(xlab);
   };

   variable lans=length(ans[0,*]);
   variable lfp = length(free_par);
   variable col = qualifier("color",1);
   variable ichain=[nc-1:lans-1:lfp*nw];
   variable y=ans[i,ichain];
   variable x=[0:length(ichain)-1];
   variable transf = qualifier("transform",&unity);
   ifnot(qualifier_exists("oplot"))
   { 
     plot(x,@transf(y),col);
   }
   oplot(x,@transf(y),col);
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try { require( "xfig" ); }
catch AnyError: {};

#ifeval (__get_reference("xfig_plot_new")!=NULL)

#ifnexists xfig_chain_plot
public define xfig_chain_plot()
{
   if(_NARGS==1)
   {
      variable nu;
      (nu) = ();
   }  
   else if(_NARGS==4)
   {
      variable i, nc, nw, ans;
      (i, nc, nw, ans) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  xfig_chain_plot(i, nc, nw, chain; xrange={# or NULL,# or NULL}, 
                               yrange={# or NULL,# or NULL},
                               xlabel=String, ylabel=String, 
                               color=[#], transform=&func); 
  xfig_chain_plot(nu; xrange={# or NULL,# or NULL},
                               yrange={# or NULL,# or NULL},
                               xlabel=String, ylabel=String);

  Xfig version of chain_plot. See the appropriate help for further details.
  This function allows to plot multiple chains with different colors,
  if nc and col are given as arrays.
  The function will return a xfig plot struct.
  All qualifiers are passed to the xfig_plot.plot function.

  For only one argument (e.g. the acceptance fraction),
  this function plots the according array.
`); %}}}
      return;
   }

   variable xr,yr;
   variable col = qualifier("color",[1]);
   variable xlab = qualifier("xlabel",NULL);
   variable ylab = qualifier("ylabel",NULL);
   variable P = xfig_plot_new();
   if(_NARGS==4) {
     connect_points(0);
     if(ylab==NULL)
     {
       ylab = get_par_info(free_par[i]).name;
     }
     if(xlab==NULL)
     {
       xlab = "Iteration";
     }
     variable lans=length(ans[0,*]);
     variable lfp = length(free_par);
     xr = qualifier("xrange",{0,lans/lfp/nw});
     yr = qualifier("yrange",{0,1});
     P.world(xr[0],xr[1],yr[0],yr[1]);
     variable k;
     ifnot (typeof(nc)==Array_Type) nc = [nc];
     _for k (0,length(nc)-1,1) {
       variable ichain=[nc[k]-1:lans-1:lfp*nw];
       variable y=ans[i,ichain];
       variable x=[0:length(ichain)-1];
       y=y[where(y!=0)];
       x=x[where(y!=0)];          
       variable transf = qualifier("transform",&unity);
       P.plot(x,@transf(y) ;; struct_combine(__qualifiers, struct { color = (typeof(col)==Array_Type ? col[k] : col) }));
     }
%     return P;
   }
  
   if(_NARGS==1) {
     if(ylab==NULL)
     {
       ylab = "Acceptance fraction";
     }
     if(xlab==NULL)
     {
       xlab = "Iteration";
     }     
     xr = qualifier("xrange",{0,length(nu)});
     yr = qualifier("yrange",{0,1});
     P.world(xr[0],xr[1],yr[0],yr[1]);
     if (typeof(col)==Array_Type) col = col[0];
     variable idx = where(nu!=0);
     P.plot([1:length(nu[idx])],nu[idx] ;; struct_combine(__qualifiers, struct { color = col }));     
   }

   P.xlabel(xlab);
   P.ylabel(ylab);
   return P;

}
#endif

  

#endif


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If we are going to overplot a contour, we must use the same grid
% binning as last time!  (And the same transform function!)
private variable onbinx=50, onbiny=50, oxlo, oxhi, oylo, oyhi, otransfx, otransfy;

#ifnexists chain_hist2d
public define chain_hist2d()
{
   if(_NARGS==3)
   {
      variable i, j, y;
      (i, j, y) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  chain_hist2d(i, j, chain; fcut=#, lcut=#, nbinx=#, nbiny=#, 
                            transform_x=&func_x, transform_y=&func_y, 
                            function_x=&func_x, function_y=&func_y,  
                            func_vars_x = [Array], func_vars_y = [Array],
                            xlo=#, xhi=#, ylo=#, yhi=#,
                            color=Array, linestyle=Array,
                            levels=Array, xlabel=String, ylabel=String,
                            oplot, best_fit, return_cstruct); 

  Plot contours for a normalized, 2D histogram of the MCMC chain for
  the i^th (x) and j^th (y) free parameters. Throw out the first fcut
  fraction of the data (default=0.2) and the last lcut fraction of the
  data (default=0), and use nbinx/nbiny bins for each axis
  (default=50), which are histogrammed over the min/max values of the
  chain, or the input xlo/xhi/ylo/yho values.  

  The x- and y-axes can be transformed by passing functions as
  references via the transform_x and transform_y qualifiers,
  respectively (e.g., transform_x=&log10).

  Alternatively, multivariate functions can be passed via the
  function_x/y qualifiers, with func_vars_x/y being arrays indicating
  which of the i/j free parameters serve as variables.  These choices
  supercede the parameter indices input, i/j, and any transforms
  input.  Each function should expect input of the form x[i,k] and
  y(j,l), where the i/j index is for the free parameters of the fit,
  and the k/l index is a specific instance from the chain.

  Contours will use the pgplot color index, color (default=1=black),
  linestyle, line_style (default=1=solid), and will be placed at the
  probability values in the array, levels (default=[0.68,0.9,0.95]).
  Use custom X/Y labels with the x/ylabel qualifiers.  Overplot on an
  existing contour if the oplot qualifier is set. (This also forces
  nbinx/nbiny and xlo/xhi/ylo/yhi and transform_x/transform_y values
  to be the same as the previously plotted values.)  If the best_fit
  qualifier is input, then an x-mark is placed at the location
  indicated by the parameter file (presumed to be the best-fit), not
  at the peak of the marginalized probability distribution.

  If the return_cstruct qualifier is set, the structure containing the
  plotted contours will be returned with the additional field levels
  containing the minimum number of walkers (in one histogram bin)
  corresponding to the given confidence levels.
`); %}}}
      return;
   }
   connect_points(-1);
   variable nbinx = qualifier("nbinx",50);
   variable nbiny = qualifier("nbiny",50);
   variable fcut = qualifier("fcut",0.2);
   variable lcut = qualifier("lcut",0);
   variable levels = qualifier("levels",[0.68,0.9,0.95]);
   variable pg_col = [qualifier("color",1)];
   variable pg_lstyle = [qualifier("linestyle",1)];
   pg_lstyle = [qualifier("line_style",pg_lstyle)];
   variable xlab=qualifier("xlabel",NULL);
   variable ylab=qualifier("ylabel",NULL);
   variable transfx = qualifier("transform_x",&unity);
   variable transfy = qualifier("transform_y",&unity);
   variable xbest, ybest;

   if(qualifier_exists("oplot"))
   {
      nbinx = onbinx; nbiny = onbiny; transfx = otransfx; transfy = otransfy;
   }
   else
   {
      onbinx = nbinx; onbiny = nbiny; otransfx = transfx; otransfy = transfy;
   }

   if(xlab==NULL)
   {
      xlabel("\\fr "+get_par_info(free_par[i]).name);
   }
   else
   {
      xlabel(xlab);
   }

   if(ylab==NULL)
   {
      ylabel("\\fr "+get_par_info(free_par[j]).name);
   }
   else
   {
      ylabel(ylab);
   }

   variable xgrid, ygrid, xlo, xhi, ylo, yhi, ycut;
   ycut = length(y[i,*]);
   ycut = [[int(fcut*ycut):ycut-1-int(lcut*ycut)]];

   variable xkeep, ykeep, func, func_vars;
   if( qualifier_exists("function_x") && qualifier_exists("func_vars_x") ) 
   {
      func = qualifier("function_x",&unity);
      func_vars = qualifier("func_vars_x",i);
      xkeep = @func(y[func_vars,ycut]);
   }
   else
   {
      xkeep = @transfx(y[i,ycut]);
   }
   if( qualifier_exists("function_y") && qualifier_exists("func_vars_y") ) 
   {
      func = qualifier("function_y",&unity);
      func_vars = qualifier("func_vars_y",j);
      ykeep = @func(y[func_vars,ycut]);
   }
   else
   {
      ykeep = @transfy(y[j,ycut]);
   }

   if(qualifier_exists("oplot"))
  {
      xlo = oxlo;
      xhi = oxhi;
  }
  else
  {
      xlo = min(xkeep);
      xhi = max(xkeep);
      xlo = 1.*qualifier("xlo",xlo);
      xhi = 1.*qualifier("xhi",xhi);
      oxlo = xlo;
      oxhi = xhi;
  }
   xgrid = [0:nbinx]*(xhi-xlo)/(nbinx-1) + xlo; 

   if(qualifier_exists("oplot"))
   {
      ylo = oylo;
      yhi = oyhi;
   }
   else
   {
      ylo = min(ykeep);
      yhi = max(ykeep);
      ylo = 1.*qualifier("ylo",ylo);
      yhi = 1.*qualifier("yhi",yhi);
      oylo = ylo;
      oyhi = yhi;
   }
   ygrid = [0:nbiny]*(yhi-ylo)/(nbiny-1) + ylo; 

   variable n = histogram2d( ykeep, xkeep, ygrid, xgrid );
   variable icull=[0:nbinx-1];
   variable ii;
   _for ii (1,nbiny-1,1) { icull = [icull, ii*(nbinx+1)+[0:nbinx-1]]; }
   variable nn = n[icull];
   reshape (nn,[nbiny,nbinx]);

   variable n1d = @nn;
   reshape(n1d,[length(n1d)]);
   
   variable isort = array_sort(n1d);
   n1d = reverse(n1d[isort]);
   variable nprob = cumsum(n1d);
   nprob /= nprob[-1];

   isort = array_sort(levels);
   levels = levels[isort];

   variable iw, lev, clevels = Double_Type[0], plevels= Double_Type[0];

   foreach lev (levels)
   {
      iw = where(nprob >= lev);
      if(length(iw)>0) 
      {
         clevels = [clevels,n1d[iw[0]]];
         plevels = [plevels,lev];
      }
   }

%   print("Plotting contour levels:");
%   print(plevels);

   variable wfmn = wherefirstmax(nn);
   variable iy = wfmn/nbinx;
   variable ix = wfmn mod nbinx;

   point_style(5); 
   variable ogpo = get_plot_options;
   variable pconf = struct{width, type, color};
   variable cstruct = struct{chisqr, px, py, best, px_best, py_best};
   cstruct.best=0;
   cstruct.px = conf_grid(i,xgrid[1],xgrid[nbinx],nbinx);
   cstruct.py = conf_grid(j,ygrid[1],ygrid[nbiny],nbiny);
   cstruct.chisqr = nn;
   cstruct.px_best = xgrid[0];
   cstruct.py_best = ygrid[0];
   pconf.width = ogpo.line_width;
   pconf.type = pg_lstyle;
   pconf.color = pg_col;
   xlin; ylin; 
   if(qualifier_exists("oplot"))
   {
      % Weird ISIS thing, where the color doesn't seem to get set in
      % an oplot_contour, so do the point first to set the color.
      if(qualifier_exists("best_fit"))
      {
         xbest = @transfx(get_par(free_par[i]));
         ybest = @transfy(get_par(free_par[j]));
         connect_points(0); oplot(xbest,ybest,pg_col);
      }
      else
      {
         connect_points(0); oplot((xgrid[ix]+xgrid[ix+1])/2,(ygrid[iy]+ygrid[iy+1])/2,pg_col);
      }
      connect_points(1); oplot_conf(cstruct,pconf,clevels);
   }
   else
   {
      connect_points(1); plot_conf(cstruct,pconf,clevels);

      if(qualifier_exists("best_fit"))
      {
         xbest = @transfx(get_par(free_par[i]));
         ybest = @transfy(get_par(free_par[j]));
         connect_points(0); oplot(xbest,ybest,pg_col);
      }
      else
      {
         connect_points(0); oplot((xgrid[ix]+xgrid[ix+1])/2,(ygrid[iy]+ygrid[iy+1])/2,pg_col);
      }
   }
   set_plot_options(ogpo);

   if (qualifier_exists("return_cstruct")) {
     return struct_combine(cstruct, struct { levels = clevels });
   }
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If we are going to overplot a contour, we must use the same grid
% binning as last time!  (And the same transform function!)
private variable onbinx_vec=50, onbiny_vec=50, oxlo_vec, oxhi_vec, oylo_vec, oyhi_vec;

#ifnexists vec_hist2d
public define vec_hist2d()
{
   if(_NARGS==2)
   {
      variable yi, yj;
      (yi, yj) = ();
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  vec_hist2d(chain_vector1, chain_vector1; fcut=#, lcut=#, nbinx=#, nbiny=#, 
                            xlo=#, xhi=#, ylo=#, yhi=#, color=#,  line_style=#,
                            levels=Array, xlabel=String, ylabel=String, oplot,
                            return_cstruct); 

  Plot contours for a normalized, 2D histogram for two vectors from
  the MCMC chain.  Throw out the first fcut fraction of the data
  (default=0.2) and the last lcut fraction of the data (default=0),
  and use nbinx/nbiny bins for each axis (default=50), which are
  histogrammed over the min/max values of the chain, or the input
  xlo/xhi/ylo/yho values.

  Contours will use the pgplot color index, color (default=1=black),
  linestyle, line_style (default=1=solid), and will be placed at the
  probability values in the array, levels (default=[0.68,0.9,0.95]).
  Use custom X/Y labels with the x/ylabel qualifiers.  Overplot on an
  existing contour if the oplot qualifier is set. (This also forces
  nbinx/nbiny and xlo/xhi/ylo/yhi and transform_x/transform_y values
  to be the same as the previously plotted values.)

  If the return_cstruct qualifier is set, the function will return the 
  structure used to plot the confidence contours.
`); %}}}
      return;
   }
   connect_points(-1);
   variable nbinx = qualifier("nbinx",50);
   variable nbiny = qualifier("nbiny",50);
   variable fcut = qualifier("fcut",0.2);
   variable lcut = qualifier("lcut",0);
   variable levels = qualifier("levels",[0.68,0.9,0.95]);
   variable pg_col = qualifier("color",1);
   variable pg_lstyle = qualifier("linestyle",1);
   pg_lstyle = qualifier("line_style",pg_lstyle);
   variable xlab=qualifier("xlabel","\\frX");
   variable ylab=qualifier("ylabel","\\frY");
   variable xbest, ybest;

   if(qualifier_exists("oplot"))
   {
      nbinx = onbinx_vec; nbiny = onbiny_vec; 
   }
   else
   {
      onbinx_vec = nbinx; onbiny_vec = nbiny; 
   }

   xlabel(xlab);
   ylabel(ylab);

   variable xgrid, ygrid, xlo, xhi, ylo, yhi, ycut;
   ycut = length(yi);
   ycut = [[int(fcut*ycut):ycut-1-int(lcut*ycut)]];

   variable xkeep = yi[ycut];
   variable ykeep = yj[ycut];

   if(qualifier_exists("oplot"))
   {
      xlo = oxlo_vec;
      xhi = oxhi_vec;
   }
   else
   {
      xlo = min(xkeep);
      xhi = max(xkeep);
      xlo = 1.*qualifier("xlo",xlo);
      xhi = 1.*qualifier("xhi",xhi);
      oxlo_vec = xlo;
      oxhi_vec = xhi;
   }
   xgrid = [0:nbinx]*(xhi-xlo)/(nbinx-1) + xlo; 

   if(qualifier_exists("oplot"))
   {
      ylo = oylo_vec;
      yhi = oyhi_vec;
   }
   else
   {
      ylo = min(ykeep);
      yhi = max(ykeep);
      ylo = 1.*qualifier("ylo",ylo);
      yhi = 1.*qualifier("yhi",yhi);
      oylo_vec = ylo;
      oyhi_vec = yhi;
   }
   ygrid = [0:nbiny]*(yhi-ylo)/(nbiny-1) + ylo; 

   variable n = histogram2d( ykeep, xkeep, ygrid, xgrid );
   variable icull=[0:nbinx-1];
   variable ii;
   _for ii (1,nbiny-1,1) { icull = [icull, ii*(nbinx+1)+[0:nbinx-1]]; }
   variable nn = n[icull];
   reshape (nn,[nbiny,nbinx]);

   variable n1d = @nn;
   reshape(n1d,[length(n1d)]);
   
   variable isort = array_sort(n1d);
   n1d = reverse(n1d[isort]);
   variable nprob = cumsum(n1d);
   nprob /= nprob[-1];

   isort = array_sort(levels);
   levels = levels[isort];

   variable iw, lev, clevels = Double_Type[0], plevels= Double_Type[0];

   foreach lev (levels)
   {
      iw = where(nprob >= lev);
      if(length(iw)>0) 
      {
         clevels = [clevels,n1d[iw[0]]];
         plevels = [plevels,lev];
      }
   }

%   print("Plotting contour levels:");
%   print(plevels);

   variable wfmn = wherefirstmax(nn);
   variable iy = wfmn/nbinx;
   variable ix = wfmn mod nbinx;

   point_style(5); 
   variable ogpo = get_plot_options;
   variable pconf = struct{width, type, color};
   variable cstruct = struct{chisqr, px, py, best, px_best, py_best};
   cstruct.best=0;
   cstruct.px = conf_grid(1,xgrid[1],xgrid[nbinx],nbinx);
   cstruct.py = conf_grid(2,ygrid[1],ygrid[nbiny],nbiny);
   cstruct.chisqr = nn;
   cstruct.px_best = xgrid[0];
   cstruct.py_best = ygrid[0];
   pconf.width = ogpo.line_width;
   pconf.type = pg_lstyle;
   pconf.color = pg_col;

   xlin; ylin; 
   if(qualifier_exists("oplot"))
   {
      % Weird ISIS thing, where the color doesn't seem to get set in
      % an oplot_contour, so do the point first to set the color.

      connect_points(0); oplot((xgrid[ix]+xgrid[ix+1])/2,(ygrid[iy]+ygrid[iy+1])/2,pg_col);
      connect_points(1); oplot_conf(cstruct,pconf,clevels);  %oplot_contour(nn,clevels);
   }
   else
   {
      connect_points(1); plot_conf(cstruct,pconf,clevels); 
      % plot_contour(nn,1,xgrid[[1:nbinx]],ygrid[[1:nbiny]],clevels);
      connect_points(0); oplot((xgrid[ix]+xgrid[ix+1])/2,(ygrid[iy]+ygrid[iy+1])/2,pg_col);
   }
   set_plot_options(ogpo);

   if(qualifier_exists("return_cstruct")) return cstruct;
}
#endif

#ifnexists plot_covar
public define plot_covar()
{ 
   if(_NARGS==2)
   {
      variable args = __pop_list(_NARGS);
   }
   else
   {
      () = fprintf(fs,  "%s\n", %{{{
`
  id=plot_covar([pars], chain; fcut=#, lcut=#, nbin=Array, 
                        prange={{lo_0,hi_0},{lo_1,hi_1},...}, 
                        hrange={{lo_0,hi_0},{lo_1,hi_1},...}, 
                        plabel=String_Array, hlabel=String_Array, 
                        pmin=Array, pmax=Array, 
                        levels=Array, color=Array, linestyle=Array,
                        device=String, window=#, 
                        viewport=#, size=#, aspect=#,
                        charsize=#, point_size=#, best_fit);

  Create an NxN grid where along the diagonal 1D probability
  histograms are plotted, and the lower left corner contains the
  N(N-1)/2 2D probability contours for the results of an MCMC chain.
  The return value id is the identifier for the opened pgplot device.

  Parameters to be plotted are passed in the [pars] array, and are
  numbered by their order in the chain (i.e., [0,1,2, ...]).  

  Throw out the first fcut fraction of the data (default=0.2) and the last
  lcut fraction of the data (default=0).  Use nbin[i] bins for
  each parameter axis (default=50). 

  pmin[i], pmax[i] retain the same meaning as in chain_hist(), except
  are now arrays, with each value corresponding to the input
  parameters.

  prange[i] is a List of Lists, with each entry being the range over
  which each parameter will have its histogram/contour plotted.
  (Defaults are pmin[i]/pmax[i], or the min/max value of the parameter
  in the portion of the chain that is used.)

  hrange[i] is a List of Lists, with each entry being the range over
  which each parameter *historgram* will have its y-axis plotted
  (Defaults are to autoscale.)

  plabel[i] is a String_Array of the axis labels associated with each
  parameter.  (Default is to take the parameter name from the model
  parameter file, which is presumed to have been loaded.)

  hlabel[i] is a String_Array of the y-axis labels associated with each
  parameter for the histogram plots (default="\\fr Probability").

  Contours will be plotted for the probability levels given by the
  levels array (default=[0.68,0.9,0.95]).

  Histograms/contours will use the pgplot color indices in the color
  array, with the first entry being for the probability histograms,
  and the remaining entries are for the contour levels
  (default=1=black for all). The linestyle array serves the same
  purpose for the lines (default=1=solid).

  The device qualifier indicates the pgplot device (default="/xw").
  Use, for example, device="covar_plot.ps/vcps" to create a postscript
  plot.  (A postscript plot will later have to be closed via
  close_plot(id).)  The window=id qualifier allows an existing pgplot
  device to be reused, but *not* have its basic grid structure (i.e.,
  number of parameters plotted) changed.

  viewport is the standard pgplot definition *for each column*
  (default=[0.1,0.95,0.1,0.98]).  size and aspect are the inputs to
  the resize(size,aspec) command *for the whole* window (defaults=20,
  1, respectively).

  charsize and point_size are the usual pgplot commands (default=1).
`); %}}}
      return;
   }
   variable pars = args[0];
   variable chain = args[1];
   variable n = length(pars);
   variable pr={}, pr_tmp = qualifier("prange",{{NULL,NULL}});
   variable hr={}, hr_tmp = qualifier("hrange",{{NULL,NULL}});
   variable plab=String_Type[n], 
            plab_tmp = qualifier("plabel",String_Type[0]);
   variable hlab=String_Type[n], 
            hlab_tmp = [qualifier("hlabel","\\fr Probability")];
   variable pbin=Integer_Type[n], 
            pbin_tmp=[qualifier("nbin",Integer_Type[0])];
   variable fcut = qualifier("fcut",0.2);
   variable lcut = qualifier("lcut",0);
   variable l = length(chain[0,*]);
   variable ikeep = [int(fcut*l):l-1-int(lcut*l)];
   variable pmin=Double_Type[n], pmin_tmp = [qualifier("pmin",Double_Type[0])];
   variable pmax=Double_Type[n], pmax_tmp = [qualifier("pmax",Double_Type[0])];
   variable levels = qualifier("levels",[0.68,0.9,0.95]);
   variable col = [qualifier("color",1)];
            if(length(col)<length(levels)+1)
            { 
               col=[col[[0:length(col)-1]],
                    Integer_Type[length(levels)-length(col)+1]+col[-1]];
            }
   variable pg_lstyle = [qualifier("linestyle",1)];
            pg_lstyle = [qualifier("line_style",pg_lstyle)];
            if(length(pg_lstyle)<length(levels)+1)
            {
               pg_lstyle=[pg_lstyle[[0:length(pg_lstyle)-1]],
                          Integer_Type[length(levels)-length(pg_lstyle)+1]
                          +pg_lstyle[-1]];
            }
   variable csize = qualifier("charsize",1);
   variable psize = qualifier("point_size",1);
            psize = qualifier("pointsize",psize);
   variable device=qualifier("device","/xw");
   variable size=qualifier("size",20);
   variable aspect=qualifier("aspect",1);
   variable id;

   % Simple error checks
   if( typeof(pr_tmp)!=List_Type && 
      (typeof(pr_tmp[0])!=List_Type && pr_tmp[0]!=NULL) )
   { message("\n prange must be a List of Lists \n"); return; }
   if( typeof(hr_tmp)!=List_Type && 
      (typeof(hr_tmp[0])!=List_Type && hr_tmp[0]!=NULL) )
   { message("\n hrange must be a List of Lists \n"); return; }
   if( typeof(plab_tmp)!=Array_Type && 
      (typeof(plab_tmp[0])!=String_Type && plab_tmp[0]!=NULL) )
   { message("\n plabel must be Array of Strings \n"); return; }
   if( typeof(hlab_tmp)!=Array_Type && 
      (typeof(hlab_tmp[0])!=String_Type && hlab_tmp[0]!=NULL) )
   { message("\n hlabel must be Array of Strings \n"); return; }
   if( typeof(pbin_tmp)!=Array_Type && 
      (typeof(pbin_tmp[0])!=Integer_Type && pbin_tmp[0]!=NULL) )
   { message("\n nbin must be an Array of Integers\n"); return; }
   if( typeof(pmin_tmp)!=Array_Type && 
      (typeof(pmin_tmp[0])!=Double_Type && pmin_tmp[0]!=NULL) )
   { message("\n pmin must be an Array of Doubles\n"); return; }
   if( typeof(pmax_tmp)!=Array_Type && 
      (typeof(pmax_tmp[0])!=Double_Type && pmax_tmp[0]!=NULL) )
   { message("\n pmax must be an Array of Doubles\n"); return; }

   variable i, j, iarray;
   
   % Set the defaults
   _for i (0,n-1,1)
   {
      try{
         if(pmin_tmp[i]!=NULL){ pmin[i]=pmin_tmp[i]; }
         else{ pmin[i] = min(chain[i,ikeep]);}
      }
      catch AnyError: {
         pmin[i] = min(chain[i,ikeep]);
      }
      try{
         if(pmax_tmp[i]!=NULL){ pmax[i]=pmax_tmp[i]; }
         else{ pmax[i] = max(chain[i,ikeep]);};
      }
      catch AnyError: {
         pmax[i] = max(chain[i,ikeep]);
      }
      try{ list_append(pr,pr_tmp[i]); }
      catch AnyError: { list_append(pr,{NULL,NULL}); }
      if(length(pr[i])!=2) { pr[i] = {pmin[i],pmax[i]}; }
      if(pr[i][0]==NULL) { pr[i][0]=pmin[i]; }
      if(pr[i][1]==NULL) { pr[i][1]=pmax[i]; }
      try{ list_append(hr,hr_tmp[i]); }
      catch AnyError: { list_append(hr,{NULL,NULL}); }
      try{ 
         if(plab_tmp[i]==NULL)
         {  plab[i]="\\fr "+get_par_info(free_par[pars[i]]).name; }
         else { plab[i]=plab_tmp[i]; }
      }
      catch AnyError: { 
         plab[i]="\\fr "+get_par_info(free_par[pars[i]]).name; 
      }
      try{
         hlab[i]=hlab_tmp[i];
      }
      catch AnyError: {
         hlab[i]="\\fr Probability";
      }
      try{ 
         if(pbin_tmp[i]==NULL)
         {  pbin[i]=50; }
         else { pbin[i]=50; }
      }
      catch AnyError: { 
         pbin[i]=50;
      }
   }

   if(qualifier_exists("window"))
   {
      id = qualifier("window");
      window(id);
   }
   else
   {
      id=open_plot(device,n,1);
   }
   variable ov=qualifier("viewport",[0.1,0.95,0.1,0.98]);
   resize(size,aspect);
   charsize(csize);
   point_size(psize);
%
   variable ov_strct = struct{xmin,xmax,ymin,ymax};
   ov_strct.xmin=ov[0];
   ov_strct.xmax=ov[1];
   ov_strct.ymin=ov[2];

   _for i (0,n-1,1)
   {
      iarray=Integer_Type[n-i]+1;
      multiplot(iarray);
      ov_strct.ymax=ov[2]+(n-i)/(1.*n)*(ov[3]-ov[2]);

      set_outer_viewport(ov_strct);

      _for j (0,n-i-1,1)
      {
         mpane(j+1);
         if(j==0)
         {
            chain_hist(pars[i],chain;fcut=fcut,lcut=lcut,nbin=pbin[i],
               pmin=pmin[i],pmax=pmax[i],xrange=pr[i],yrange=hr[i],
               color=col[0],linestyle=pg_lstyle[0],ylabel=hlab[i]);
         }
         else
         {
            if(qualifier_exists("best_fit"))
            {
               chain_hist2d(pars[i],pars[j+i],chain;fcut=fcut,lcut=lcut,
                  nbinx=pbin[i],nbiny=pbin[j+i],xlo=pr[i][0],
                  xhi=pr[i][1],ylo=pr[j+i][0],yhi=pr[j+i][1],
                  color=col[[1:length(col)-1]],
                  linestyle=pg_lstyle[[1:length(pg_lstyle)-1]],
                  levels=levels,xlabel=plab[i],ylabel=plab[j+i]);
            } 
            else
            {
               chain_hist2d(pars[i],pars[j+i],chain;fcut=fcut,lcut=lcut,
                  nbinx=pbin[i],nbiny=pbin[j+i],xlo=pr[i][0],
                  xhi=pr[i][1],ylo=pr[j+i][0],yhi=pr[j+i][1],
                  color=col[[1:length(col)-1]],
                  linestyle=pg_lstyle[[1:length(pg_lstyle)-1]],
                  levels=levels,xlabel=plab[i],ylabel=plab[j+i]);
            } 
         }
      }
   }
   return id;
}
#endif
