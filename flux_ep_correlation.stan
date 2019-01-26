functions {


  vector background_model(vector observed_counts, vector background_counts, vector background_error, vector expected_model_counts) {

    int N = num_elements(expected_model_counts);
    
    vector[N] MB = background_counts + expected_model_counts;
    vector[N] s2 = square(background_error);

    vector[N] b = 0.5 * (sqrt(square(MB) - 2 * s2 .* (MB - 2 * observed_counts) + square(s2))
	       + background_counts - expected_model_counts - s2);
    return b;
    

  }

  
  real pgstat(vector observed_counts, vector background_counts, vector background_error, vector expected_model_counts, int[] idx_background_zero, int[] idx_background_nonzero) {

    int N = num_elements(expected_model_counts);
    vector[N] log_likes;


    
    vector[N] s2 = square(background_error);

    vector[N] b = background_model(observed_counts, background_counts, background_error, expected_model_counts);

    vector[N] factorial_term = expected_model_counts - lgamma(observed_counts +1 );
    
    //print(expected_model_counts);

 
    
    log_likes[idx_background_nonzero] = (-square(b[idx_background_nonzero] - background_counts[idx_background_nonzero]) ./ (2 * s2[idx_background_nonzero])
			+ observed_counts[idx_background_nonzero] .* log(b[idx_background_nonzero] + expected_model_counts[idx_background_nonzero])
			- b[idx_background_nonzero] - factorial_term[idx_background_nonzero]
			- 0.5 * log(2 * pi()) - log(background_error[idx_background_nonzero]));	


    for (n in 1:num_elements(idx_background_zero)) {

    
    log_likes[idx_background_zero[n]] = lmultiply(observed_counts[idx_background_zero[n]], expected_model_counts[idx_background_zero[n]]) - factorial_term[idx_background_zero[n]];	
    }
    
    
    return sum(log_likes);

  }
  
  

  real ggrb_int_pl(real alpha, real beta, real ec, real emin, real emax) {

    real pre = pow(alpha - beta, alpha - beta) * exp(beta - alpha) / pow(ec, beta);

    if (beta !=-2) {
    
      return pre/(2.+beta) * (pow(emax, 2+ beta) - pow(emin, 2+ beta));
    }

    else {

      return pre * log(emax/emin);
    }
  }

  
  real ggrb_int_cpl(real alpha, real ec, real emin, real emax) {

    real i1 = gamma_q(2 + alpha, emin / ec) * tgamma(2 + alpha);
    real i2 = gamma_q(2 + alpha, emax / ec) * tgamma(2 + alpha);

    return -square(ec) * (i2-i1);

    
  }

  
  vector differential_flux(vector energy, real flux, real alpha, real beta, real epeak, real emin, real emax) {
    
    real ec;

    real esplit;
    
    real erg2keV = 6.24151e8;
    real intflux;
    real norm;
    vector[num_elements(energy)] out;


    if (alpha !=-2.) {

      ec = epeak / (2. + alpha);
      
    }

    else {

      ec = epeak/.0001;
      
    }
    

    esplit = (alpha - beta) * ec;
    
    if ((emin<= esplit)  &&  (esplit <=emax)) {
      
      intflux = (ggrb_int_cpl(alpha, ec, emin, esplit) + ggrb_int_pl(alpha, beta, ec, esplit, emax));
      
    }
    
    else if (esplit < emin) {
      
      intflux = ggrb_int_pl(alpha, beta, ec, esplit, emax);
      
    }

    else {


      print(emin);
      print(esplit);
      print(emax);

    }

	
    norm = flux * erg2keV / intflux;


    for (n in 1:num_elements(energy)) {
      
      if (energy[n] < esplit) {
	out[n] = pow(energy[n] / ec, alpha) * exp(-energy[n] / ec);
      }
      
      else {
	out[n] = pow(alpha - beta, alpha - beta) * exp(beta - alpha) * pow(energy[n] / ec, beta);
      }

      
    }

    return norm * out;
    
  }
  
  vector integral_flux(vector ebounds_lo, vector ebounds_hi, real flux, real alpha, real  beta, real epeak, real emin, real emax) {

    return (ebounds_hi - ebounds_lo) / 6.0
      .* (differential_flux(ebounds_lo, flux, alpha, beta, epeak, emin, emax)
	  + 4 * differential_flux(0.5*(ebounds_lo + ebounds_hi), flux,  alpha,  beta, epeak, emin, emax)
	  + differential_flux(ebounds_hi, flux, alpha, beta, epeak, emin, emax));
      
  }
  
  
}

data {

  int<lower=1> N_intervals;
  int max_n_echan;
  int max_n_chan;

  int<lower=0> N_dets[N_intervals]; // number of detectors poer data type
  int<lower=0> N_chan[N_intervals, max(N_dets)]; // number of channels in each detector
  int<lower=0> N_echan[N_intervals,  max(N_dets)]; // number of energy side channels in each detector

  int grb_id[N_intervals];
  int N_grbs;
  
  vector[max_n_echan] ebounds_hi[N_intervals, max(N_dets)];
  vector[max_n_echan] ebounds_lo[N_intervals, max(N_dets)];
  


  vector[max_n_chan] observed_counts[N_intervals, max(N_dets)];
  vector[max_n_chan] background_counts[N_intervals, max(N_dets)];
  vector[max_n_chan] background_errors[N_intervals, max(N_dets)];

  int idx_background_zero[N_intervals, max(N_dets), max_n_chan];
  int idx_background_nonzero[N_intervals, max(N_dets), max_n_chan];
  int N_bkg_zero[N_intervals,max(N_dets)];
  int N_bkg_nonzero[N_intervals,max(N_dets)];
  
  real exposure[N_intervals, max(N_dets)];

  matrix[max_n_echan, max_n_chan] response[N_intervals, max(N_dets)];


  
  int mask[N_intervals, max(N_dets), max_n_chan];
  int N_channels_used[N_intervals,max(N_dets)];

  vector[N_intervals] dl;
  vector[N_intervals] z;


  int N_gen_spectra;
  vector[N_gen_spectra] model_energy;

  int N_correlation;
  vector[N_correlation] model_correlation;
  
  
}

transformed data {
  vector[N_intervals] dl2 = square(dl); 

  real emin = 10.;
  real emax = 1E6;
    


  /* vector[N_intervals] emin; */
  /* vector[N_intervals] emax; */

  /* emin = 10. ./ (1+z); */
  /* emax = 1.E5 ./ (1+z); */
  

  
}



parameters {

  vector<lower=-1.8, upper=1.>[N_intervals] alpha;
  vector<lower=-6., upper=-2.>[N_intervals] beta;
  vector<lower=1, upper=1E4>[N_intervals] epeak;
  //vector<lower=0>[N_intervals] energy_flux;

  vector[N_grbs] gamma_offset;
  real gamma_mu;
  real<lower=0> gamma_sigma;


  vector[N_grbs] delta_offset;
  real delta_mu;
  real<lower=0> delta_sigma;

  
  

  
}

transformed parameters {

  
  vector[N_grbs] gamma;
  vector[N_grbs] delta;
  vector[N_intervals] log_energy_flux;
  vector[max_n_chan] expected_model_counts[N_intervals, max(N_dets)];
 
  gamma = gamma_mu + gamma_offset * gamma_sigma;
  delta = delta_mu + 52  + delta_offset * delta_sigma;

  
  log_energy_flux = delta[grb_id] + gamma[grb_id] .* log10(epeak .* (1+z)/100.) - 2*log10(dl[grb_id]) - log10(4*pi());   
  
  
  // compute the folded counts
  
    for (n in 1:N_intervals) {

      for (m in 1:N_dets[n]) {	
	
	expected_model_counts[n,m,:N_chan[n,m]] = ((to_row_vector(integral_flux(ebounds_lo[n, m, :N_echan[n, m]],
										ebounds_hi[n, m, :N_echan[n, m]],
										10^log_energy_flux[n],
										alpha[n],
										beta[n],
										epeak[n],
										emin,
										emax)) * response[n, m,:N_echan[n,m],:N_chan[n,m]]) * exposure[n,m])';
	
      }
    }
  
  
}


model {


  
  alpha ~ normal(-1,.5);
  beta ~ normal(-3,1);
  epeak ~ normal(500.,500.);




  gamma_mu ~ normal(1.5, 1);
  delta_mu ~ normal(0,4);
  
  gamma_offset ~ std_normal();
  delta_offset ~ std_normal();

  gamma_sigma ~ normal(0, 5);
  delta_sigma ~ normal(0, 5);

  /* gamma ~ normal(0,10); */
  /* delta_raw ~ normal(0,3); */
  
  for (n in 1:N_intervals) {

    for (m in 1:N_dets[n]) {
      
      
      target += pgstat(observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
		       background_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
		       background_errors[n, m, mask[n,m,:N_channels_used[n,m]]],
		       expected_model_counts[n,m, mask[n,m,:N_channels_used[n,m]]],
		       idx_background_zero[n,m, :N_bkg_zero[n,m]],
		       idx_background_nonzero[n,m, :N_bkg_nonzero[n,m]]  );
      
    }
    
  }
  
  
  
}


generated quantities {

  vector[N_gen_spectra] vfv_spectra[N_intervals];
  vector[max_n_chan] count_ppc[N_intervals, max(N_dets)];

  vector[N_correlation] correlations[N_grbs];

  for (n in 1:N_intervals) {

    vfv_spectra[n] =square(model_energy) .* differential_flux(model_energy, 10^log_energy_flux[n], alpha[n], beta[n], epeak[n], emin, emax);


      
 
    

    for (m in 1:N_dets[n]) {

        vector[N_channels_used[n,m]] ppc_background = background_model(observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
								       background_counts[n, m, mask[n,m,:N_channels_used[n,m]]],
								       background_errors[n, m, mask[n,m,:N_channels_used[n,m]]],
								       expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]]);
	
	vector[N_channels_used[n,m]] rate = ppc_background + expected_model_counts[n, m, mask[n,m,:N_channels_used[n,m]]] ;
	for (i in 1:N_channels_used[n,m]) {

	  if (rate[i]>2^30) {
	    

	    count_ppc[n,m,i] = 0;
	    
	  }
	  
	  else {

	    
	    count_ppc[n,m,i] = poisson_rng( rate[i] );
	    
	  }
	  
	}
	
    }
    
    
    
  
  }    
  
  
  for (n in 1:N_grbs) {

    correlations[n] = delta[n] + gamma[n] * model_correlation;
    
  }
  
  
}

  
