functions {

  real pgstat(vector observed_counts, vector background_counts, vector background_error, vector expected_model_counts, int[] idx_background_zero, int[] idx_background_nonzero) {

    int N = num_elements(expected_model_counts);
    vector[N] log_likes;


    
    
    
    vector[N] MB = background_counts + expected_model_counts;
    vector[N] s2 = square(background_error);

    vector[N] b = 0.5 * (sqrt(square(MB) - 2 * s2 .* (MB - 2 * observed_counts) + square(s2))
	       + background_counts - expected_model_counts - s2);


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

    if (intflux < 0) {
      print({alpha, beta, epeak, flux, emin, emax, ec, esplit});

    }
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
  
  
  vector convole_drm(matrix expected_flux, matrix response ) {
    
    return rows_dot_product(expected_flux,response); 
  }
}

data {

  
  int<lower=1> N_shards;
  
  int<lower=1> N_intervals;
  int max_n_echan;
  int max_n_chan;

  int<lower=0> N_dets[N_intervals]; // number of detectors poer data type
  int<lower=0> N_chan[N_intervals, max(N_dets)]; // number of channels in each detector
  int<lower=0> N_echan[N_intervals,  max(N_dets)]; // number of energy side channels in each detector

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

  
  
}

transformed data {
  vector[N_intervals] dl2 = square(dl); 
  vector[N_intervals] emin;
  vector[N_intervals] emax;

  emin = 10. ./ (1+z);
  emax = 1.E7 ./ (1+z);
  

  
}



parameters {

  vector<lower=-2., upper=1.>[N_intervals] alpha;
  vector<upper=-2.>[N_intervals] beta;
  vector<lower=1, upper=1E5>[N_intervals] epeak;
  //vector<lower=0>[N_intervals] energy_flux;

  real gamma;
  real delta_raw;


  
}

transformed parameters {

  vector[N_intervals] energy_flux;
  vector[N_intervals] log_luminosity;
  real log_delta;

  log_delta = delta_raw + 52;
  
  log_luminosity = log_delta + gamma*log10( (1+z).*epeak );

  for (n in 1:N_intervals) {
    energy_flux[n] = 10^log_luminosity[n] / (4*pi()*dl2[n]);
  }


  // compute the folded counts
  
  

  
}


model {


  
  alpha ~ normal(-1,1);
  beta ~ normal(-3,1);
  epeak ~ normal(500.,500.);
  //  energy_flux ~ normal(1E-6,1E-2);

  gamma ~ normal(0,10);
  delta_raw ~ normal(0,3);
  
  for (n in 1:N_intervals) {

    for (m in 1:N_dets[n]) {
      
      vector[N_echan[n,m]] tmp_flux = integral_flux(ebounds_lo[n, m, :N_echan[n, m]], ebounds_hi[n, m, :N_echan[n, m]], energy_flux[n], alpha[n], beta[n], epeak[n], emin[n], emax[n]);


 
      
      row_vector[N_chan[n,m]] expected_model_counts = (to_row_vector(tmp_flux) * response[n, m,:N_echan[n,m],:N_chan[n,m]]) * exposure[n,m];

      vector[N_channels_used[n,m]] selected_counts = observed_counts[n, m, mask[n,m,:N_channels_used[n,m]]]; 
      vector[N_channels_used[n,m]] selected_background = background_counts[n, m, mask[n,m, :N_channels_used[n,m]] ]; 
      vector[N_channels_used[n,m]] selected_background_errors = background_errors[n, m, mask[n,m,:N_channels_used[n,m]]];
      row_vector[N_channels_used[n,m]] selected_expectation = expected_model_counts[mask[n,m,:N_channels_used[n,m]]];


      
      target += pgstat(selected_counts, selected_background, selected_background_errors, selected_expectation', idx_background_zero[n,m, :N_bkg_zero[n,m]], idx_background_nonzero[n,m, :N_bkg_nonzero[n,m]]  );
      
    }

  }

  
  
}


generated quantities {
  
}


