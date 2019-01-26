import stan_utility
import pystan
import numpy as np
model = stan_utility.compile_model('gbm_sum.stan', model_name='test_gbm2')

data = pystan.read_rdump('alpha_data.R')

N_gen_spectra = 100
model_energy = np.logspace(0,5,N_gen_spectra)
data['N_gen_spectra'] = N_gen_spectra
data['model_energy'] = model_energy

warmup = 1000
iter = 100

total = warmup + iter

chains = 8

fit = model.sampling(
    data=data,
    iter=total,
    warmup=warmup,
    chains=chains,
    n_jobs=chains,
    control=dict(max_treedepth=13,
    #             adapt_delta=0.9
    ),
    seed=1234)

stan_utility.stanfit_to_hdf5(fit, '/data/jburgess/stan_fits/gbm_stan_fit_small.h5')
