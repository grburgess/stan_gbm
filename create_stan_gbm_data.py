from astropy.cosmology import WMAP9 as cosmo
from threeML import *
import pystan
from glob import glob
import pandas as pd


def get_tr_datalists(grb,pha):
    
   
    
    print(grb)
    
    dets = np.unique([f.split('/')[-1].split('_')[1]   for f in  pha])
    
    with fits.open(pha[0]) as f:
        
        n_bins = int(f[1].header['NAXIS2'])
    
    
     
    
    datalists = []

    set = False
    for i in range(n_bins):


        
        pi = []
        
        for det in dets:
        
            obs = '%s/prepared_pha_files/%s_%s_time-resolved.pha'%(grb,grb,det)
            bak = '%s/prepared_pha_files/%s_%s_time-resolved_bak.pha'%(grb,grb,det)
            rsp = '%s/prepared_pha_files/%s_%s_time-resolved.rsp'%(grb,grb,det)
        
            oo = OGIPLike(det,observation=obs,background=bak,response=rsp,verbose=False, spectrum_number=i+1)
            
            if det[0] =='b':
                
                oo.set_active_measurements('250-30000')
                oo.use_effective_area_correction(.2, 1.8)

            elif det=='LLE':

                oo.set_active_measurements('40000-300000')
                oo.use_effective_area_correction(.5, 1.5)

                print('using lle!!!!!!!!!!!!')
                

            else:
                if not set:

                    set = True

                else:
                    oo.use_effective_area_correction(.5, 1.5)
                
                oo.set_active_measurements('8.1-30','35-900')
            
            pi.append(oo)




        #data.fix_effective_area_correction(1.)
            
        datalists.append(pi)
        
    return datalists



all_files = glob(os.path.join('GRB*','prepared*'  ))
single_pulse = filter(lambda x: '/prepared_pha_files_' not in x, all_files)
pha_files = [glob(os.path.join(p,'*.pha') ) for p in single_pulse]


all_data_sets = []
grb_name = []
grb_number = []

for n, pha in enumerate(pha_files):

    grb = pha[0].split('/')[0]

    
    
    this_data_set=get_tr_datalists(grb,pha=pha)

    grb_name.extend([grb]*len(this_data_set))
    grb_number.extend([n+1]*len(this_data_set))
    
    all_data_sets.extend(this_data_set)



data_dict = {}

max_n_dets = max([ len(ds) for ds in all_data_sets])
max_n_chans = max(np.max( [ [ len(d.observed_counts)  for d in  ds] for ds in all_data_sets]))
max_n_echans = max(np.max( [ [ d.response.matrix.shape[1]  for d in  ds] for ds in all_data_sets]))

n_intervals = len(all_data_sets)

data_dict['N_intervals'] = n_intervals


observed_counts = np.zeros((n_intervals,max_n_dets,max_n_chans))
background_counts = np.zeros((n_intervals,max_n_dets,max_n_chans))
background_errors = np.zeros((n_intervals,max_n_dets,max_n_chans))

idx_background_zero = np.zeros((n_intervals,max_n_dets,max_n_chans))
idx_background_nonzero = np.zeros((n_intervals,max_n_dets,max_n_chans))
n_bkg_zero = np.zeros((n_intervals,max_n_dets))
n_bkg_nonzero = np.zeros((n_intervals,max_n_dets))


responses = np.zeros((n_intervals,max_n_dets,max_n_echans,max_n_chans))

exposures = np.zeros((n_intervals,max_n_dets))
n_echan = np.zeros((n_intervals,max_n_dets))
n_chan = np.zeros((n_intervals,max_n_dets))
n_dets = [ len(ds) for ds in all_data_sets]



masks = np.zeros((n_intervals,max_n_dets,max_n_chans))
n_channels_used = np.zeros((n_intervals,max_n_dets))
grb_id = np.zeros(n_intervals)
ebounds_lo = np.zeros((n_intervals, max_n_dets, max_n_echans))
ebounds_hi = np.zeros((n_intervals, max_n_dets, max_n_echans))

cbounds_lo = np.zeros((n_intervals, max_n_dets, max_n_chans))
cbounds_hi = np.zeros((n_intervals, max_n_dets, max_n_chans))


dl = np.zeros(n_intervals)
z = np.zeros(n_intervals)

redshifts = pd.read_json('redshifts.json')



total_number_of_channels_used = 0
grb_counter = 0

print(len(grb_number))

for i,ds in enumerate(all_data_sets):

    z_idx = redshifts['trigger'] == grb_name[i] 

    
    z[i] = float(redshifts['z'][z_idx])
    dl[i] = cosmo.luminosity_distance(z[i]).to('cm').value
    
    
    for j, aa in enumerate(ds):
        
        observed_counts[i,j,:len(aa.observed_counts)] = aa.observed_counts
        background_counts[i,j,:len(aa.background_counts)] = aa.background_counts
        
        idx_background_zero[i,j,:sum(aa.background_counts==0)] = np.where(aa.background_counts==0)[0] +1
        idx_background_nonzero[i,j,:sum(aa.background_counts>0)] = np.where(aa.background_counts>0)[0] +1
        n_bkg_zero[i,j] = sum(aa.background_counts==0)
        n_bkg_nonzero[i,j] = sum(aa.background_counts>0)
        
        background_errors[i,j,:len(aa.background_count_errors)] = aa.background_count_errors
        
        responses[i, j, :aa.response.matrix.shape[1] ,:len(aa.background_count_errors)] = aa.response.matrix.T
        
        
        
        this_mask = np.where(aa.mask)[0] +1
        
        masks[i,j,:len(this_mask)] = this_mask
        
        idx_background_zero[i,j,:sum(aa.background_counts[aa.mask]==0)] = np.where(aa.background_counts[aa.mask]==0)[0] +1
        idx_background_nonzero[i,j,:sum(aa.background_counts[aa.mask]>0)] = np.where(aa.background_counts[aa.mask]>0)[0] +1
        n_bkg_zero[i,j] = sum(aa.background_counts[aa.mask]==0)
        n_bkg_nonzero[i,j] = sum(aa.background_counts[aa.mask]>0)
        
        
        n_channels_used[i,j] = len(this_mask)
        n_chan[i,j] = len(aa.observed_counts)
        n_echan[i,j] = aa.response.matrix.shape[1]
        ebounds_lo[i,j,:aa.response.matrix.shape[1]] = aa.response.monte_carlo_energies[:-1]
        ebounds_hi[i,j,:aa.response.matrix.shape[1]] = aa.response.monte_carlo_energies[1:]
        cbounds_lo[i,j,:aa.response.matrix.shape[0]] = aa.response.ebounds[:-1]
        cbounds_hi[i,j,:aa.response.matrix.shape[0]] = aa.response.ebounds[1:]
        exposures[i,j] = aa.exposure
        
        total_number_of_channels_used += sum(aa.mask)

    
    
        
data_dict['N_all'] = int(total_number_of_channels_used)
data_dict['N_grbs'] = len(np.unique(grb_number))
data_dict['max_n_echan'] = max_n_echans
data_dict['max_n_chan'] = max_n_chans
data_dict['N_dets'] = n_dets
data_dict['N_chan'] = n_chan.astype(int)
data_dict['N_echan'] = n_echan.astype(int)
data_dict['observed_counts'] = observed_counts
data_dict['background_counts'] = background_counts
data_dict['idx_background_nonzero'] = idx_background_nonzero.astype(int)
data_dict['idx_background_zero'] = idx_background_zero.astype(int)
data_dict['N_bkg_zero'] = n_bkg_zero.astype(int)
data_dict['N_bkg_nonzero'] = n_bkg_nonzero.astype(int)
data_dict['object_idx'] = np.array(grb_number).astype(int)
data_dict['background_errors'] = background_errors
data_dict['ebounds_lo'] = ebounds_lo
data_dict['ebounds_hi'] = ebounds_hi
data_dict['cbounds_lo'] = cbounds_lo
data_dict['cbounds_hi'] = cbounds_hi
data_dict['exposure'] = exposures
data_dict['response'] = responses
data_dict['mask'] = masks.astype(int)
data_dict['N_channels_used'] = n_channels_used.astype(int)
data_dict['grb_id'] = np.array(grb_number).astype(int)
data_dict['dl'] = dl
data_dict['z'] = z


pystan.stan_rdump(data_dict,'all_data.R')
