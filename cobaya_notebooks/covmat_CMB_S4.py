## Functions for calculating covariance matrix for CMB-S4 experiment 
from collections import defaultdict
import numpy as np
import fishchips
from scipy import linalg
import fishchips.experiments
import fishchips.cmb_lensing
import math

def δ(α, β):
    return np.array(α==β, dtype=float)

def cov_2pt(signal_dict, noise_dict, bin_left, bin_right=None, 
            f_sky=0.6, N_splits=2,
            instruments=('Planck', 'Planck', 'Planck', 'Planck'), 
            observables=('T', 'T', 'T', 'T'), 
            seasons=(0,0,0,0),
            verbose=True):
    """Compute a covariance matrix for a 2-pt correlation function.
    
    This is an implementation of equation B2 in arXiv:1610.02360.
    
    signal_dict : dict
        Must have tuple keys for all two-point (T/E/B/κ) pairs
        i.e. a key could be ('T', 'E').
    noise_dict : dict
        Must have tuple keys in the format
            (season, instrument, season, instrument, T/E/B/κ, T/E/B/κ)
        If you have only one season, set every season value to 0.
    """
    
    # if bin_right is not specified, bin_right[n] = bin_left[n+1]
    if bin_right == None: 
        Δbin = bin_left[-1] - bin_left[-2] # last bin-width
        bin_right = np.hstack((np.array(bin_left[1:]), [bin_left[-1]+Δbin]))
        
    A, B, C, D = instruments
    W, X, Y, Z = observables
    α, β, γ, τ = seasons
    
    N_s = N_splits
    n_b = len(bin_left)
    
    # assume zero signal and infinite noise for non-specified terms
    S_b = defaultdict(lambda:0.0, signal_dict)
    N_b = defaultdict(lambda:np.inf, noise_dict)
    
    ν_b = np.zeros(n_b)
    for i, (bl, br) in enumerate(zip(bin_left, bin_right)):
        ν_b[i] = np.sum(2 * np.arange(bl, br) + 1) * f_sky
    
    # deliberately not PEP8, I think this is more readable
    sum1 = S_b[W,Y] * S_b[X,Z] + S_b[W,Z] * S_b[X,Y]
    
    sum2 = (
        S_b[W,Y] * δ(β,τ) * N_b[β,B,τ,D,X,Z] + 
        S_b[X,Z] * δ(α,γ) * N_b[α,A,γ,C,W,Y] +
        S_b[W,Z] * δ(β,γ) * N_b[β,B,γ,C,X,Y] + 
        S_b[X,Y] * δ(α,τ) * N_b[α,A,τ,D,W,Z]
    )
    
    sum3 = (
        δ(α,γ) * δ(β,τ) * N_b[α,A,γ,C,W,Y] * N_b[β,B,τ,D,X,Z] +
        δ(β,γ) * δ(α,τ) * N_b[α,A,τ,D,W,Z] * N_b[β,B,γ,C,X,Y]
    )
    
    prefactor_3 = (1/ν_b) * (
        (N_s**2 - N_s * (δ(α,β) + δ(γ,τ)) + N_s * δ(α,β) * δ(γ,τ)) /
        (N_s**4 - N_s**3 * (δ(α,β) + δ(γ,τ)) + N_s**2 * δ(α,β) * δ(γ,τ))
    )
    
    return (1/ν_b) * sum1 + (1/(N_s * ν_b)) * sum2 + prefactor_3 * sum3

def fiducial_model(omega_b = 0.0224, omega_cdm = 1e-22, omega_dmeff = 0.1197, m_dmeff = 1,npow_dmeff = 0, sigma_dmeff = 0,
                   theta_s = 1.0416e-2, tau_reio = 0.07, A_s = 2.2e-9, n_s = 0.96, N_ur = 2.0328, 
                   N_ncdm = 1, m_ncdm = 0.06, T_ncdm = 0.71611):
    fiducial_params = {
        'omega_b': omega_b, 'omega_cdm': omega_cdm, 'omega_dmeff': omega_dmeff, 'm_dmeff': m_dmeff,
        'npow_dmeff' : npow_dmeff, 'sigma_dmeff': sigma_dmeff,
        '100*theta_s': 1e2*theta_s, 'tau_reio': tau_reio,
        'A_s': A_s, 'n_s': n_s,
        'N_ur': N_ur, 'N_ncdm': N_ncdm, 'm_ncdm': m_ncdm, 'T_ncdm': T_ncdm,}
    l_min = 2
    l_max = 2500

    modules_path = '/home/zequnl/src/cobaya_modules'

    info_fiducial = {
        'params': fiducial_params,
        'likelihood': {'one': None},
        'theory': {
            'classy' : 
                   {'extra_args': 
                    {  'output': 'tCl,pCl,lCl,mPk',
                        'l_max_scalars': 5000,
                        'lensing': 'yes',
                        #'omega_cdm':1e-22,
                        #'m_dmeff':1,
                        #'npow_dmeff' : 0,


                        'tight_coupling_trigger_tau_c_over_tau_k':0.,
                        'tight_coupling_trigger_tau_c_over_tau_h':0.,
                        'reionization_optical_depth_tol': 1e-07,
#                         'tol_background_integration': 1e-8,
#                         'tol_perturb_integration': 1e-8,
#                         'tol_thermo_integration': 1e-9,
                        'perturb_sampling_stepsize':0.01,
                        'k_max_tau0_over_l_max' : 6,
                        'P_k_max_h/Mpc' : 5.,
                        'gauge' : 'synchronous',
                        'k_per_decade_for_pk' : 100} } 

        }
    }
    from cobaya.model import get_model
    model_fiducial = get_model(info_fiducial)

    # Declare our desired theory product
    # (there is no cosmological likelihood doing it for us)
    
    # model_fiducial.likelihood.theory.needs(Cl={'tt': l_max, 'ee': l_max})

    # Compute and extract the CMB power spectrum
    # (In muK^-2, without l(l+1)/(2pi) factor)
    # notice the empty dictionary below: all parameters are fixed
    
    # model_fiducial.logposterior({})
    # Cls = model_fiducial.likelihood.theory.get_Cl(ell_factor=False)

    # Our fiducial power spectrum
    
    # Cl_est_tt = Cls['tt'][:l_max+1]
    # Cl_est_ee = Cls['ee'][:l_max+1]
    # Cl_est_pp = Cls['pp'][:l_max+1]
    # Ell = Cls['ell'][:l_max+1]
    # Cl_est_kk = []
    # for i in Ell:
        # if i < 2:
            # Cl_est_kk.append(0)
        # else:
            # Cl_est_kk.append(1/4*(math.factorial(i+2)/math.factorial(i-2))*Cl_est_pp[i])
    return model_fiducial #, np.array([Cls, Cl_est_tt, Cl_est_ee, Cl_est_pp, Cl_est_kk])
    
def noise_Planck_Pol(l_max = 2500):
    
    exper = fishchips.experiments.get_PlanckPol_combine()
    Nltt = np.concatenate((exper[0].noise_T, exper[1].noise_T[31:], exper[2].noise_T[101:]))[:l_max+1]
    Nlee = np.concatenate((exper[1].noise_P, exper[2].noise_P[101:]))[:l_max+1]
    Nlkk = fishchips.cmb_lensing.CMB_Lensing_Only().noise_k[:l_max+1]
    
    return np.array([Nltt, Nlee, Nlkk])

def noise_CMB_S4(l_max = 2500):
    exper = fishchips.experiments.get_S4()
    Nltt = exper[1].noise_T[:l_max+1]
    Nlee = np.concatenate((exper[0].noise_P, exper[1].noise_P[301:3001], exper[2].noise_P[3001:]))[:l_max+1]
    Nlkk = fishchips.experiments.get_S4_Lensing_Only().noise_k[:l_max+1]
    
    return np.array([Nltt, Nlee, Nlkk])


def calculate_covmat(experiment, l_min = 2, l_max = 2500):
    model_fiducial = fiducial_model();
    Cltt = model_fiducial.likelihood.theory.get_Cl(ell_factor=False)['tt'][:l_max+1]
    Clee = model_fiducial.likelihood.theory.get_Cl(ell_factor=False)['ee'][:l_max+1]
    Clte = model_fiducial.likelihood.theory.get_Cl(ell_factor=False)['te'][:l_max+1]
    Clpp = model_fiducial.likelihood.theory.get_Cl(ell_factor=False)['pp'][:l_max+1]
    Ell = model_fiducial.likelihood.theory.get_Cl(ell_factor=False)['ell'][:l_max+1]
    Clkk = []
    for i in Ell:
        if i <= 2:
            Clkk.append(0)
        else:
            Clkk.append(1/4*(math.factorial(i+2)/math.factorial(i-2))*Clpp[i-2])
    Clkk = np.array(Clkk) 
    Cltk = np.zeros(len(Ell))
    Clke = np.zeros(len(Ell))
    
    if experiment == "Planck":
        noise_array = noise_Planck_Pol()
    elif experiment == "CMB-S4":
        noise_array = noise_CMB_S4()
        
    Nltt = noise_array[0]
    Nlee = noise_array[1]
    Nlkk = noise_array[2]
    Nlte = np.zeros(l_max+1)[:l_max+1]
    Nltk = np.zeros(l_max+1)[:l_max+1]
    Nlke = np.zeros(l_max+1)[:l_max+1]
    
    # calculate the components with temperature (T) and polarization (E)
    covTTTT = cov_2pt(
        {('T', 'T'): Cltt},
        {(0, 'S4', 0, 'S4', 'T', 'T'): Nltt},
        bin_left=np.arange(l_max+1),
        observables=('T', 'T', 'T', 'T'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covTTTE = cov_2pt(
        {('T', 'T'): Cltt, ('E', 'E'): Clee,
         ('T', 'E'): Clte, ('E', 'T'): Clte},
        {(0, 'S4', 0, 'S4', 'T', 'T'): 
         Nltt, 
        (0, 'S4', 0, 'S4', 'E', 'E'): 
         Nlee,
        (0, 'S4', 0, 'S4', 'T', 'E'): 
         Nlte,
        (0, 'S4', 0, 'S4', 'E', 'T'): 
         Nlte
        },
        bin_left=np.arange(l_max+1),
        observables=('T', 'T', 'T', 'E'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covTETT = covTTTE

    covTTEE = cov_2pt(
        {('T', 'T'):
          Cltt, 
         ('E', 'E'): 
          Clee,
         ('T', 'E'): 
          Clte,
         ('E', 'T'): 
          Clte
        },
        {(0, 'S4', 0, 'S4', 'T', 'T'): 
         Nltt, 
        (0, 'S4', 0, 'S4', 'E', 'E'): 
         Nlee,
        (0, 'S4', 0, 'S4', 'T', 'E'): 
         Nlte,
        (0, 'S4', 0, 'S4', 'E', 'T'): 
         Nlte
        },
        bin_left=np.arange(l_max+1),
        observables=('T', 'T', 'E', 'E'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covEETT = covTTEE


    covTETE = cov_2pt(
        {('T', 'T'):
          Cltt, 
         ('E', 'E'): 
          Clee,
         ('T', 'E'): 
          Clte,
         ('E', 'T'): 
          Clte
        },
        {(0, 'S4', 0, 'S4', 'T', 'T'): 
         Nltt, 
        (0, 'S4', 0, 'S4', 'E', 'E'): 
         Nlee,
        (0, 'S4', 0, 'S4', 'T', 'E'): 
         Nlte,
        (0, 'S4', 0, 'S4', 'E', 'T'): 
         Nlte
        },
        bin_left=np.arange(l_max+1),
        observables=('T', 'E', 'T', 'E'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )

    covTEEE = cov_2pt(
        {('T', 'T'):
          Cltt, 
         ('E', 'E'): 
          Clee,
         ('T', 'E'): 
          Clte,
         ('E', 'T'): 
          Clte
        },
        {(0, 'S4', 0, 'S4', 'T', 'T'): 
         Nltt, 
        (0, 'S4', 0, 'S4', 'E', 'E'): 
         Nlee,
        (0, 'S4', 0, 'S4', 'T', 'E'): 
         Nlte,
        (0, 'S4', 0, 'S4', 'E', 'T'): 
         Nlte
        },
        bin_left=np.arange(l_max+1),
        observables=('T', 'E', 'E', 'E'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covEETE = covTEEE

    covEEEE = cov_2pt(
        {('E', 'E'): Clee},
        {(0, 'S4', 0, 'S4', 'E', 'E'): Nlee},
        bin_left=np.arange(l_max+1),
        observables=('E', 'E', 'E', 'E'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    
    # Compute components with lensing (K)
    covKKKK = cov_2pt(
        {('K', 'K'): Clkk},
        {(0, 'S4', 0, 'S4', 'K', 'K'): Nlkk},
        bin_left=np.arange(l_max+1),
        observables=('K', 'K', 'K', 'K'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )

    covTTKK = cov_2pt(
        {('T', 'T'):
          Cltt, 
         ('K', 'K'): 
          Clkk,
         ('T', 'K'): 
          Cltk,
         ('K', 'T'): 
          Cltk
        },
        {(0, 'S4', 0, 'S4', 'T', 'T'): 
         Nltt, 
        (0, 'S4', 0, 'S4', 'K', 'K'): 
         Nlkk,
        (0, 'S4', 0, 'S4', 'T', 'K'): 
         Nltk,
        (0, 'S4', 0, 'S4', 'K', 'T'): 
         Nltk
        },
        bin_left=np.arange(l_max+1),
        observables=('T', 'T', 'K', 'K'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covKKTT = covTTKK

    covTEKK = cov_2pt(
        {('T', 'T'): Cltt, ('K', 'K'): Clkk, ('E', 'E'): Clee,
         ('T', 'E'): Clte, ('E', 'T'): Clte,
         ('T', 'K'): Cltk, ('K', 'T'): Cltk, 
         ('E', 'K'): Clke, ('K', 'E'): Clke, 
        },
        {(0, 'S4', 0, 'S4', 'T', 'T'): 
         Nltt, 
        (0, 'S4', 0, 'S4', 'E', 'E'): 
         Nlee,
        (0, 'S4', 0, 'S4', 'K', 'K'): 
         Nlkk,
        (0, 'S4', 0, 'S4', 'T', 'E'): 
         Nlte,
        (0, 'S4', 0, 'S4', 'E', 'T'): 
         Nlte,
        (0, 'S4', 0, 'S4', 'E', 'K'): 
         Nlke,
        (0, 'S4', 0, 'S4', 'K', 'E'): 
         Nlke,
        (0, 'S4', 0, 'S4', 'K', 'T'):
         Nltk,
        (0, 'S4', 0, 'S4', 'T', 'K'): 
         Nltk
        },
        bin_left=np.arange(l_max+1),
        observables=('T', 'E', 'K', 'K'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covKKTE = covTEKK

    covEEKK = cov_2pt(
        {('E', 'E'):
          Clee, 
         ('K', 'K'): 
          Clkk,
         ('E', 'K'): 
          Clke,
         ('K', 'E'): 
          Clke
        },
        {(0, 'S4', 0, 'S4', 'K', 'K'): 
         Nlkk, 
        (0, 'S4', 0, 'S4', 'E', 'E'): 
         Nlee,
        (0, 'S4', 0, 'S4', 'K', 'E'): 
         Nlke,
        (0, 'S4', 0, 'S4', 'E', 'K'): 
         Nlke
        },
        bin_left=np.arange(l_max+1),
        observables=('E', 'E', 'K', 'K'), 
        instruments=('S4', 'S4', 'S4', 'S4'), 
        seasons=(0,0,0,0) # season is always 0 here
    )
    covKKEE = covEEKK

    cov = np.block([[np.diag(covTTTT[l_min:]), np.diag(covTTTE[l_min:]), np.diag(covTTEE[l_min:]), np.diag(covTTKK[l_min:])],
                   [np.diag(covTETT[l_min:]), np.diag(covTETE[l_min:]), np.diag(covTEEE[l_min:]), np.diag(covTEKK[l_min:])], 
                   [np.diag(covEETT[l_min:]), np.diag(covEETE[l_min:]), np.diag(covEEEE[l_min:]), np.diag(covEEKK[l_min:])], 
                   [np.diag(covKKTT[l_min:]), np.diag(covKKTE[l_min:]), np.diag(covKKEE[l_min:]), np.diag(covKKKK[l_min:])]])

    invcov = linalg.cho_solve( linalg.cho_factor(cov), b = np.identity(cov.shape[0]) )
    sign_of_logdet, logdet = np.linalg.slogdet(cov)
    return invcov, sign_of_logdet, logdet



