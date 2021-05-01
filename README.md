# Fisher analysis for constraining dark matter+baryon interaction in early universe. 

The purpose for this repo is to constrain cosmolocial parameters the CMB simulated in the Boltzman code CLASS. In particular I have been constraining the 6 LCDM parameters in addition the the early-universe dark matter-baryon scattering cross section using a modified version of CLASS. With the intention to forcast sensitivities of future CMB experiment by plotting upper limits for crassections in mass ranging between 10keV - 1TeV (See: upperlimit_crossection.ipynb)

<b>compute_fisher_dmb</b> computes the fisher matrix for a given experiment with required input of dark matter (DM) mass, velocity dependence, and relagtive bulk velocity which uses simulated cosmological observables from CLASS.

Each of these notebooks (for given experiment) contains details on the equations used and their references. Each one has input paramters that are easily modified for velocity dependance (n_power/n_dmb), relative bulk velocity (relative_bulk_velocity/Vrel_dmb), dark matter mass (dm_mass/m_dmb), sky coverage (f_sky) etc.

The fiducal values for the LCDM paramters are also easily changed, as well as which one to include/not include.

Each setting above requres adjustment of the step size used in numerically differenting the C_ls in the fisher matrix (see details in notebooks).

To numerically differentiation the Cl's with respect to the parameters you want to contrain its important to choose the right stepsize. One that is small enough to converge but also large enough not to hit numerical noise. See 'convergence_check.ipynb'

## Formalism
 

[Wu et al. 2014](https://arxiv.org/abs/1402.4108)


## Experimental noise
Each experiment has their own sets of noise curves:
##Sources: 

<b>Simons observatory</b>
Lensing noise: https://github.com/simonsobs/so_noise_models/tree/master/LAT_comp_sep_noise/v3.1.0
Temperature and Polarization: https://github.com/simonsobs/so_noise_models

<b>CMB STAGE-4</b>
https://github.com/xzackli/fishchips-public

<b>Cosmic variance limited</b>

<b>Planck</b>
https://github.com/xzackli/fishchips-public
