# Fisher analysis for constraining dark matter+baryon interaction in early universe. 

<b>compute_fisher_dmb</b> computes the fisher matrix for a given experiment with required input of dark matter (DM) mass, velocity dependence, and relagtive bulk velocity which uses simulated cosmological observables from CLASS.

## Formulism
 
$$ F_{ij} = \sum_{\ell} \frac{2 \ell + 1}{2} f_{\mathrm{sky}} \mathrm{Tr}\,\left( \mathbf{C}_{\ell}^{-1} \frac{\partial \mathbf{C}_{\ell}}{\partial \theta_i} \mathbf{C}_{\ell}^{-1} \frac{\mathbf{C}_{\ell}}{\partial \theta_j} \right)$$
[Wu et al. 2014](https://arxiv.org/abs/1402.4108)


## Experimental noise
Each experiment has their own sets of noise curves:
##Sources: 

<b>Simons observatory</b>

<b>CMB STAGE-4</b>

<b>Cosmic variance limited</b>

<b>Planck</b>
