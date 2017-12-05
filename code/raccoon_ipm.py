from __future__ import division
import numpy as np
import scipy.stats as stats

def survival_fxn(age, current_load, params):
    """
    Probability of survival given age and parasite load.

    Parameters
    ----------
    age : int
        Age in months of a raccoon
    current_load : int
        Number of parasites infecting a raccoon
    params : dict
        Dictionary of necessary parameters

    Returns
    -------
    : survival probability
    """

    surv_prob = _age_dependent_survival(age, params['baby_death'], params['old_death']) * \
                _load_dependent_survival(current_load, params['death_thresh'], params['patho']) * \
                (1 - params['random_death'])

    return(surv_prob)

def _growth_fxn(next_load, current_load, age, eggs, params):
    """
    Probability of worm load transitioning from current_load to next_load

    Given raccoon age and eggs in the environment.

    Parameters
    ----------
    next_load : int
        Number of parasites at the next time step
    current_load : int
        Raccoon parasites at the current time step
    age : int
        Age, in months, of a raccoon
    eggs : float
        Number of eggs in the environment at a given time
    params : dict
        Dictionary containing parameters for the growth function

    Returns
    -------
    : float
        Probability of of next_load given parameters

    """

    diff_worms = next_load - current_load

    worms_lost = np.arange(current_load + 1)
    prob_worms_lost = _worms_lost(worms_lost, current_load, params)

    prob_worms_gained = _worms_gained(diff_worms + worms_lost,
                                      eggs, age, current_load, params)

    prob_next_load = np.sum(prob_worms_lost*prob_worms_gained)
    return(prob_next_load)

growth_fxn = np.vectorize(_growth_fxn) # Vectorized version

def _worms_lost(worms_lost, current_load, params):
    """ 
    Worms lost due to worm death 

    Draw probabilities from binomial distribution
    """
    prob_worms_lost = stats.binom.pmf(worms_lost, 
                                      n=current_load, 
                                      p=1-np.exp(-params['worm_death_rate']))
    return(prob_worms_lost)


def _worms_gained(worms, eggs, age, current_load, params):
    """
    Pick up worms from heterogeneously distributed latrines OR rodents

    Parameters
    ----------
    worms : array
        Possible values of worms that are gained
    eggs : float
        Current eggs in the environment
    age : int
        Current age of the raccoon in months
    current_load : int
        Current worm load of the raccoon
    params : dict
        Dictionary holding the model parameters
    """

    worms = np.atleast_1d(worms)
    inds = worms >= 0
    prob_inf = _worm_infection_prob(age, current_load, params)

    # Gain worms from environment
    if age < params['age_resistance']:

        # This reduces the mean eggs encountered and acquired by adjusting the mean.
        zeros = np.sum(~inds)
        prob_worms_gained = np.r_[np.array([0]*zeros), nbinom_pmf(worms[inds], 
                                mu=params['egg_contact']*eggs*prob_inf, 
                                k=params['k_latrine'])]

    else: # Above a certain age, raccoons get worms from rodents or raccoons

        worm_vects = [np.arange(x + 1) for x in worms[inds]]

        # Call the NBD once for rodents and environment...speeds things up a lot
        environ_probs = nbinom_pmf(worm_vects[-1], 
                            mu=params['egg_contact']*eggs*prob_inf, 
                            k=params['k_latrine'])

        rodent_probs = nbinom_pmf(worm_vects[-1], 
                            mu=params['mean_rodent']*params['larval_infect']*params['rodent_contact'], 
                            k=params['k_rodent'])

        prob_worms_gained = [0 for x in worms[~inds]] + \
                            [np.sum(environ_probs[tworms]*rodent_probs[tworms[::-1]])
                                    for tworms in worm_vects]

    return(np.array(prob_worms_gained))

def _worm_infection_prob(age, load, params):
    """
    Probability of a contacted worm egg establishing as a worm in a raccoon
    """

    # Exponential decline in infectivity with load
    prob_inf = params['infect'] * np.exp(-params['resist'] * load)

    if age <= params['age_resistance']:
        prob_inf = prob_inf
    else:
        # Exponential decline in susceptibility with age after some age
        prob_inf = prob_inf * np.exp(-params['age_immunity'] * \
                                        (age - params['age_resistance']))

    return(prob_inf)

def _age_dependent_survival(age, baby_death, old_death):
    """
    Age-dependent survival probability
    """

    if age == 0:
        surv_prob = 1 - baby_death
    else:
        surv_prob = (1 - np.min([1, old_death*age**2]))

    return(surv_prob)

age_dependent_survival = np.vectorize(_age_dependent_survival)


def _load_dependent_survival(load, death_thresh, patho):
    """
    Parasite load-dependent survival fxn
    """

    surv_prob = np.exp(death_thresh + patho*np.log(load + 1)) / \
                        (1 + np.exp(death_thresh + patho*np.log(load + 1)))
    return(surv_prob)


def nbinom_pmf(x, mu, k):
    """ Replace nan with 1...annoying that scipy doesn't do this"""

    vals = np.atleast_1d(stats.nbinom._pmf(x, n=k, p=(k / (mu + k))))
    vals[np.isnan(vals)] = 1 
    return(vals)
    
