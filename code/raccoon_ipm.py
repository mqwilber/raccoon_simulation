from __future__ import division
import numpy as np
import scipy.stats as stats
import macroeco.models as md

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

def growth_fxn(next_load, current_load, age, eggs, params):
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
    prob_worms_lost = stats.binom.pmf(worms_lost, 
                                      n=current_load, 
                                      p=1-np.exp(-params['worm_death_rate']))

    prob_worms_gained = _worms_gained(diff_worms + worms_lost, params,
                                      eggs, age, current_load)

    prob_next_load = np.sum(prob_worms_lost*prob_worms_gained)
    return(prob_next_load)

def _worms_gained(worms, params, eggs, age, current_load, 
                            max_encounter=10000):
    """
    Pick up worms from heterogeneously distributed latrines
    """

    if eggs != 0:

        # This reduces the mean eggs encountered and acquired by adjusting the mean.
        prob_inf = _worm_infection_prob(age, current_load, params)
        prob_worms_gained = md.nbinom.pmf(worms, 
                                mu=params['egg_contact']*eggs*prob_inf, 
                                k_agg=params['k_latrine'])

        # Exact summation...very slow and we don't need it.
        # prob_inf = _worm_infection_prob(age, current_load, params)
        # eggs_cont_vects = [np.arange(worm, max_encounter + 1) for worm in worms]
        # worm_probs = [md.nbinom.pmf(ec, mu=params['egg_contact']*eggs, k_agg=params['k_latrine'])*stats.binom.pmf(worm, n=ec, p=prob_inf) 
        #               for worm, ec in zip(worms, eggs_cont_vects)]

        # prob_worms_gained = np.nan_to_num([np.sum(wp) for wp in worm_probs])

    else:
        # If there are no eggs, only can encounter 0 eggs with prob = 1
        prob_worms_gained = (worms == 0).astype(np.int)

    return(prob_worms_gained)

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
        # TODO FIX THIS NO INFECTION LEADS TO NAN
        prob_inf = prob_inf * np.exp(-params['age_immunity'] * age) - params['age_resistance']))

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

    surv_prob = np.exp(death_thresh + patho*np.log(load + 1)) / (1 + np.exp(death_thresh + patho*np.log(load + 1)))
    return(surv_prob)



    
