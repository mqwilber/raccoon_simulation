from __future__ import division
import numpy as np
import scipy.stats as stats
import multiprocessing as mp

def get_growth_matricies(params, max_age=100, max_load=100, processes=4):
    """
    Calculate age-dependent growth matrices with parallel processing

    Parameters
    ----------
    params : dict
        Dictionary with model parameters
    max_age : int
        Number of age-dependent growth matrices to calculate
    max_load : int
        Maximum parasite load to consider in matrices
    processes : int
        Number of cores to use

    Returns
    -------
    : list
        List of length max_age of parasite growth matrices
    """

    pool = mp.Pool(processes=processes)
    nxt, cur = np.meshgrid(range(max_load + 1), range(max_load + 1))
    Gmat = []

    for age in range(max_age):
    
        # This is what we want to parallelize over
        results = [pool.apply_async(_growth_fxn_mp, args=(samp_id, tnxt, tcur, age, 10000, params)) for 
                    samp_id, (tnxt, tcur) in enumerate(zip(nxt.ravel(), cur.ravel()))]

        # Format results
        results = [p.get() for p in results]
        results.sort()
        gvals = np.array(zip(*results)[1])

        Gmat.append(np.reshape(gvals, (len(nxt), len(nxt))).T)

    return(Gmat)

def get_survival_matrix(params, max_age=100, max_load=100):
    """
    Leslie matrix for age-dependent survival on the weekly time scale

    Parameters
    ----------
    params : dict
        Dictionary of parameters
    max_age : int
        Maximum number of ages in Leslie Matrix

    Returns
    -------
    : np.array of shape (max_load + 1,  max_age)
        Matrix where each column is the survival probability for a given age
        class after accounting for load, age, and random death.
    """

    # Load-dependent survival
    load_surv = _load_dependent_survival(np.arange(max_load + 1), params)

    # Age-dependent survival
    age_surv = age_dependent_survival(np.arange(max_age), params)

    # Age by load matrix survival matrix
    S = np.outer(load_surv, age_surv)
    S = S * (1 - params['random_death'])

    return(S)

def get_repro_matrix(params, pop_size, max_age=240):
    """
    Get the monthly reproduction matrix

    Parameters
    ----------
    params : dict
        Dictionary of parameters
    pop_size : float 
        Raccoon population size before reproduction
    max_age : int
        Maximum number of ages in Leslie Matrix

    Returns
    -------
    : np.array
        Leslie matrix with reproduction values on the first row and 1's on the
        lower off diagonal.
    

    """
    repro_vals = [reproduction_fxn(age, pop_size, params) for age in np.arange(max_age)]
    R = np.diag(np.repeat(1, max_age - 1), k=-1).astype(np.float)
    R[0, :] = repro_vals
    return(R)

def reproduction_fxn(age, pop_size, params):
    """
    Age-dependent expected number offspring 

    """
    age = np.atleast_1d(age)

    if age >= params['repro_age']:
        return(np.exp(-params['ricker_effect']*pop_size) * params['litter_size'])
    else:
        return(0)


def survival_fxn(age, current_load, params):
    """
    Probability of survival given age and parasite load. Weekly survival.

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

    surv_prob = _age_dependent_survival(age, params) * \
                _load_dependent_survival(current_load, params) * \
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

def _growth_fxn_mp(samp_id, next_load, current_load, age, eggs, params):
    return((samp_id, _growth_fxn(next_load, current_load, age, eggs, params)))

def _worms_lost(worms_lost, current_load, params):
    """ 
    Worms lost due to worm death 

    Draw probabilities from binomial distribution
    """
    prob_worms_lost = stats.binom._pmf(worms_lost, 
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

def _age_dependent_survival(age, params):
    """
    Age-dependent survival probability
    """

    if age == 0:
        surv_prob = 1 - params['baby_death']
    else:
        surv_prob = np.exp(-params['weibull_lambda']*(age / 240.)**params['weibull_alpha'])

    return(surv_prob)

age_dependent_survival = np.vectorize(_age_dependent_survival)


def _load_dependent_survival(load, params):
    """
    Parasite load-dependent survival fxn
    """

    surv_prob = np.exp(params['death_thresh'] + params['patho']*np.log(load + 1)) / \
                        (1 + np.exp(params['death_thresh'] + params['patho']*np.log(load + 1)))
    return(surv_prob)


def nbinom_pmf(x, mu, k):
    """ Replace nan with 1...annoying that scipy doesn't do this"""

    vals = np.atleast_1d(stats.nbinom._pmf(x, n=k, p=(k / (mu + k))))
    vals[np.isnan(vals)] = 1 
    return(vals)

def full_R0(params, age_struct):
    """
    Calculate macroparasite R0 using approach of Diekmann

    Parameters
    ----------
    params : dict
        Parameters of the model
    age_struct : array-like
        Age structure of the raccoon population

    """

    # Two parasite life stages: 
    # 1. Eggs 
    # 2. Worms in raccoons

    para_stages = 2
    num_ages = len(age_struct)

    full_M = np.zeros((num_ages*para_stages, num_ages*para_stages))
    full_D = np.zeros((num_ages*para_stages, num_ages*para_stages))

    # Parasite eggs to Larvae matrix
    # TODO: Is it params['egg_contact']*racsofage or 1 -  np.exp(-params['egg_contact']*racsofage)?
    # Seem like the latter makes more sense.
    trans_probs = [_worm_infection_prob(age, 0, params)*(1 - np.exp(-params['egg_contact']*racsofage)) for racsofage, age in zip(age_struct, range(num_ages))]
    M_LE = np.diag(trans_probs)
    M_EL = np.diag(np.repeat(params['worm_repro'], num_ages))

    # Duration in stage matrix
    stay_egg = (np.exp(-params['egg_death_rate']) * (1 - np.sum(trans_probs)))
    D_EE = np.diag(np.repeat(1 - stay_egg, num_ages))

    stay_larv = np.array([survival_fxn(age, 0, params)*(np.exp(-params['worm_death_rate'])) for age in range(num_ages)])
    D_LL = np.diag(1 - stay_larv)

    full_M[num_ages:, :num_ages] = M_LE
    full_M[:num_ages, num_ages:] = M_EL
    full_D[:num_ages, :num_ages] = D_EE
    full_D[num_ages:, num_ages:] = D_LL

    # These two approaches are equivalent
    R0_mat = np.linalg.matrix_power(np.dot(full_M, np.linalg.inv(full_D)), 2)
    R0_mat_other = np.dot(np.dot(M_LE, np.linalg.inv(D_EE)), np.dot(M_EL, np.linalg.inv(D_LL)))

    return((R0_mat, R0_mat_other))

    

def relative_R0(params, max_age=100):
    """
    Calculate the relative R0 for model
    """

    R0s = []
    for age in range(max_age):

        host_surv_age_para = survival_fxn(age, 1, params)
        para_surv = np.exp(-params['worm_death_rate']) # Worm survival over 1 week
        prob_inf = _worm_infection_prob(age, 0, params)

        # TODO: This is not yet correct...still working on it
        term1 = params['worm_repro'] / (1 - (host_surv_age_para*para_surv))
        term2 = 1 / (1 - np.exp(-params['egg_death_rate'])*np.exp(-params['egg_contact']))
        term3 = np.exp(-params['egg_contact']) * prob_inf

        rel_R0 = term1 * term2 * term3
        R0s.append(rel_R0)

    return(R0s)

