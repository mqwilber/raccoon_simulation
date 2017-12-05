from __future__ import division

from numpy.testing import (TestCase, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, assert_, assert_raises)

import numpy as np 
from raccoon_ipm import *
from raccoon_ipm import _worms_gained, _worm_infection_prob, _growth_fxn
from scipy.stats import nbinom
import yaml


"""
Dsecription
------------
Unittests for ipm raccoon methods

"""

with open('ipm_params.yml') as f:
    # use safe_load instead load
    def_params = yaml.safe_load(f)

class TestIPM(TestCase):

    def test_baby_survival(self):

        # Test correct baby survival is returned
        baby_surv = survival_fxn(0, 0, def_params)
        assert_almost_equal(baby_surv, (1 - def_params['baby_death'])* (1 - def_params['random_death']))


    def test_adult_survival(self):

        # Adults at 20 years should die
        adult_surv = survival_fxn(12*20, 0, def_params)
        assert_equal(1 - adult_surv, 1)


    def test_worms_gained(self):

        new_params = def_params.copy()
        new_params['egg_contact'] = 0

        # Prob should be 0 if there is no egg contact
        prob1 = _worms_gained(5, 10000, 0, 5, new_params)
        assert_equal(prob1, 0)

        # Prob should be non 0 if they are above age_resistance
        age = def_params['age_resistance'] + 1
        prob2 = _worms_gained(5, 10000, age, 5, new_params)
        assert_equal(prob2 != 0, True)

        # Check that prob is only from rodent
        new_params['age_immunity'] = 100000 # Set arbitrarily high
        prob3 = _worms_gained(5, 10000, age, 5, new_params)
        true_prob = nbinom_pmf(5, mu=new_params['mean_rodent']*new_params['larval_infect']*new_params['rodent_contact'], 
                                  k=new_params['k_rodent'])
        assert_equal(prob3, true_prob)
    
    def test_worm_infection_prob(self):

        # Check the infection probability with no load
        prob_inf = _worm_infection_prob(0, 0, def_params)
        assert_equal(prob_inf, def_params['infect'])

    def test_growth_fxn(self):

        # Check that loads sum to 1
        marginal_prob = np.sum([_growth_fxn(x, 5, 0, 10000, def_params) for x in range(200)])
        assert_almost_equal(marginal_prob, 1)