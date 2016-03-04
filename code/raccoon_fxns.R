## Functions for Raccoon Simulation

death_probability = function(load, beta, alpha){
    # Calculate death probability given some worm load
    # PIHM.  This is calculated as 1 - survival probability

    # Parameters
    # -----------
    # load : int
    #     Worm load >= 0
    # alpha : float
    #     Pathogenicity
    # beta : float
    #   "Threshold" above which mortality is likely. Depends on alpha
    #
    # Returns
    # -------
    # : float
    #    Probability of death

    # TODO: IF LOAD IS 0 THROW AN ERROR.


    prob_death = 1 - exp(beta + log(load + 1)*alpha) / (1 + exp(beta + log(load + 1)*alpha))

    #prob_death = beta + load * alpha

    return(prob_death)
    #return(0)
}

kill_my_raccoon = function(worm_load, death_thresh, patho){
    # Kill raccoon based on worm load and intrinsic mortality

    tr = runif(1)
    alive_now = tr > death_probability(worm_load, death_thresh, patho)

    return(alive_now)

}

kill_raccoon_worms = function(previous_chorts, death_thresh, death_slope){
    # Function to kill worms in raccoon based on worm age
    # Assuming a logistic function of survival prob of worms with cohort
    # age and then killing the worms based on a draw from a binomial

    ages = length(previous_cohorts):1
    surv_probs = 1 / (1 + exp(-1*(death_thresh + death_slope * ages)))
    new_cohort = rbinom(length(previous_cohorts), previous_cohorts,
                                                 surv_probs)

}

give_birth = function(age_now, time,
                        month_at_repro,
                        first_repro_age, litter_size){
    # Decide how many babies are produced by a raccoon. Check if the month
    # is right and if the raccoon is old enough.

    repro = 0

    if(time %% month_at_repro == 0){ # If the month is right

        if(age_now >= first_repro_age){ # If the age is right

            repro = litter_size

        }

    }

    return(repro)

}


pick_up_eggs = function(eprob, emean, infect, resist, load){
    # Function to pick up eggs. Depends on eprob (encounter_probability),
    # emean (mean number of eggs contacted), infect (infectivity)

    # Exponential decline of infectivity.
    infect_red = infect * exp(-resist * load)

    return(round(eprob * rpois(1, emean) * infect_red))

}

update_arrays = function(time, new_babies, new_babies_vect,
                                           initial_age_vector,
                                           raccoon_dead_alive_array,
                                           raccoon_worm_array,
                                           age_array){

    # Function takes in the various arrays used in the raccoon simulation
    # and updates them based on the new babies that were born
    #
    # Returns
    # -------
    # : a list with all of the updated arrays

    new_babies_vect[time] = new_babies

    ## Extend the arrays to account for births ##
    initial_age_vector = c(initial_age_vector, rep(0, new_babies))

    new_alive_babies = array(NA, dim=c(TIME_STEPS + 1, new_babies))
    new_worm_babies = array(NA, dim=c(TIME_STEPS + 1, new_babies))
    new_age_babies = array(NA, dim=c(TIME_STEPS + 1, new_babies))

    # # Set the baby arrays
    new_alive_babies[time, ] = 1 # All new babies are alive
    new_worm_babies[time, ] = 0 # All new babies have 0 worms
    new_age_babies[time, ] = 0 # All new babies have age 0

    age_array = cbind(age_array, new_age_babies)

    raccoon_dead_alive_array = cbind(raccoon_dead_alive_array,
                                                new_alive_babies)

    raccoon_worm_array = cbind(raccoon_worm_array, new_worm_babies)

    return(list(new_babies_vect=new_babies_vect,
                initial_age_vector=initial_age_vector,
                age_array=age_array,
                raccoon_dead_alive_array=raccoon_dead_alive_array,
                raccoon_worm_array=raccoon_worm_array))

}