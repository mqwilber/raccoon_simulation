## Functions for Raccoon Simulation

death_probability = function(load, beta, alpha){
    # Calculate death probability given some worm load
    # PIHM

    # Parameters
    # -----------
    # load : int
    #     Worm load >= 0
    # alpha : float
    #     Pathogenicity
    # beta : float
    #   Death probability in absence of parasites. Between 0 and 1
    #
    # Returns
    # -------
    # : float
    #    Probability of death

    # TODO: IF LOAD IS 0 THROW AN ERROR
    # TODO: If beta is not a probability throw and ERROR

    # CHANGE TO LOGISTIC
    prob_death = beta + load * alpha

    return(prob_death)
}

kill_my_raccoon = function(worm_load, death_prob, patho){
    # Kill raccoon based on worm load and intrinsic mortality

    tr = runif(1)
    alive_now = tr > death_probability(worm_load, death_prob, patho)

    return(alive_now)

}


pick_up_eggs = function(eprob, emean, infect){
    # Function to pick up eggs. Depends on eprob (encounter_probability),
    # emean (mean number of eggs contacted), infect (infectivity)

    return(round(eprob * rpois(1, emean) * infect))

}