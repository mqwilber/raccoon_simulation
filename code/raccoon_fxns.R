## Functions for Raccoon Simulation

library(ggplot2)
library(reshape2)
library(data.table)

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

intrinsic_death_fxn = function(age, intrinsic_death_rate, baby_death){

    # TODO: UPDATE THE AGE_DEPENDENT SURVIVAL FUNCTION
    return(baby_death * exp(-intrinsic_death_rate * age))
}

kill_my_raccoon = function(worm_load, age, death_thresh, patho, intrinsic_death,
                                baby_death, random_death_prob){
    # Kill raccoon based on worm load and intrinsic mortality

    tr = runif(1)

    # Probaility of not dying from worms * not dying from age * not dying
    # from random
    surv_prob = (1 - death_probability(worm_load, death_thresh, patho)) *
                        (1 - intrinsic_death_fxn(age, intrinsic_death,
                            baby_death)) * (1 - random_death_prob)
    alive_now = tr > (1 - surv_prob)

    return(alive_now)

}

kill_raccoon_worms = function(previous_chorts, death_thresh, death_slope){
    # Function to kill worms in raccoon based on worm age
    # Assuming a logistic function of survival prob of worms with cohort
    # age and then killing the worms based on a draw from a binomial

    updated_cohort = previous_cohorts
    non_na_ind = !is.na(previous_cohorts) # get non-nas

    ages = length(previous_cohorts):1
    surv_probs = 1 / (1 + exp(-1*(death_thresh + death_slope * ages)))

    new_cohort = rbinom(sum(non_na_ind), previous_cohorts[non_na_ind],
                                                 surv_probs[non_na_ind])
    updated_cohort[non_na_ind] = new_cohort
    return(updated_cohort)

}

give_birth = function(age_now, time, tot_racs,
                        month_at_repro,
                        first_repro_age, litter_size, beta){
    # Decide how many babies are produced by a raccoon. Check if the month
    # is right and if the raccoon is old enough. Then chose the

    repro = 0

    if(time %% month_at_repro == 0){ # If the month is right

        if(age_now >= first_repro_age){ # If the age is right

            repro_prob = exp(-(beta*tot_racs)) # Think about this beta function
            repro = rbinom(1, litter_size, repro_prob)
        }

    }

    return(repro)

}


pick_up_eggs = function(emean, ek, infect, resist, prev, load, egg_decay,
                            eprob_param){
    # Function to pick up eggs. Depends on eprob (encounter_probability),
    # emean (mean number of eggs contacted),
    # ek: negative binomial k
    # infect: infectivity
    # resit: Acquired immunity
    # load: worm load at time t - 1 for a given rac
    # prev: Prevalence of worms in pop for all past time points
    # egg_decay: Egg decay rate
    # eprob_param: parameter for determining encounter probability

    # Exponential decline of infectivity.
    infect_red = infect * exp(-resist * load)

    # Encounter and get eggs on face and get infected with eggs
    eprob = get_eprob(prev, egg_decay, eprob_param)
    print(eprob)
    # eprob = 1
    new_eggs = rbinom(1, 1, eprob) * rbinom(1, rnbinom(1, size=ek, mu=emean), infect_red)
    return(new_eggs)

}

get_eprob = function(prev_vector, egg_decay, eprob_param){
    # Calculating our encounter probability from previous prevalences

    # Calculate egg decay probability
    weights = exp(-egg_decay * 1:length(prev_vector))

    adjusted_prev = weights[length(prev_vector):1] * prev_vector
    metric = sum(adjusted_prev)

    # if there is no past infection eprob = 0
    if(metric == 0){
        return(0)
    } else{
        metric = log(metric)
        eprob = exp(eprob_param[1] + eprob_param[2] * metric) /
                        (1 + exp(eprob_param[1] + eprob_param[2] * metric))
        return(eprob)
    }


}

assign_human_contacts = function(num_racs){
    # Assign probabilities of human contacts
    return(runif(num_racs))
}

update_arrays = function(time, new_babies, new_babies_vect,
                                           initial_age_vector,
                                           raccoon_dead_alive_array,
                                           raccoon_worm_array,
                                           age_array, infra_worm_array,
                                           human_array,
                                           babies_at_this_time_vect,
                                           time_steps){

    # Function takes in the various arrays used in the raccoon simulation
    # and updates them based on the new babies that were born
    #
    # Returns
    # -------
    # : a list with all of the updated arrays

    new_babies_vect[time] = new_babies

    ## Extend the arrays to account for births ##
    initial_age_vector = c(initial_age_vector, rep(0, new_babies))

    new_alive_babies = array(NA, dim=c(time_steps + 1, new_babies))
    new_worm_babies = array(NA, dim=c(time_steps + 1, new_babies))
    new_age_babies = array(NA, dim=c(time_steps + 1, new_babies))

    # # Set the baby arrays
    new_alive_babies[time, ] = 1 # All new babies are alive
    new_worm_babies[time, ] = 0 # All new babies have 0 worms
    new_age_babies[time, ] = 0 # All new babies have age 0

    age_array = cbind(age_array, new_age_babies)

    # Add new babies behaviors onto human array. Only for the mothers that
    # are alive (i.e. not NAs)
    ind = !is.na(babies_at_this_time_vect)
    human_array = c(human_array, rep(human_array[ind],
                                    babies_at_this_time_vect[ind]))

    raccoon_dead_alive_array = cbind(raccoon_dead_alive_array,
                                                new_alive_babies)

    raccoon_worm_array = cbind(raccoon_worm_array, new_worm_babies)

    # Extend the worm_infra_array
    current_racs = length(infra_worm_array)
    for(k in 1:new_babies){
        infra_worm_array[[k + current_racs]] =
                            array(NA, dim=c(time_steps + 1, time_steps + 1))

        # Assign 0 worms to new racs
        infra_worm_array[[k + current_racs]][time, time] = 0
    }

    return(list(new_babies_vect=new_babies_vect,
                initial_age_vector=initial_age_vector,
                age_array=age_array,
                raccoon_dead_alive_array=raccoon_dead_alive_array,
                raccoon_worm_array=raccoon_worm_array,
                infra_worm_array=infra_worm_array,
                human_array=human_array))

}

## Summary functions ##

get_prevalence = function(raccoon_worm_array){
    # Get prevalence for worm array
    # Input : the raccoon worm array
    # Returns : prevalence for all time points

    prev_fxn = function(x) {
        return(sum(x > 0, na.rm=T) / sum(!is.na(x)))
    }
    return(apply(raccoon_worm_array, 1, prev_fxn))
}


## Functions to summarize simulation output ##

plot_pop_traj = function(raccoon_dead_alive_array){
    # plot population trajectories
    pop_trun = rowSums(raccoon_dead_alive_array, na.rm=T)
    dat = data.frame(list(population=pop_trun, time=1:length(pop_trun)))
    ggplot(dat, aes(x=time, y=population)) + geom_point() + geom_line() +
        ylim(c(0, max(dat$population)))
}

plot_age_hist = function(age_array, range){

    # Plot the raccoon age distribution at any time point with violin plots

    # age_array: age array from simulation
    # range: indices of the time points at which to make histograms of age
    # distributions. i.e. 30:45

    trun_age = data.frame(t(age_array[range, ]))
    names(trun_age) = range
    stacked_data = stack(trun_age)
    stacked_data_trun = stacked_data[!is.na(stacked_data$values), ]
    ggplot(stacked_data_trun, aes(as.factor(ind), values)) + geom_violin()

}

plot_worm_dist = function(raccoon_worm_array, range){

    # Plot the raccoon age distribution at any time point with violin plots

    # raccoon_worm_array: age array from simulation
    # range: indices of the time points at which to make histograms of age
    # distributions. i.e. 30:45

    trun_age = data.frame(t(raccoon_worm_array[range, ]))
    names(trun_age) = range
    stacked_data = stack(trun_age)
    stacked_data_trun = stacked_data[!is.na(stacked_data$values), ]
    ggplot(stacked_data_trun, aes(as.factor(ind), values)) + geom_boxplot()

}

age_intensity_full = function(raccoon_worm_array, age_array, range){
    # Plot full age intensity profile
    # range specifies which time range to calculate the age-intensity for

    flat_age = as.vector(age_array[range, ])
    flat_worms = as.vector(raccoon_worm_array[range, ])
    full_dat = as.data.table(data.frame(list(age=flat_age, worms=flat_worms)))
    means = full_dat[, list(mean_worms=mean(worms, na.rm=T)), by="age"]
    means = means[!is.na(means$age), ]
    ggplot(means, aes(x=age, y=mean_worms)) + geom_point() + geom_line()
}

worm_traj = function(raccoon_worm_array){
    # Plot trajectories of individual raccoons

    sdata = stack(data.frame(raccoon_worm_array))
    sdata$time = rep(1:dim(raccoon_worm_array)[1], dim(raccoon_worm_array)[2])
    ggplot(sdata, aes(x=time, y=values, color=ind)) + geom_line() + geom_point()
}