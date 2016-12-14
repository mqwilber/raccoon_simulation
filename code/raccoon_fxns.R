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

    # TODO: REVISIT THIS.  ARE WE DOUBLE COUNTING DEATH
    return(baby_death * exp(-intrinsic_death_rate * age))
}

senescence_fxn = function(age, old_death){
    # Probability of dying of old age
    prob_sen = min(c(1, old_death * age^2))
    return(prob_sen)
} 

kill_my_raccoon = function(worm_load, age, overlap, death_thresh, patho, intrinsic_death,
                                baby_death, random_death_prob, old_death,
                                cull=FALSE){
    # Kill raccoon based on worm load and intrinsic mortality
    # 
    # Parameters
    # ----------
    # worm_load : int, worm load for a given raccoon
    # age : int, age of given raccoon
    # overlap : float, between 0 and 1. Human overlap of a given raccoon 
    #            where 1 is high overlap
    # ... : parameters defined in `raccoon_parameters.R`
    # cull: bool, TRUE = cull raccoon, FALSE = don't kill

    tr = runif(1)

    # Probaility of not dying from worms * not dying from age * not dying
    # from random
    surv_prob = (1 - death_probability(worm_load, death_thresh, patho)) *
                        (1 - intrinsic_death_fxn(age, intrinsic_death,
                            baby_death)) * (1 - random_death_prob) * 
                            (1 - senescence_fxn(age, old_death))
                            #(1 - cull_strategy(cull_params, age, overlap))

    alive_now = (tr > (1 - surv_prob)) * (!cull)

    return(alive_now)

}

get_cull_indices = function(cull_params, raccoon_dead_alive_vect, 
                                age_vect, human_risk_through_time_vect){
    # Function for updating cull_params each time step to get per capita
    # raccoon survival prob based on quota of raccoon to kill each time point
    # 
    # Possible strategies are:
    # 1. `random`: Cull everybody equally
    #   - Additional parameters: `cull_prob`: probability of being culled
    # 2. `age`: Only cull juveniles i.e. less than a year old
    #   - Additional parameters: `cull_prob`: probability of being culled
    # 3. `human`: Only cull individuals with a human overlap greater than 
    #              `overlap_threshold`
    #   - Additional parameters: `quota`: absolute number of raccoons to kill
    #                            `overlap_threshold`: Between 0 and 1,
    #                             threshold above which you are culled
    #
    # Parameters
    # ----------
    # cull_params : list or NULL, cull params, see cull_strategy
    # raccoon_dead_alive_vect: vector, vector specifying whether raccoons 
    #      are dead or alive at time t - 1 
    # age_vect : vector, specifying the ages of alive raccoons at year time
    #                   t - 1
    # human_risk_through_time_vect: vector, specifies the human overlap value
    #                   of alive raccoons at time t - 1.
    # Returns
    # --------
    # : indices of raccoons to cull
    #


    if(is.null(cull_params)){
        return(integer(0)) # return an empty vector
    }

    if(cull_params$strategy == "random"){

        pop_inds = which((raccoon_dead_alive_vect != 0) & 
                                        (!is.na(raccoon_dead_alive_vect)))

    } else if(cull_params$strategy == "age"){

        pop_inds = which((age_vect < 12) & (!is.na(age_vect)))
      
    } else if(cull_params$strategy == "human"){

        if(is.null(cull_params$overlap_threshold)){
            stop("Provide overlap_threshold")
        }

        pop_inds = which(human_risk_through_time_vect > cull_params$overlap_threshold)

    } else{
        stop(paste(cull_params$strategy, "is not a recognized strategy. Try random, age, or human"))
    }


    # Check if there are any raccoons to cull.  Only cull quota or pop_inds
    if(length(pop_inds) > 1){
        cull_inds = sample(pop_inds, min(c(length(pop_inds), cull_params$quota)))
    } else{
        cull_inds = pop_inds # Empty vector
    }

    return(cull_inds)

}


kill_raccoon_worms = function(previous_cohorts, death_thresh, death_slope, 
                got_bait=0){
    # Function to kill worms in raccoon based on worm age
    # Assuming a logistic function of survival prob of worms with cohort
    # age and then killing the worms based on a draw from a binomial

    updated_cohort = previous_cohorts
    non_na_ind = !is.na(previous_cohorts) # get non-nas

    ages = length(previous_cohorts):1
    surv_probs = 1 / (1 + exp(-1*(death_thresh + death_slope * ages)))

    new_cohort = rbinom(sum(non_na_ind), previous_cohorts[non_na_ind],
                                                 surv_probs[non_na_ind])

    # Kill worms in cohort based on getting bait
    new_cohort = new_cohort * (1 - got_bait)

    updated_cohort[non_na_ind] = new_cohort
    return(updated_cohort)

}

picked_up_bait = function(overlap, worm_control_params=NULL){
    # Function which determines whether or not a raccoon picks up 
    # anti-helminthic bait. worm_control_params is a list that contains the 
    # necessary parameters
    # Worm_control_params must contain a slot 'strategy' that specifies the
    # strategy that is being used to bait. If NULL, no bating is done.
    #
    # Possible strategies are:
    # 1. `random`: Bait everybody equally
    #   - Additional parameters: `distribution`: distribution of bait in the 
    #                                        environment between 0 and 1
    # 2. `human`: Only bait individuals with a human overlap greater than 
    #              `overlap_threshold`
    #   - Additional parameters: `distribution`: Distribution of bait between 0 and 1
    #                            `overlap_threshold`: Between 0 and 1,
    #                             threshold above which you are baited
    #


     # Clear worms based on strategy
    if(!is.null(worm_control_params)){

        if(worm_control_params$strategy == "random"){

            got_bait = rbinom(1, 1, worm_control_params$distribution)

        } else if(worm_control_params$strategy == "human"){

            if(is.null(worm_control_params$overlap_threshold)){
                stop("Provide overlap_threshold")
            }

            got_bait = ifelse(overlap > worm_control_params$overlap_threshold, 
                    rbinom(1, 1, worm_control_params$distribution), 0)

        } else {
            stop(paste(worm_control_params$strategy, "is not a recognized strategy. Try random or human"))

        }
    } else{
        got_bait = 0
    }

    return(got_bait)

}

give_birth = function(age_now, time, tot_racs,
                        month_at_repro,
                        first_repro_age, litter_size, beta,
                        birth_control_params=NULL){
    # Decide how many babies are produced by a raccoon. Check if the month
    # is right and if the raccoon is old enough. Then chose the

    repro = 0

    if(time %% month_at_repro == 0){ # If the month is right

        if(age_now >= first_repro_age){ # If the age is right

            birth_event = birth_control_strategy(birth_control_params)
            repro_prob = birth_event * exp(-(beta*tot_racs)) # Think about this beta function
            repro = rbinom(1, litter_size, repro_prob)
        }

    }

    return(repro)

}

birth_control_strategy = function(birth_control_params){
    # Strategies for birth control. `birth_control_params` must have the 
    # name `strategy` which specifies which birth control strategy to use.
    # 
    # 
    # Possible strategies are:
    # 1. `random` : All individuals are susceptible to birth control
    #   - Additional parameters: `distributions` - between 0 and 1 where
    #   0 is no distribution of birth control and 1 is complete distribution
    #   of birth control.
    #
    # Parameters
    # ----------
    # birth_control_params : list, see above
    #
    # Returns
    # --------
    # 0 or 1 determining whether or not it enter breeding cycle

    if(is.null(birth_control_params$strategy)){
        return(1)
    }

    if(birth_control_params$strategy == "random"){
        birth_event = rbinom(1, 1, 1 - birth_control_params$distribution)
    }

    return(birth_event)
}



pick_up_eggs = function(emean, ek, infect, resist, prev, load, egg_decay,
                            eprob_param){
    # Function to pick up eggs. Depends on eprob (encounter_probability),
    # emean (mean number of eggs contacted),
    # ek: negative binomial k
    # infect: infectivity
    # resist: Acquired immunity
    # load: worm load at time t - 1 for a given rac
    # prev: Prevalence of worms in pop for all past time points
    # egg_decay: Egg decay rate
    # eprob_param: parameter for determining encounter probability


    # Exponential decline of infectivity.
    infect_red = infect * exp(-resist * load)

    # Encounter and get eggs on face and get infected with eggs
    eprob = get_eprob(prev, egg_decay, eprob_param)
    # eprob = 1
    new_eggs = rbinom(1, 1, eprob) * rbinom(1, rnbinom(1, size=ek, mu=emean), infect_red)
    return(new_eggs)

}

pick_up_rodents = function(mouse_worm_mean, mouse_worm_agg, 
                        rodent_encounter_prob, larval_worm_infectivity,
                        prev_vector, egg_decay, eprob_param){
    # Function to pick up rodents from rodent pool.  Rodent 
    # encounters worm, picks up worms from a negative binomial, and 
    # then worms establish with some prob
    #
    # mouse_worm_mean : mean number of larval worms acquired
    # mouse_worm_agg : aggregation of larval worms
    # rodent_encounter_prob : probability of raccoon encountering rodent
    # larval_worm_infectivity : Probability of larval worm establishing
    #
    # Returns
    # -------
    # : number of larval worms acquired

    # Update rodent mean and variance
    mouse_mean_k = update_rodent_mean_var(prev_vector, egg_decay, eprob_param, 
                                mouse_worm_mean, mouse_worm_agg)

    larval_worms = rnbinom(1, mu=mouse_mean_k$mean, size=mouse_mean_k$k) # Larval worms per mouse


    new_worms = rbinom(1, 1, rodent_encounter_prob) * # Encounter with mice
                rbinom(1, larval_worms, larval_worm_infectivity) # Worms establishing

    return(new_worms)

}

update_rodent_mean_var = function(prev_vector, egg_decay, eprob_param, 
                                mouse_worm_mean, mouse_worm_agg){

    eprob = get_eprob(prev_vector, egg_decay, eprob_param)

    new_mean = eprob * mouse_worm_mean

    # From Shaw and Dobson 1995 TPL
    # TODO: Change intercept to match empirically observed k/variance value from
    # rodents.
    logvar = 1.551*log10(new_mean) + 1.098

    # Convert to NBD k
    new_k = (new_mean)^2 / (10^logvar - new_mean)

    return(list(mean=new_mean, k=new_k))

}

get_eprob = function(prev_vector, egg_decay, eprob_param){
    # Calculating our encounter probability from previous prevalences

    # Calculate egg decay probability
    # TODO: REVIST WEIGHTS SHOULD THIS START AT 0? (i.e. 0:(length(prev_vector) - 1))
    weights = exp(-egg_decay * 1:length(prev_vector))

    # Weight the the past infections prevalences and then sum them.
    # In general, when this metric is 1 or above we'd expect a high
    # probability of encounter probability
    adjusted_prev = weights[length(prev_vector):1] * prev_vector
    metric = sum(adjusted_prev)

    # if there is no past infection eprob = 0
    if(metric == 0){
        return(0)
    } else{ # The metric is not zero, it determines infection prob
            # according to a logistic function.

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
                                           age_array, all_worms_infra_array,
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

    for(i in 1:length(all_worms_infra_array)){

        current_racs = length(all_worms_infra_array[[i]])
        for(k in 1:new_babies){
            all_worms_infra_array[[i]][[k + current_racs]] =
                                array(NA, dim=c(time_steps + 1, time_steps + 1))

            # Assign 0 worms to new racs
            all_worms_infra_array[[i]][[k + current_racs]][time, time] = 0
        }
    }

    return(list(new_babies_vect=new_babies_vect,
                initial_age_vector=initial_age_vector,
                age_array=age_array,
                raccoon_dead_alive_array=raccoon_dead_alive_array,
                raccoon_worm_array=raccoon_worm_array,
                all_worms_infra_array=all_worms_infra_array,
                human_array=human_array))

}


combine_worm_arrays = function(all_worm_arrays, raccoon_dead_alive_array){
    # Takes in multiple independent worm arrays and combines them into a 
    # total worm array. Used for distinguishing between worms from rodents
    # and worms from eggs.
    #
    # all_worm_arrays : list of different infra_worm_arrays from different
    #                   possible sources. Could be from arbitrarily many sources
    #

    # Create matrices in lists
    convert_lists = lapply(seq_along(all_worm_arrays), 
                        function(x) simplify2array(all_worm_arrays[[x]]))

    # 4 dimensional matrix
    converted_matrix = simplify2array(convert_lists)

    # Sum across different worm arrays
    sum_worms = apply(converted_matrix, 1:3, sum, na.rm=TRUE)

    # Convert back to lists
    sum_worms_list = lapply(1:dim(sum_worms)[3], function(i) sum_worms[, , i])

    # Replace upper diags with NAs
    sum_worms_list = lapply(sum_worms_list, function(x){x[upper.tri(x)] = NA; return(x)})

    # Replace 0s with NA if raccoon did not yet exist
    for(j in 1:length(sum_worms_list)){

        # If raccooon doesn't yet exist set to NA
        inds = which(is.na(raccoon_dead_alive_array[, j]))
        sum_worms_list[[j]][inds, ] = NA
        sum_worms_list[[j]][, inds] = NA

        # If raccoon is dead set to NA
        dead = which(raccoon_dead_alive_array[, j] == 0)
        sum_worms_list[[j]][dead, ] = NA

    }

    return(sum_worms_list)


}

get_tot_worm_array_from_infra = function(infra_worm_array, raccoon_dead_alive_array){
    # Convert infra-worm array to total raccoon worm array

    total_worm_array = do.call(cbind, lapply(infra_worm_array, 
                                function(x) apply(x, 1, sum, na.rm=T)))
    total_worm_array[raccoon_dead_alive_array == 0] = NA
    total_worm_array[is.na(raccoon_dead_alive_array)] = NA

    return(total_worm_array)

}

get_human_risk_metric = function(human_risk_through_time,
                                        raccoon_worm_array,
                                        egg_decay){
    # Calculating the human risk metric through time
    # This is using a weighted measure of human risk based on the past
    # worm loads in the population
    #
    # Parameters
    # human_risk_through_time: list containing the vectors that have human
    #                           risks for each raccoon
    # raccoon_worm_array: Raccoon worm array
    # egg_decay: parameter determining egg decay rate in the environment

    # Number of raccoons at end of simulation
    tot_rac = length(human_risk_through_time[[length(human_risk_through_time)]])

    # Function to grow all arrays to same size
    growth_array = function(x, size){
        num_grow = size - length(x)
        return(c(x, array(NA, dim=num_grow)))
    }

    # Matrix of humans risks
    new_risk = lapply(human_risk_through_time, growth_array, tot_rac)
    new_risk = do.call(rbind, new_risk)

    tot_time = length(human_risk_through_time)

    overlap_array = raccoon_worm_array * new_risk
    unweighted_risk_array = rowSums(overlap_array, na.rm=T)

    weighted_risk_array = array(NA, dim=length(unweighted_risk_array))

    for(time in 1:tot_time){

        vals = 0:(time - 1)
        weights = exp(-egg_decay * vals)
        weighted_hr = sum(unweighted_risk_array[1:time] * weights[time:1])
        weighted_risk_array[time] = weighted_hr

    }

    return(list(weighted=weighted_risk_array, unweighted=unweighted_risk_array))

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
    ggplot(stacked_data_trun, aes(as.factor(ind), values)) + geom_boxplot()

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

age_prevalence_plot = function(raccoon_worm_array, age_array, range){
    # Plot full age intensity profile
    # range specifies which time range to calculate the age-intensity for

    flat_age = as.vector(age_array[range, ])
    flat_worms = as.vector(raccoon_worm_array[range, ])
    full_dat = as.data.table(data.frame(list(age=flat_age, worms=flat_worms)))

    prev_fxn = function(x){
        no_na_x = x[!is.na(x)]
        return(sum(no_na_x != 0) / length(no_na_x))
    }

    prevs = full_dat[, list(prev_worms=prev_fxn(worms)), by="age"]
    prevs = prevs[!is.na(prevs$age), ]
    ggplot(prevs, aes(x=age, y=prev_worms)) + geom_point() + geom_line()
}

age_intensity_plot = function(raccoon_worm_array, age_array, range){
    # Plot full age intensity profile
    # range specifies which time range to calculate the age-intensity for

    flat_age = as.vector(age_array[range, ])
    flat_worms = as.vector(raccoon_worm_array[range, ])
    full_dat = as.data.table(data.frame(list(age=flat_age, worms=flat_worms)))
    ggplot(full_dat, aes(x=age, y=worms)) + geom_point() + geom_smooth()

}


worm_traj = function(raccoon_worm_array){
    # Plot trajectories of individual raccoons

    sdata = stack(data.frame(raccoon_worm_array))
    sdata$time = rep(1:dim(raccoon_worm_array)[1], dim(raccoon_worm_array)[2])
    tplot = ggplot(sdata, aes(x=time, y=values, color=ind)) + 
                    geom_line() + geom_point() +
                    theme(legend.position="none")
    return(tplot)
}


## Get parameters and arrays ##

get_simulation_parameters = function(...){
    # Get parameters for the simulation
    #
    # You can pass keyword arguments to this function to reassign the 
    # default parameters

    ## Raccoon parameters

    INIT_NUM_RACCOONS = 500
    #DEATH_PROB = 22 # Natural death
    DEATH_THRESHOLD = 22 #22 # beta in death_probability fxn
    PATHOGENICITY = -4.2 # alpha in death_probability fxn
    BABY_DEATH = 1 - (.52^(1/7))^4 # Probability of dying as baby
    INTRINSIC_DEATH_RATE = 0.33 # Age related death rate
    RANDOM_DEATH_PROB = 0.01 # Lower bound to death prob
    OLD_DEATH = (1 / (20 * 12)^2) # Above 20 years old the raccoon dies
    AGE_EGG_RESISTANCE = 4 # Age above which raccoons no longer pick up eggs
    RODENT_ENCOUNTER_PROB = 0.5 # Monthly probability of encountering a rodent
    INIT_WORMS = 10 # Initial number of non-rodent worms in raccoons

    ## Rodent parameters
    MOUSE_WORM_MEAN = 3.49 # Abundance of worms in peromyscus estimated from Sara's data
    MOUSE_WORM_AGG = 0.22 # Aggregation of worms (k parameter) in peromyscus
    LARVAL_WORM_INFECTIVITY = 0.25 # Probability of larval worm establishing


    ## Age parameters
    FIRST_REPRO_AGE = 10 # months
    LITTER_SIZE = 4 # Maximum litter size
    MONTH_AT_REPRO = 12
    DISPERSAL_AGE = 12 # months

    ## Ricker function for density-dependent recruitment of new babies
    K_CAPACITY = 150 # "Carrying" capacity for raccoons. Need to figure out what
                    # what is determining carrying capacity for the deterministic
                    # version of this model
    DI_NEW_BABY_DEATH_RATE = 0 # Density independent new baby death
    BIRTH_RATE = log(LITTER_SIZE) # Gives birth rate. Little r in ricker function
    BETA = BIRTH_RATE / K_CAPACITY # TODO: CHECK WHY I DID THIS?

    # PARASITE PARAMETERS
    # ENVIRONMENTAL_POOL = 0
    # EGG_PRODUCTION_PER_WORM = 100 # eggs per month (note this is very wrong)
    ENCOUNTER_MEAN = 500 # Mean number of encountered eggs
    ENCOUNTER_K = 1 # Aggregation parameter of the NBD

    # A set of parameters that determines the
    # probability of encounter given a weighted history of prevalence. The first
    # parameter dictates at what level of the log metric eprob is close to 1 and
    # the second parameter dictates how quickly eprob goes to 1.
    # TODO: WE WILL NEED TO FIDDLE WITH THESE
    ENCOUNTER_PARAMS = c(2, 5)

    INFECTIVITY = 0.02 # Probability of infectivity
    RESISTANCE = 0.03 # How quickly a raccoon gains resistance based on previous load
    EGG_DECAY = 0.3 # Rate of egg decay such that 2%-3% chance of survival after year. TODO: Revisit this

    # See fit_param.R for how we got these values
    WORM_SURV_TRESH = 4.7104 #/ 2 # Threshold parameter of worm survival probability
    WORM_SURV_SLOPE = -0.9446 #/ 2 # Slope of worm surv probability

    # Time parameters: Each time step is a month
    TIME_STEPS = 100

    params=list(INIT_NUM_RACCOONS=INIT_NUM_RACCOONS,
                DEATH_THRESHOLD=DEATH_THRESHOLD,
                PATHOGENICITY=PATHOGENICITY,
                BABY_DEATH=BABY_DEATH,
                INTRINSIC_DEATH_RATE=INTRINSIC_DEATH_RATE,
                RANDOM_DEATH_PROB=RANDOM_DEATH_PROB,
                OLD_DEATH=OLD_DEATH,
                AGE_EGG_RESISTANCE=AGE_EGG_RESISTANCE,
                RODENT_ENCOUNTER_PROB=RODENT_ENCOUNTER_PROB,
                INIT_WORMS=INIT_WORMS,
                MOUSE_WORM_MEAN=MOUSE_WORM_MEAN,
                MOUSE_WORM_AGG=MOUSE_WORM_AGG,
                LARVAL_WORM_INFECTIVITY=LARVAL_WORM_INFECTIVITY,
                FIRST_REPRO_AGE=FIRST_REPRO_AGE,
                LITTER_SIZE=LITTER_SIZE,
                MONTH_AT_REPRO=MONTH_AT_REPRO,
                DISPERSAL_AGE=DISPERSAL_AGE,
                K_CAPACITY=K_CAPACITY,
                DI_NEW_BABY_DEATH_RATE=DI_NEW_BABY_DEATH_RATE,
                BIRTH_RATE=BIRTH_RATE,
                BETA=BETA,
                ENCOUNTER_MEAN=ENCOUNTER_MEAN,
                ENCOUNTER_K=ENCOUNTER_K,
                ENCOUNTER_PARAMS=ENCOUNTER_PARAMS,
                INFECTIVITY=INFECTIVITY,
                RESISTANCE=RESISTANCE,
                EGG_DECAY=EGG_DECAY,
                WORM_SURV_TRESH=WORM_SURV_TRESH,
                WORM_SURV_SLOPE=WORM_SURV_SLOPE,
                TIME_STEPS=TIME_STEPS)

    new_params = list(...)

    # Update parameter list with new params
    if(length(new_params) != 0){

        for(name in names(new_params)){
            params[[name]] = new_params[[name]]
        }

    }

    # Update if any parameters changed
    params$BIRTH_RATE = log(params$LITTER_SIZE) # Gives birth rate. Little r in ricker function
    params$BETA = params$BIRTH_RATE / params$K_CAPACITY # TODO: CHECK WHY I DID THIS?

    return(params)

}


get_init_arrays = function(prms){
    # Initial arrays for raccoon simulation

    # Set up the raccoon arrays

    raccoon_dead_alive_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$INIT_NUM_RACCOONS))
    initial_age_vector = rep(10, prms$INIT_NUM_RACCOONS)
    age_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$INIT_NUM_RACCOONS))
    age_array[1, ] = initial_age_vector
    human_array = assign_human_contacts(prms$INIT_NUM_RACCOONS)
    human_risk_through_time = list()
    human_risk_through_time[[1]] = human_array


    new_babies_vect = array(NA, dim=prms$TIME_STEPS + 1)
    new_babies_vect[1] = 0

    # Initialize all arrays
    # Seed a raccoon with some worms
    # raccoon_worm_array[1, ] = 10 # Seeding worms
    raccoon_dead_alive_array[1, ] = 1 # All raccoons are alive

    # Set up worm arrays.  This array hold the cohort of worms for each raccoon
    # so that we can track the age-dependent death of the worms in the raccoons
    infra_mouse_worm_array = lapply(1:prms$INIT_NUM_RACCOONS,
                    function(x) array(NA, dim=c(prms$TIME_STEPS + 1, prms$TIME_STEPS + 1)))

    infra_nonmouse_worm_array = lapply(1:prms$INIT_NUM_RACCOONS,
                    function(x) array(NA, dim=c(prms$TIME_STEPS + 1, prms$TIME_STEPS + 1)))


    # # Initialize worm arrays
    # for(i in 1:length(infra_worm_array)){
    #     infra_worm_array[[i]][1, 1] = raccoon_worm_array[1, i]

    # }

    # Initialize non-rodent worm array
    for(i in 1:length(infra_nonmouse_worm_array)){
        infra_nonmouse_worm_array[[i]][1, 1] = prms$INIT_WORMS

    }

    # Initialize rodent worm array
    for(i in 1:length(infra_mouse_worm_array)){
        infra_mouse_worm_array[[i]][1, 1] = 0

    }

    all_worms_infra_array = list(infra_nonmouse_worm_array, infra_mouse_worm_array)

    # Make total worm arrays
    # infra_worm_array = combine_worm_arrays(all_worms_infra_array, raccoon_dead_alive_array)
    raccoon_worm_array = get_tot_worm_array_from_infra(infra_nonmouse_worm_array, raccoon_dead_alive_array) 


    init_arrays = list(raccoon_dead_alive_array=raccoon_dead_alive_array,
                       initial_age_vector=initial_age_vector,
                       age_array=age_array,
                       human_array=human_array,
                       human_risk_through_time=human_risk_through_time,
                       new_babies_vect=new_babies_vect,
                       infra_mouse_worm_array=infra_mouse_worm_array,
                       infra_nonmouse_worm_array=infra_nonmouse_worm_array,
                       all_worms_infra_array=all_worms_infra_array,
                       raccoon_worm_array=raccoon_worm_array)

    return(init_arrays)
}