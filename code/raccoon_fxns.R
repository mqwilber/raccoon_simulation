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

    prob_death = 1 - exp(beta + log(load + 1)*alpha) / (1 + exp(beta + log(load + 1)*alpha))

    return(prob_death)
}

baby_death_fxn = function(age, baby_death){
    # Gives the probability of dying in the first month of life
    #
    # Parameters
    # ----------
    # age : int, age of raccoon
    # baby_death : float, probability of dying at age 0
    #
    # Returns
    # --------
    # : float, probability of death 

    if(age == 0){
        return(baby_death)
    } else{
        return(0) # Don't die due to baby death
    }
        
}

senescence_fxn = function(age, old_death){
    # Probability of dying of old age
    # 
    # Parameters
    # ----------
    # age : int, age of raccoon
    # old_death: float, dictates how quickly death prob. increases with age
    #
    # Returns
    # -------
    # : float, probability of dying of old age

    prob_sen = min(c(1, old_death * age^2))
    return(prob_sen)
} 

kill_my_raccoon = function(worm_load, age, overlap, death_thresh, patho,
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
    # ... : parameters defined in `get_simulation_parameters`
    # cull: bool, TRUE = cull raccoon, FALSE = don't cull
    #
    # Returns
    # -------
    # : int, 0/FALSE (dead) or 1/TRUE (alive)

    tr = runif(1)

    # Probability of not dying from worms * not dying from baby death * 
    # not dying from old age * not dying from random death
    surv_prob = (1 - death_probability(worm_load, death_thresh, patho)) *
                        (1 - baby_death_fxn(age, baby_death)) * 
                            (1 - random_death_prob) * 
                            (1 - senescence_fxn(age, old_death))

    alive_now = (tr > (1 - surv_prob)) * (!cull)

    return(alive_now)

}

get_trap_indices = function(management_params, raccoon_dead_alive_vect,
                                age_vect, human_overlap_through_time_vect){
    # Obtain indices of raccoons to be trapped in a time step
    # 
    # Possible strategies are:
    # 1. `random`: Keep all raccoons that are trapped
    #   - Additional parameters: `quota`: Number of traps set out per time step
    # 2. `age`: Only keep trapped that are < a certain age i.e. less than a year old
    #   - Additional parameters: `quota`: Number traps per time step
    #   -                        `age`: age in months below which you keep the 
    #                                   raccoons you trap
    #                            `success`: Probability of trapping success
    # 3. `human`: Only trap in regions with a certain human overlap. 
    #               individuals with a human overlap greater than 
    #              `overlap_threshold`
    #   - Additional parameters: `quota`: Number of traps to set out
    #                            `overlap_threshold`: Between 0 and 1,
    #                             The region in which you are trapping
    #
    # Parameters
    # ----------
    # cull_params : list or NULL, cull params, see cull_strategy
    # raccoon_dead_alive_vect: vector, vector specifying whether raccoons 
    #      are dead or alive at time t - 1 
    # age_vect : vector, specifying the ages of alive raccoons at year time
    #                   t - 1
    # human_overlap_through_time_vect: vector, specifies the human overlap value
    #                   of alive raccoons at time t - 1.
    # Returns
    # --------
    # : indices of raccoons to cull
    #

    if(is.null(management_params)){
        return(integer(0)) # return an empty vector
    }

    if(management_params$strategy == "random" | management_params$strategy == "age"){

        # Only choosing currently alive raccoons to trap.
        pop_inds = which((raccoon_dead_alive_vect != 0) & 
                                        (!is.na(raccoon_dead_alive_vect)))

        # 20 km^2 is the size of the whole world. And we adjusted K_CAPACITY to
        # give realistic densities in this world given variation in human overlap
        trapping_area = 20 

    } else if(management_params$strategy == "human"){

        if(is.null(management_params$overlap_threshold)){
            stop("Provide overlap_threshold")
        }

        # NOTE: Using a continuous overlap with discrete zones. Probably should
        # chose overlap_threshold at 0.0, 0.1, 0.2 ... 0.9.
        pop_inds = which(human_overlap_through_time_vect > management_params$overlap_threshold)

        # 20 km^2 is the size of the whole world. And we adjusted K_CAPACITY to
        # give realistic densities in this world given variation in human overlap
        # Only looking at a portion of the 20 km^2 world so scale accordingly
        trapping_area = (1 - management_params$overlap_threshold) * 20

    } else{
        stop(paste(management_params$strategy, "is not a recognized strategy. Try random, age, or human"))
    }


    # Check if there are any raccoons to cull.
    if(length(pop_inds) > 1){

        rac_density = length(pop_inds) / trapping_area

        # Functional response for relating trap probability to raccoon density (Type I)
        # Number is from a regression model fit to capture data Sara collected.
        # Data in trapping.csv and script saved as fit_trapping_function.R
        trap_prob = 1 - exp(-(rac_density * 0.0030545))

        # Traps only catch a random number of raccoons.
        total_caught = rbinom(1, management_params$quota, trap_prob)

        # Cull either quota or how ever many individuals are there < quota. 
        trap_inds = sample(pop_inds, min(c(length(pop_inds), total_caught)))

        # Throw out raccoons that do not conform to age requirement
        if(management_params$strategy == "age"){

            trap_inds = trap_inds[age_vect[trap_inds] < management_params$age]
        }

    } else{
        trap_inds = pop_inds # Empty vector
    }

    return(trap_inds)

}

get_bait_indices = function(worm_control_params, raccoon_dead_alive_vect,
                            human_overlap_through_time_vect) {
    # Obtain indices of raccoons to be baited in a time step
    # 
    # Possible strategies are:
    # 1. `random`: Bait everybody equally
    #   - Additional parameters: `quota`: Number of raccoons to bait per time step
    # 2. `human`: Only bait individuals with a human overlap greater than 
    #              `overlap_threshold`
    #   - Additional parameters: `quota`: absolute number of bait to distribute
    #                            `overlap_threshold`: Between 0 and 1,
    #                             Where the bait is distributed.
    #
    # Parameters
    # ----------
    # worm_control_params : list or NULL
    # raccoon_dead_alive_vect: vector, vector specifying whether raccoons 
    #      are dead or alive at time t - 1 
    # human_overlap_through_time_vect: vector, specifies the human overlap value
    #                   of alive raccoons at time t - 1.
    # Returns
    # --------
    # : indices of raccoons to bait
    #


    if(is.null(worm_control_params)){
        return(integer(0))
    }

    if(worm_control_params$strategy == "random"){

        pop_inds = which((raccoon_dead_alive_vect != 0) & 
                                        (!is.na(raccoon_dead_alive_vect)))

    } else if(worm_control_params$strategy == "human"){

        pop_inds = which(human_overlap_through_time_vect > 
                                    worm_control_params$overlap_threshold)

    } else {
        stop(paste(worm_control_params$strategy, 
                    "is not a recognized strategy. Try random or human"))
    }


    # Check if there are any raccoons to bait
    if(length(pop_inds) > 1){

        # 0.4 is default probability of initial bait not disappearing
        # This 0.4 is based on Smyser 2015
        bait_remaining = rbinom(1, worm_control_params$quota, 0.4)
        worm_control_inds = unique(sample(pop_inds, bait_remaining, replace=TRUE))

    } else {

        worm_control_inds = pop_inds

    }

    return(worm_control_inds)

}

kill_raccoon_worms = function(previous_cohorts, death_thresh, death_slope, 
                got_bait=FALSE){
    # Function to kill worms in raccoon based on worm age
    # Assuming a logistic function of survival prob of worms with cohort
    # age and then killing the worms based on a draw from a binomial
    #
    # Parameters
    # ----------
    # previous_cohorts : vector, worms in raccoon at time - 1 time points
    # death_thresh, death_slope : see `get_simulation_parameters`
    # got_bait : bool, TRUE or FALSE.  If raccoon gets bait worms are dead
    #
    # Return
    # ------
    # : vector, updated cohorts based on dying worms

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

give_birth = function(age_now, time, tot_racs, repro_able,
                        month_at_repro,
                        first_repro_age, litter_size, beta){
    # Decide how many babies are produced by a raccoon. Check if the month
    # is right and if the raccoon is old enough.
    #
    # Parameters
    # ----------
    # age_now : int, age of raccoon in months
    # time : int, current time
    # tot_racs : int, total number of raccoons in population
    # repro_able : int, 0 if rac can't reproduce and 1 if it can reproduce
    # month_at_repro, first_repro_age, litter_size, beta: see `get_simulation_parameters`
    #
    # Returns
    # -------
    # : int, how many babies are produced


    repro = 0 # No reproduction if below conditions aren't met

    if(time %% month_at_repro == 0){ # If the month is right

        if(age_now >= first_repro_age){ # If the age is right

            if(repro_able == 1){ # If the raccoon has never had birth control

                repro_prob = exp(-(beta*tot_racs)) # Ricker function from Encyclopedia of Theoretical Ecology pg 634
                repro = rbinom(1, litter_size, repro_prob)
            }
        }

    }

    return(repro)

}


pick_up_eggs = function(emean, ek, infect, resist, eggs_environ, load,
                            egg_contact_param){
    # Function for juvenile raccoons to pick up eggs
    #
    # Parameters
    # -----------
    # emean, ek, infect, resist, egg_contact_param : see `get_simulation_parameters`
    # load: int, worm load at time t - 1 for a given raccoon
    # eggs_environ: float, zone-dependent eggs in environment
    #
    # Returns
    # -------
    # : int, number of new eggs/worms acquired


    # Exponential decline of infectivity.
    infect_red = infect * exp(-resist * load)

    # Encounter and get eggs on face and get infected with eggs
    eprob = get_eprob(eggs_environ, egg_contact_param)

    new_eggs = rbinom(1, 1, eprob) * rbinom(1, rnbinom(1, size=ek, mu=emean), infect_red)
    return(new_eggs)

}

pick_up_rodents = function(mouse_worm_mean,  mouse_worm_k,
                        rodent_encounter_prob, larval_worm_infectivity){
    # Function to pick up rodents from rodent pool.  Rodent 
    # encounters worm, picks up worms from a negative binomial, and 
    # then worms establish with some prob
    #
    # Parameters
    # ----------
    # mouse_worm_mean, larval_worm_infectivity
    #       egg_contact_param, rodent_encounter_prob : see `get_simulation_parameters`
    # mouse_worm_k : Aggregation of worms in mice.
    #
    # Returns
    # -------
    # : number of larval worms acquired

    # Update rodent mean and variance

    if(mouse_worm_k <= 0 || is.na(mouse_worm_k)){
        # When mean is small enough that the moment estimator for k is < 0, set
        # larval worms to 0.
        larval_worms = 0
    } else{ 
        larval_worms = rnbinom(1, mu=mouse_worm_mean, size=mouse_worm_k) # Larval worms per mouse
    }

    new_worms = rbinom(1, 1, rodent_encounter_prob) * # Encounter with mice
                rbinom(1, larval_worms, larval_worm_infectivity) # Worms establishing

    return(new_worms)

}

update_rodent_mean_var = function(eggs_environ, egg_contact_param, 
                                mouse_worm_mean){
    # Function for updating rodent worm mean and rodent worm aggregation
    # based on environmental egg load. Using a base empirical mean rodent
    # worm load. 
    #
    # Parameters
    # ----------
    # eggs_environ: float, zone-dependent eggs in environment
    # egg_contact_param, mouse_worm_mean : see `get_simulation_parameters`
    #
    # Return
    # ------
    # : list, updated mean and k based on environmental egg load in a zone
    #
    # Note: If new_mean is 0, the variance is undefined, but rnbinom will handle
    # this correctly and return 0 if mean = 0.

    eprob = get_eprob(eggs_environ, egg_contact_param)

    new_mean = eprob * mouse_worm_mean

    # Convert to NBD k
    new_k = update_rodent_k(new_mean)

    return(list(mean=new_mean, k=new_k))

}

update_rodent_k = function(rodent_mean){
    # Use Taylors Law to get new aggregation k given a mean.

    # From Shaw and Dobson 1995 TPL
    logvar = 1.551*log10(rodent_mean) + 1.098

    # Convert to NBD k
    new_k = (rodent_mean)^2 / (10^logvar - rodent_mean)

    return(new_k)

}

get_eprob = function(eggs_environ, egg_contact_param){
    # Calculating our encounter probability from previous prevalences
    #    
    # Parameters
    # ---------- 
    # eggs_environ: float, zone-dependent eggs in environment
    # egg_contact_param: see `get_simulation_parameters`
    #
    # Return
    # ------
    # : float, encounter probability

    log_eggs_environ = log10(eggs_environ + 1)
    eprob = 1 - exp(-egg_contact_param * log_eggs_environ)
    return(eprob)

}

assign_human_contacts = function(num_racs){
    # Assign probabilities of human contacts
    return(runif(num_racs))
}

assign_zone_dispersal = function(k_capacity, zone_density){
    # Assign raccoon dispersal based on carrying capacity and zone density
    # Specifically, take the log ratio between k_capacity and current zone
    # density for all zones.  If all ratios are equal, raccoon has an equal
    # chance of dispersing to all zones.  If they are unequal, assign weights
    # based on scaled ratios and the inverse normal CDF. 
    #
    # Parameters
    # ----------
    # k_capacity : vector, carrying capacity constants that are proportional
    #               to carrying capacity for a given zone
    # zone_density : vector, same length and k_capacity. Gives the raccoon
    #                densities in a zone
    #
    # Returns
    # -------
    # : float, a predicted overlap threshold for a raccoon weighted by 
    #         the space availability in all zones.

    zones = length(k_capacity)

    # Add 0.1 to correct if no raccoons are in the zone
    ratios = log(k_capacity / (zone_density + 0.1))

    # If all ratios are the same, equal probability of dispersal
    if(length(unique(ratios)) == 1){

        weights = rep(1 / zones, zones)

    } else{

        ratios = scale(ratios) # Scale to make ratios comparable
        probs = pnorm(ratios) # Inverse CDF of Normal
        weights = probs / sum(probs) # Standardize to sum to one.

    }

    # Assign zone dispersal
    dispersal_zone = sample(1:zones, 1, prob=weights)
    human_overlap = runif(1, min=dispersal_zone/zones - 1/zones, max=dispersal_zone/zones)

    #return(list(weights, ratios))
    return(human_overlap)
}

get_raccoon_zone = function(overlap, zones){
    # Converts continuous raccoon overlap into discrete zone label 
    #
    # Parameters
    # ----------
    # overlap : vector or float, overlap of raccoon(s)
    # zones : int, number of zones
    #
    # Returns
    # -------
    # : vector or int, the zone corresponding to continuous parameter overlap

    # Group raccoons by zones
    breaks = seq(0, 1, len=zones + 1)
    zone_labels = .bincode(overlap, breaks)

    return(zone_labels)

}

get_zone_density = function(overlap, zones){
    # Gets the population density in each zone
    #
    # Parameters
    # ----------
    # overlap : vector, overlap of raccoon(s)
    # zones : int, number of zones
    # 
    # Returns
    # -------
    # : vector, of length zones with density of raccoons in each zone

    zone_labels = get_raccoon_zone(overlap, zones)
    racs_per_zone = table(zone_labels) # This drops NAs

    # If zone is missing...append to include 0 zones
    if(length(racs_per_zone) != zones){

        missing_labels = setdiff(1:zones, as.integer(names(racs_per_zone)))
        missing_vals = rep(0, length(missing_labels))
        names(missing_vals) = missing_labels
        racs_per_zone = c(racs_per_zone, missing_vals) # Could be slow
    }

    # Return ordered by zone
    return(racs_per_zone[order(as.numeric(names(racs_per_zone)))])
}

assign_egg_production = function(raccoon_worm_vect, human_vect, zones){
    # From the raccoon worm array vector at time t and the human array
    # vector at time t, assign the egg production per zone
    #
    # Parameters
    # ----------
    # raccoon_worm_vect : vector, worms in raccoons at time t
    # human_vect : vector, overlap values for raccoons
    # zones : int,  number of zones
    #
    # Returns
    # -------
    # : vector, egg_production in all zones at time t

    # Group raccoons by zones
    zone_labels = get_raccoon_zone(human_vect, zones)
    # breaks = seq(0, 1, len=zones + 1)
    # zone_labels = .bincode(human_vect, breaks)

    # Calc. number infected in each zone, dropping NAs
    dt = data.table(zone_num=zone_labels, pa=(raccoon_worm_vect > 0))
    grouped_data = dt[!is.na(zone_num), list(pa_sums=sum(pa)), by=c("zone_num")]

    missing_labels = setdiff(1:zones, grouped_data$zone_num)

    # Only do this is we are missing some labels
    if(length(missing_labels) != 0){

        dt_missing = data.table(zone_num=missing_labels, pa_sums=0)
        grouped_data = rbind(grouped_data, dt_missing)

    }

    grouped_data = grouped_data[order(zone_num), ]

    return(grouped_data$pa_sums)

}

get_cumulative_egg_load = function(time, eggproduction_array, egg_decay){
    # Getting time-dependent cumulative egg production per zone
    # 
    # Parameters
    # ----------
    # time : int, time
    # eggproduction_array : array, the num infected raccoons in each zone by time
    #       dim = (time, zones)
    # egg_decay: float, parameter specifying egg decay rate
    #
    # Returns
    # -------
    # : vector
    #       Vector with the remaining cumulative egg load across zones at time
    #       Has length == zones.


    weights = exp(-egg_decay * time:1)

    if(time == 1){
        eggs_remaining_vector = eggproduction_array[1:time, ] * weights
    } else{
        eggs_remaining_vector = apply(eggproduction_array[1:time, ], 2, 
                                                function(x) sum(x*weights))
    }
    return(eggs_remaining_vector)
 
}

update_arrays = function(time, new_babies, new_babies_vect,
                                           initial_age_vector,
                                           raccoon_dead_alive_array,
                                           raccoon_worm_array,
                                           age_array, infra_worm_array,
                                           human_vect,
                                           repro_able_vect,
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
    human_vect = c(human_vect, rep(human_vect[ind],
                                    babies_at_this_time_vect[ind]))

    repro_able_vect = c(repro_able_vect, rep(1, new_babies))

    raccoon_dead_alive_array = cbind(raccoon_dead_alive_array,
                                                new_alive_babies)

    raccoon_worm_array = cbind(raccoon_worm_array, new_worm_babies)

    # Extend the infra_worm_array
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
                human_vect=human_vect,
                repro_able_vect=repro_able_vect))

}

get_human_risk_metric = function(eggproduction_array, egg_decay){
    # Calculating the human risk metric through time
    # This is using a weighted measure of human risk based on the past
    # worm loads in the population
    #
    # Parameters
    # -----------
    # eggproduction_array: array of zone based past egg production
    # egg_decay: parameter determining egg decay rate in the environment

    # Number of raccoons at end of simulation

    full_time = nrow(eggproduction_array)
    zones = ncol(eggproduction_array)

    eggs_remaining = array(NA, dim=c(full_time, zones))

    for(time in 1:full_time){
        eggs_remaining[time, ] = get_cumulative_egg_load(time, eggproduction_array, 
                                        egg_decay)
    }

    # Calculate human risk through time
    human_overlap = seq(0, 1, len=zones + 1)[2:(zones + 1)]
    human_risk = apply(eggs_remaining, 1, function(x) sum(x*human_overlap))

    return(human_risk)
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

get_mean_intensity = function(raccoon_worm_array){
    # Get mean intensity (not including 0s) for worm array
    # Input : the raccoon worm array
    # Returns : prevalence for all time points 

    intensity_fxn = function(x) {
        return(sum(x[x > 0], na.rm=T) / sum(!is.na(x[x > 0])))
    }
    return(apply(raccoon_worm_array, 1, intensity_fxn))

}


## Functions to summarize simulation output ##

plot_pop_traj = function(raccoon_dead_alive_array){
    # plot population trajectories
    pop_trun = rowSums(raccoon_dead_alive_array, na.rm=T)
    dat = data.frame(list(population=pop_trun, time=1:length(pop_trun)))
    ggplot(dat, aes(x=time, y=population)) + geom_point() + geom_line() +
        ylim(c(0, max(dat$population)))
}

plot_age_hist = function(age_array, range, cut=12){

    # Plot the raccoon age distribution at any time point with violin plots

    # age_array: age array from simulation
    # range: indices of the time points at which to make histograms of age
    # distributions. i.e. 30:45

    trun_age = data.frame(t(age_array[range, ]))
    names(trun_age) = range
    stacked_data = stack(trun_age)
    stacked_data_trun = stacked_data[!is.na(stacked_data$values), ]
    ggplot(stacked_data_trun, aes(as.factor(ind), values)) + geom_violin() + 
            theme_bw()

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
    tplot = ggplot(means, aes(x=age, y=mean_worms)) + geom_point() + geom_line()
    return(tplot)
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

plot_prev = function(raccoon_worm_array){
    # Plot whole world prevalence through time

    prev_time = rowSums(raccoon_worm_array > 0, na.rm=T) / 
                            rowSums(raccoon_worm_array > -1, na.rm=T)
    plot_dat = data.frame(time=1:length(prev_time), prev=prev_time)
    prev_plot = ggplot(plot_dat, aes(x=time, y=prev)) + geom_line()
    return(prev_plot)
}

age_intensity_plot = function(raccoon_worm_array, age_array, range){
    # Plot full age intensity profile
    # range specifies which time range to calculate the age-intensity for

    flat_age = as.vector(age_array[range, ])
    flat_worms = as.vector(raccoon_worm_array[range, ])
    full_dat = as.data.table(data.frame(list(age=flat_age, worms=flat_worms)))
    ggplot(full_dat, aes(x=age, y=worms)) + geom_point() + geom_smooth()

}


plot_worm_traj = function(raccoon_worm_array){
    # Plot trajectories of individual raccoons

    sdata = stack(data.frame(raccoon_worm_array))
    sdata$time = rep(1:dim(raccoon_worm_array)[1], dim(raccoon_worm_array)[2])
    tplot = ggplot(sdata, aes(x=time, y=values, color=ind)) + 
                    geom_line() + geom_point() +
                    theme(legend.position="none")
    return(tplot)
}

plot_worm_traj_total = function(raccoon_worm_array){
    # Plot trajectories of individual raccoons

    all_worms = rowSums(raccoon_worm_array, na.rm=TRUE)
    sdata = data.frame(time=1:length(all_worms), worms=all_worms)
    tplot = ggplot(sdata, aes(x=time, y=worms)) + 
                    geom_line() + geom_point()
    return(tplot)
}




## Get parameters and arrays ##

get_simulation_parameters = function(...){
    # Get parameters for the simulation
    #
    # You can pass keyword arguments to this function to reassign the 
    # default parameters.  The keyword parameters must match the parameters
    # defined below
    #
    # Returns
    # -------
    # list : a list of all the simulation parameters

    ## Raccoon parameters
    INIT_NUM_RACCOONS = 500
    DEATH_THRESHOLD = 22 #22 # beta in death_probability fxn
    PATHOGENICITY = -4.2 # alpha in death_probability fxn

    # Probability of dying as baby. 
    BABY_DEATH = 1 - (.52^(1/7))^4 # From data at age 0 in fit_param.R
    RANDOM_DEATH_PROB = 0.01 # Lower bound to death prob
    OLD_DEATH = (1 / (20 * 12)^2) # Above 20 years old the raccoon dies
    AGE_EGG_RESISTANCE = 4 # Age above which raccoons no longer pick up eggs
    RODENT_ENCOUNTER_PROB = 0.5 # Monthly probability of encountering a rodent
    INIT_WORMS = 10 # Initial number of non-rodent worms in raccoons

    ## Rodent parameters
    MOUSE_WORM_MEAN = 3.49 # Abundance of worms in peromyscus estimated from Sara's data
    LARVAL_WORM_INFECTIVITY = 0.25 # Probability of larval worm establishing


    ## Age parameters
    FIRST_REPRO_AGE = 10 # months
    LITTER_SIZE = 3 # Maximum litter size
    MONTH_AT_REPRO = 12
    DISPERSAL_AGE = 12 # months

    # How many discrete zones in our egg zone array
    ZONES = 10

    ## Ricker function for density-dependent recruitment of new babies
    K_CAPACITY = 200*(1:ZONES / sum(1:ZONES)) # rep(200 / ZONES, ZONES) # "Carrying" capacity for raccoons.

    BIRTH_RATE = log(LITTER_SIZE) # Gives birth rate. Little r in ricker function
    BETA = BIRTH_RATE / K_CAPACITY # From Encyclopedia of Theoretical Ecology, pg. 634.

    # PARASITE PARAMETERS
    ENCOUNTER_MEAN = 500 # Mean number of encountered eggs
    ENCOUNTER_K = 1 # Aggregation parameter of the NBD...

    # A set of parameters that determines the
    # probability of encounter given a weighted history of prevalence. The first
    # parameter dictates at what level of the log metric eprob is close to 1 and
    # the second parameter dictates how quickly eprob goes to 1.
    # TODO: WE WILL NEED TO FIDDLE WITH THESE
    EGG_CONTACT = 3 # FIDDLE WITH THIS, probably between 2-4.

    INFECTIVITY = 0.02 # Probability of infectivity
    RESISTANCE = 0.03 # How quickly a raccoon gains resistance based on previous load
    EGG_DECAY = 0.3 # Rate of egg decay such that 2%-3% chance of survival after year. 

    # See fit_param.R for how we got these values
    WORM_SURV_TRESH = 4.7104 #/ 2 # Threshold parameter of worm survival probability
    WORM_SURV_SLOPE = -0.9446 #/ 2 # Slope of worm surv probability

    # Time parameters: Each time step is a month
    TIME_STEPS = 100


    params=list(INIT_NUM_RACCOONS=INIT_NUM_RACCOONS,
                DEATH_THRESHOLD=DEATH_THRESHOLD,
                PATHOGENICITY=PATHOGENICITY,
                BABY_DEATH=BABY_DEATH,
                RANDOM_DEATH_PROB=RANDOM_DEATH_PROB,
                OLD_DEATH=OLD_DEATH,
                AGE_EGG_RESISTANCE=AGE_EGG_RESISTANCE,
                RODENT_ENCOUNTER_PROB=RODENT_ENCOUNTER_PROB,
                INIT_WORMS=INIT_WORMS,
                MOUSE_WORM_MEAN=MOUSE_WORM_MEAN,
                LARVAL_WORM_INFECTIVITY=LARVAL_WORM_INFECTIVITY,
                FIRST_REPRO_AGE=FIRST_REPRO_AGE,
                LITTER_SIZE=LITTER_SIZE,
                MONTH_AT_REPRO=MONTH_AT_REPRO,
                DISPERSAL_AGE=DISPERSAL_AGE,
                K_CAPACITY=K_CAPACITY,
                BIRTH_RATE=BIRTH_RATE,
                BETA=BETA,
                ENCOUNTER_MEAN=ENCOUNTER_MEAN,
                ENCOUNTER_K=ENCOUNTER_K,
                EGG_CONTACT=EGG_CONTACT,
                INFECTIVITY=INFECTIVITY,
                RESISTANCE=RESISTANCE,
                EGG_DECAY=EGG_DECAY,
                WORM_SURV_TRESH=WORM_SURV_TRESH,
                WORM_SURV_SLOPE=WORM_SURV_SLOPE,
                TIME_STEPS=TIME_STEPS,
                ZONES=ZONES)

    new_params = list(...)

    # Update parameter list with new params
    if(length(new_params) != 0){

        for(name in names(new_params)){
            params[[name]] = new_params[[name]]
        }

    }

    # Update if any parameters changed
    params$BIRTH_RATE = log(params$LITTER_SIZE) # Gives birth rate. Little r in ricker function
    params$BETA = params$BIRTH_RATE / params$K_CAPACITY # From Encyclopedia of Theoretical Ecology, pg. 634.

    return(params)

}


get_init_arrays = function(prms){
    # Initial arrays for raccoon simulation
    #
    # Parameters
    # ----------
    # prms : list, the list returned from `get_simulation_parameters`
    #
    # Returns
    # -------
    # : list holding the initialized arrays for the simulation

    # Set up the raccoon arrays

    raccoon_worm_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$INIT_NUM_RACCOONS))
    raccoon_dead_alive_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$INIT_NUM_RACCOONS))
    initial_age_vector = rep(24, prms$INIT_NUM_RACCOONS)
    age_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$INIT_NUM_RACCOONS))
    age_array[1, ] = initial_age_vector

    # Assign raccoons to zones based on carrying capacity of zones
    human_vect = sapply(1:prms$INIT_NUM_RACCOONS, 
            function(x) assign_zone_dispersal(prms$K_CAPACITY, rep(0, prms$ZONES)))

    repro_able_vect = rep(1, prms$INIT_NUM_RACCOONS)
    human_overlap_through_time = list()
    human_overlap_through_time[[1]] = human_vect
    eggproduction_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$ZONES))
    humanrisk_array = array(NA, dim=prms$TIME_STEPS + 1)

    # Rodent mean array
    rodent_mean_array = array(NA, dim=c(prms$TIME_STEPS + 1, prms$ZONES))
    rodent_mean_array[1, ] = prms$MOUSE_WORM_MEAN

    new_babies_vect = array(NA, dim=prms$TIME_STEPS + 1)
    new_babies_vect[1] = 0

    # Initialize all arrays
    # Seed a raccoon with some worms
    # raccoon_worm_array[1, ] = 10 # Seeding worms
    raccoon_dead_alive_array[1, ] = 1 # All raccoons are alive

    # Set up worm arrays.  
    infra_worm_array = lapply(1:prms$INIT_NUM_RACCOONS,
                    function(x) array(NA, dim=c(prms$TIME_STEPS + 1, prms$TIME_STEPS + 1)))

    # Initialize non-rodent worm array
    for(i in 1:length(infra_worm_array)){
        infra_worm_array[[i]][1, 1] = prms$INIT_WORMS

    }

    # Make total worm arrays
    raccoon_worm_array[1, ] =  sapply(1:length(infra_worm_array), 
                                        function(rac) sum(infra_worm_array[[rac]][1, ], na.rm=T))

    # Assign the initial egg production
    eggproduction_array[1, ] = assign_egg_production(raccoon_worm_array[1, ], 
                                        human_vect, prms$ZONES)

    humanrisk_array[1] = get_human_risk_metric(
                        matrix(eggproduction_array[1, ], nrow=1, ncol=prms$ZONES), 
                        prms$EGG_DECAY)

    init_arrays = list(raccoon_dead_alive_array=raccoon_dead_alive_array,
                       initial_age_vector=initial_age_vector,
                       age_array=age_array,
                       human_vect=human_vect,
                       human_overlap_through_time=human_overlap_through_time,
                       new_babies_vect=new_babies_vect,
                       infra_worm_array=infra_worm_array,
                       raccoon_worm_array=raccoon_worm_array,
                       eggproduction_array=eggproduction_array,
                       humanrisk_array=humanrisk_array,
                       repro_able_vect=repro_able_vect,
                       rodent_mean_array=rodent_mean_array)

    return(init_arrays)
}


full_simulation = function(prms, init_arrays, cull_params=NULL, 
                           birth_control_params=NULL, worm_control_params=NULL,
                           latrine_cleanup_params=NULL,
                           management_time=50,  print_it=FALSE){
    # Runs the full individual based simulation
    #
    # Parameters
    # ----------
    # prms : list, simulation parameters
    # init_arrays : initial arrays for holding results
    # print_it : bool, if TRUE prints progress, otherwise False.
    # cull_params : list or NULL, parameters for culling strategies. 
    #    i.e. list(strategy="age", quota=5, overlap_threshold=0.5)
    #    Default is NULL.
    # birth_control_params : list or NULL, parameters for birth control strategies. 
    #    i.e. list(strategy="human", quota=5, overlap_threshold=0.5)
    #    Default is NULL.
    # worm_control_params : list or NULL, TODO
    # management_time : int, time in simulation where management begins.
    #                   Default is 50.
    # 


    # Save the init arrays to the current environment
    list2env(init_arrays, envir=environment())

    # Loop through time steps
    for(time in 2:(prms$TIME_STEPS + 1)){

        if(print_it){
            print(paste("Beginning step", time - 1, "of", prms$TIME_STEPS))
        }


        new_babies = 0
        babies_at_this_time_vect = array(NA, dim=dim(raccoon_worm_array)[2])


        # Total raccoons in a give zone at the previous time step
        tot_racs_by_zone = get_zone_density(human_overlap_through_time[[time - 1]], 
                                                prms$ZONES)

        # If time is >= management time begin management
        if(time >= management_time){

            cull_indices = get_trap_indices(cull_params,
                                         raccoon_dead_alive_array[time - 1, ],
                                         age_array[time - 1, ], 
                                         human_overlap_through_time[[time - 1]])

            # Pick individuals for birth control
            birth_control_indices = get_trap_indices(birth_control_params,
                                         raccoon_dead_alive_array[time - 1, ],
                                         age_array[time - 1, ], 
                                         human_overlap_through_time[[time - 1]])

            # If picked for birth control, can't have babies
            repro_able_vect[birth_control_indices] = 0 

            worm_control_indices = get_bait_indices(worm_control_params,
                                                    raccoon_dead_alive_array[time - 1, ],
                                                    human_overlap_through_time[[time - 1]])

            # Latrine clean up
            if(!is.null(latrine_cleanup_params)){

                # Find minimum overlap zone to cleanup
                lower = ceiling(latrine_cleanup_params$overlap_threshold * prms$ZONES)

                    eggproduction_array[1:(time - 1), lower:prms$ZONES] = 
                        eggproduction_array[1:(time - 1), lower:prms$ZONES] * 0#* (1 - latrine_cleanup_params$cleanup_efficiency)
            }


        } else{
            birth_control_indices = integer(0) 
            cull_indices = integer(0)
            worm_control_indices = integer(0) 
        }

        # Calculate eggs remaining in each zone for this time step
        eggs_remaining = get_cumulative_egg_load(time - 1, eggproduction_array, 
                                            prms$EGG_DECAY)

        # Getting the rodent mean
        rodent_mean_array[time, ] = unlist(sapply(eggs_remaining, update_rodent_mean_var, 
                                    prms$EGG_CONTACT, prms$MOUSE_WORM_MEAN)[1, ])

        # Loop through raccoons but only through raccoons that are alive
        alive_inds = which(raccoon_dead_alive_array[time - 1, ] == 1)
        dead_inds = which(raccoon_dead_alive_array[time - 1, ] == 0)

        # Keep the dead raccoons dead
        raccoon_worm_array[time, dead_inds] = NA
        raccoon_dead_alive_array[time, dead_inds] = 0

        for(rac in alive_inds){

            # 1. Raccoon dies (Killing at the beginning of the time interval)

            age_now = initial_age_vector[rac] +
                        sum(raccoon_dead_alive_array[, rac], na.rm=T)

            # Get overlap and raccoon zones
            overlap_now = human_overlap_through_time[[time - 1]][rac]

            # NOTE: We could speed this up by assigning raccoons discrete zones
            zone_now = get_raccoon_zone(overlap_now, prms$ZONES)

            alive_now = kill_my_raccoon(raccoon_worm_array[time - 1, rac],
                                            age_now, overlap_now, 
                                            prms$DEATH_THRESHOLD, prms$PATHOGENICITY,
                                            prms$BABY_DEATH,
                                            prms$RANDOM_DEATH_PROB, prms$OLD_DEATH,
                                            cull=any(rac == cull_indices))

            raccoon_dead_alive_array[time, rac] = alive_now

            if(alive_now){ # Do this if raccoon is alive

                # 1. Give birth if appropriate

                age_array[time, rac] = age_now

                new_babies_now = give_birth(age_now, time, 
                                                    tot_racs_by_zone[zone_now], 
                                                    repro_able_vect[rac],
                                                    prms$MONTH_AT_REPRO,
                                                    prms$FIRST_REPRO_AGE,
                                                    prms$LITTER_SIZE, 
                                                    prms$BETA[zone_now])

                # Add new babies to array to assign human contacts later
                babies_at_this_time_vect[rac] = new_babies_now
                new_babies = new_babies + new_babies_now

                # 2. Kill old worms for each possible source of worms
                previous_cohorts = infra_worm_array[[rac]][time - 1, 1:(time - 1)]
                new_cohort = kill_raccoon_worms(previous_cohorts,
                                                    prms$WORM_SURV_TRESH,
                                                    prms$WORM_SURV_SLOPE,
                                                    got_bait=any(rac == worm_control_indices))

                # Assign previous cohort to current cohort
                infra_worm_array[[rac]][time, 1:(time - 1)] = new_cohort

                # 3. Pick up eggs or pick up rodents depending on age

                if(age_now <= prms$AGE_EGG_RESISTANCE){ # Only pick up worms from eggs

                    worms_acquired = pick_up_eggs(prms$ENCOUNTER_MEAN,
                                            prms$ENCOUNTER_K,
                                            prms$INFECTIVITY, prms$RESISTANCE,
                                            eggs_remaining[zone_now],
                                            raccoon_worm_array[time - 1, rac],
                                            prms$EGG_CONTACT)


                } else{ # Only pick up worms from rodent

                    worms_acquired = pick_up_rodents(rodent_mean_array[time - 1, zone_now],
                                                     update_rodent_k(rodent_mean_array[time - 1, zone_now]),
                                                     prms$RODENT_ENCOUNTER_PROB,
                                                     prms$LARVAL_WORM_INFECTIVITY)

                }

                infra_worm_array[[rac]][time, time] = worms_acquired
                raccoon_worm_array[time, rac] = sum(infra_worm_array[[rac]][time, ], na.rm=T)

                # }

                # 4. Disperse if raccoon is 12
                # TODO: Zone density-based on t or t - 1
                if(age_now == prms$DISPERSAL_AGE){
                    human_vect[rac] = assign_zone_dispersal(prms$K_CAPACITY, 
                                            tot_racs_by_zone)
                }

            } else { # Raccoon dead and its worms die
                raccoon_worm_array[time, rac] = NA
                human_vect[rac] = NA
                repro_able_vect[rac] = NA
            }

        }

        # Update the various vectors if babies were born
        if(new_babies > 0){

            updated_arrays = update_arrays(time, new_babies, new_babies_vect,
                                               initial_age_vector,
                                               raccoon_dead_alive_array,
                                               raccoon_worm_array,
                                               age_array, infra_worm_array,
                                               human_vect,
                                               repro_able_vect,
                                               babies_at_this_time_vect,
                                               prms$TIME_STEPS)

            new_babies_vect = updated_arrays$new_babies_vect
            initial_age_vector = updated_arrays$initial_age_vector
            age_array = updated_arrays$age_array
            raccoon_dead_alive_array = updated_arrays$raccoon_dead_alive_array
            raccoon_worm_array = updated_arrays$raccoon_worm_array
            infra_worm_array = updated_arrays$infra_worm_array
            human_vect = updated_arrays$human_vect
            repro_able_vect = updated_arrays$repro_able_vect

        }

        # Save the new human overlap and human risk arrays
        human_overlap_through_time[[time]] = human_vect
        eggproduction_array[time, ] = assign_egg_production(
                                            raccoon_worm_array[time, ],
                                            human_vect,
                                            prms$ZONES)
        humanrisk_array[time] = get_human_risk_metric(eggproduction_array[1:time, ], 
                                                             prms$EGG_DECAY)[time]

    }

    return(list(new_babies_vect=new_babies_vect, 
                initial_age_vector=initial_age_vector,
                age_array=age_array,
                raccoon_dead_alive_array=raccoon_dead_alive_array,
                raccoon_worm_array=raccoon_worm_array,
                infra_worm_array=infra_worm_array,
                human_vect=human_vect,
                repro_able_vect=repro_able_vect,
                human_overlap_through_time=human_overlap_through_time,
                eggproduction_array=eggproduction_array,
                humanrisk_array=humanrisk_array,
                rodent_mean_array=rodent_mean_array))
} # End full simulation


run_and_extract_results = function(i, quota, management_time, 
                                    cull_params, birth_control_params, 
                                    worm_control_params, latrine_cleanup_params, 
                                    time_steps){
    # Runs simulations and extracts results for a given quota. 
    # Extracts results from last 12 months of the simulation.
    #
    # Parameters
    # ----------
    # i : integer, simulation index 
    # quota : double. Use NULL unless quota means something in the simulation. 
    #           value does not influence analysis
    # management_time : int, the simulation time at which the management starts
    # cull_params , birth_control_params, worm_control_params, latrine_cleanup_params : 
    #               Management parameters.
    # time_steps : int, number of time steps in the simulation

    print(paste("Simulation", i, "for quota", quota))

    params = get_simulation_parameters(TIME_STEPS=time_steps) # Load in simulation parameters
    init_arrays = get_init_arrays(params) # Load in init arrays

    all_res = full_simulation(params, init_arrays,
                                cull_params=cull_params, 
                                birth_control_params=birth_control_params, 
                                worm_control_params=worm_control_params,
                                latrine_cleanup_params=latrine_cleanup_params, 
                                management_time=management_time)

    # Extract min, max rac pop sizes
    decrement = 24 # Use the last two years of the simulation
    pop_traj = rowSums(all_res$raccoon_dead_alive_array, na.rm=T)
    tot_steps = length(pop_traj)
    max_rac_pop = max(pop_traj[(tot_steps - decrement):tot_steps])
    min_rac_pop = min(pop_traj[(tot_steps - decrement):tot_steps])
    mean_rac_pop = mean(pop_traj[(tot_steps - decrement):tot_steps])

    # Extract min, max, worm pop sizes
    worm_pop_traj = rowSums(all_res$raccoon_worm_array, na.rm=T)
    max_worm_pop = max(worm_pop_traj[(tot_steps - decrement):tot_steps])
    min_worm_pop = min(worm_pop_traj[(tot_steps - decrement):tot_steps])
    mean_worm_pop = mean(worm_pop_traj[(tot_steps - decrement):tot_steps])

    # Human risk
    human_risk = all_res$humanrisk_array 
    max_human_risk = max(human_risk[(tot_steps - decrement):tot_steps])
    min_human_risk = min(human_risk[(tot_steps - decrement):tot_steps])
    mean_human_risk = mean(human_risk[(tot_steps - decrement):tot_steps])

    # Prevalence measures
    prev_traj = get_prevalence(all_res$raccoon_worm_array)
    max_prev = max(prev_traj[(tot_steps - decrement):tot_steps])
    min_prev = min(prev_traj[(tot_steps - decrement):tot_steps])
    mean_prev = mean(prev_traj[(tot_steps - decrement):tot_steps])

    # Mean intensity measures
    intensity_traj = get_mean_intensity(all_res$raccoon_worm_array)
    max_intensity = max(intensity_traj[(tot_steps - decrement):tot_steps])
    min_intensity = min(intensity_traj[(tot_steps - decrement):tot_steps])
    mean_intensity = mean(intensity_traj[(tot_steps - decrement):tot_steps])

    return(c(min_rac_pop, mean_rac_pop, max_rac_pop, 
             min_worm_pop, mean_worm_pop, max_worm_pop,
             min_human_risk, mean_human_risk, max_human_risk,
             min_prev, mean_prev, max_prev,
             min_intensity, mean_intensity, max_intensity))

}