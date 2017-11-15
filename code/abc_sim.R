## ABC-SMC algorithm for estimating four unknown parameters of RacoonWormModel.
## 
## Author: Mark Wilber

library(parallel)
library(ggplot2)
source("raccoon_fxns.R") # Contains functions for simulating IBM

# Uniform upper and lower priors
param_priors = list(RODENT_ENCOUNTER_PROB=c(min=0, max=1),
                    ENCOUNTER_K=c(min=0, max=1),
                    EGG_CONTACT=c(min=1e-7, max=1e-4),
                    AGE_SUSCEPTIBILITY=c(min=0.001, max=20),
                    CLEAR_PROB=c(min=0, max=1),
                    AGE_EGG_RESISTANCE=c(min=0, max=10))

compare_to_data = function(all_res, time_steps, stat_set="all", 
                    datasource="../data/formatted/raccoon_age_intensity_full.csv", lag=24){
    # Given a simulation result, extract the simulation data that matches with
    # the empirical data so they can be compared.
    #
    # Parameters
    # ----------
    # all_res : results from the IBM model simulation
    # time_steps : number of time steps the model was run for
    # lag : the number of time steps at the end of the model from which to
    #       extract raccoons to compare to data
    # datasource : File that holds the empirical data/simulated data
    # stat_set : str, "all": use mean, prev and IQR
    #                 "means": use only means
    #                 "prevs" : use only prevs
    #                 "means_iqr" : use means and iqr
    #                 "means_prevs" : use means and prevs
    #
    # Returns
    # ------
    # : means, prevs, iqrs
    #    means is the age-abundance profile of raccoons for the same age bins
    #    as the and prevs is the age-prevalence profile. iqrs is the
    #    interquartile range for the abundance in each age class. Depending
    #    on the stat_set argument, only particular combinations of these may be
    #    returned.

    ages = all_res$age_array[(time_steps + 1 - lag):(time_steps + 1), ]
    dim(ages) = NULL
    worms = all_res$raccoon_worm_array[(time_steps + 1 - lag):(time_steps + 1), ]
    dim(worms) = NULL

    age_dat = data.frame(age=ages, worm=worms)
    age_dat = age_dat[complete.cases(age_dat), ]

    obs_dat = read.csv(datasource)
    means = array(NA, dim=nrow(obs_dat))
    prevs = array(NA, dim=nrow(obs_dat))
    iqrs = array(NA, dim=nrow(obs_dat))

    for(i in 1:nrow(obs_dat)){

        lower = obs_dat$age_lower[i] 
        upper = obs_dat$age_upper[i] 
        trun_age = age_dat[age_dat$age >= lower & 
                                    age_dat$age < upper, ]

        # Sample (can lead to some problems with small popuulations)
        inds = tryCatch({

            inds = sample(1:nrow(trun_age), obs_dat$sample_size[i], replace=FALSE)

          }, error = function(err){

            # For small populations, allow for resampling of the same individual
            inds = sample(1:nrow(trun_age), obs_dat$sample_size[i], replace=TRUE)
            return(inds)

          })
        #inds = sample(1:nrow(trun_age), obs_dat$sample_size[i], replace=FALSE)

        means[i] = mean(trun_age$worm[inds])
        prevs[i] = mean(trun_age$worm[inds] > 0)
        iqrs[i] = diff(quantile(trun_age$worm[inds], c(0.25, 0.75)))

    }

    if(stat_set == "all")
        stats = c(means, prevs, iqrs)
    else if(stat_set == "means")
        stats = c(means)
    else if(stat_set == "prevs")
        stats = c(prevs)
    else if(stat_set == "means_prevs")
        stats = c(means, prevs)

    return(stats)
}

simulate_and_compare = function(i, abc_params, time_steps=100, stat_set="all",
                            datasource="../data/formatted/raccoon_age_intensity_full.csv"){
    # Simulate the raccoon model and compare it to the age intensity/prevalence
    # data.
    #
    # Parameters
    # ----------
    # i : int, and index for simulating and tracking progress
    # abc_params : matrix of N rows and 4 columns where the columns correspond
    #               to the parameters of interest drawn from their priors
    # time_steps : How many time steps to run the model.
    # stat_set : see `compare_to_data` function
    # datasource : path to the observed data
    #
    # Returns
    # -------
    # : list
    #       `summary_stats`: Summary statistics for each simulation
    #       `params`: Parameter values for each simulation

    print(paste("Working on sim", i))
    abc_params_vect = abc_params[[i]]

    params = do.call(get_simulation_parameters, c(list(TIME_STEPS=time_steps), 
                                                    as.list(abc_params_vect)))
    init_arrays = get_init_arrays(params) # Load in init arrays

    all_res = full_simulation(params, init_arrays, 
                                  cull_params=NULL, 
                                  birth_control_params=NULL,
                                  worm_control_params=NULL,
                                  latrine_cleanup_params=NULL, 
                                  management_time=1000,
                                  print_it=FALSE)

    summary_stats = compare_to_data(all_res, time_steps, stat_set=stat_set, 
                                        datasource=datasource)

    return(list(stats=summary_stats, params=abc_params_vect))

}

distances = function(stats, method="euclidean",
                    stat_set="all", 
                    datasource="../data/formatted/raccoon_age_intensity_full.csv"){
    # Get distances between stats and observed vector. Standardize and 
    # convert into pca space to deal with different scales and eliminate
    # correlation.  Distances are calculated as Euclidean distances after 
    # standardization.  Observed data is also put on the same scale as the
    # predicted data  
    #
    # Parameters
    # ----------
    # stats : Matrix of summary statistics (16 columns) that holds
    #           age-intensity/prevalence profiles from the simulations
    # method : Can be any distance method that is taken by the `dist` function in R.
    #          "euclidean" is the L2 norm and "manhattan" is the L1 norm. 
    #
    # Returns
    # -------
    # : The normalized distances between the observed data and simulated data

    # TODO: I don't think I really need to do PCA since all that is doing is
    # rotating the data.
    #
    # TODO: Should weight each vector by sample size such that the data points 
    # that have more samples get more weight.

    sims = nrow(stats)
    non_zero_inds = apply(stats, 2, sd) > 0
    stats = stats[, non_zero_inds] # drop any constant columns as they don't help distinguish between distances
    pca_stats = scale(stats)
    # pca_stats = princomp(scale(stats))$scores

    # Convert observed into PCA space...what are the problems with this approach?
    # Essentially I am asking how similar is is the predicted data to the observed 
    # data when a remove correlations and standardize.
    obs_dat = read.csv(datasource)

    if(stat_set == "all")
        obs_stat = c(obs_dat$abundance, obs_dat$prev / 100, obs_dat$iqr)[non_zero_inds]
    else if(stat_set == "means")
        obs_stat = c(obs_dat$abundance)[non_zero_inds]
    else if(stat_set == "prevs")
        obs_stat = c(obs_dat$prev / 100)[non_zero_inds]
    else if(stat_set == "means_prevs")
        obs_stat = c(obs_dat$abundance, obs_dat$prev / 100)[non_zero_inds]

    scaled_obs = t(matrix((obs_stat - attr(scale(stats), "scaled:center")) / 
                        attr(scale(stats), "scaled:scale")))
    pca_obs = scaled_obs

    # loadings = eigen(cov(scale(stats)))$vectors
    # pca_obs = scaled_obs %*% loadings

    # Distances between observed and predicted data
    dists = as.matrix(dist(rbind(pca_stats, pca_obs), method=method))[sims + 1, 1:(sims)]

    return(dists)
}

resample_and_perturb = function(params, weights, num_samps, perturb_sds, model){
    # Bootstrap resample a set of parameters and perturb them with a uniform
    # kernel.  This is the sampling and perturbation step of the ABC-SMC
    # algorithm
    #
    # Parameters
    # ----------
    # params : matrix of parameter vectors to resample
    # weights : weights on each of the parameter vectors
    # num_samps : Number of bootstrap resamples to perform
    # perturb_sds : The perturbation standard deviations
    # model: Either 1, 2, 3, 4, or 5.  
    #          1: Model with concomitant immunity starting at age 4. 
    #          2: Model with concomitant immunity gradually increasing at age 4.
    #          3: Model with concomitant immunity gradually increasing at 4 months and worm clearance.
    #          4: Model with concomitant immunity gradually increasing a fitted age and worm clearance 
    #          5: Same as model 4, but clearance probability is 0.
    #
    # Returns
    # -------
    # : a matrix of perturbed, resampled parameters

    if(nrow(params) != 0){

        inds = 1:nrow(params)
        new_param_inds = sample(inds, num_samps, replace=TRUE, prob=weights)
        new_params = params[new_param_inds, ]

        # Uniform perturbation
        new_params_perturbed = t(apply(new_params, 1, 
                function(x) x + (perturb_sds)*runif(length(x), min=-1, max=1)))

        # Normal perturbation
        # new_params_perturbed = t(apply(new_params, 1, 
        #         function(x) rnorm(length(x), mean=x, sd=perturb_sds)))

        # Ensure the newly perturbed parameters are not violating the priors. If so
        # sample them again
        new_params_perturbed = check_params(new_params_perturbed, new_params, model)
        return(new_params_perturbed)
    } else{
        return(numeric())
    }

}

check_params = function(params_perturbed, params, model){
    # Makes sure each parameter is within its proper bounds
    # If a pertubation pushed it out, set it to old parameters 
    # (TODO: Change this to simply resample)
    #
    # Parameters
    # ----------
    # params_perturbed : An array of perturbed parameters
    # params : An array of unperturbed parameters
    # model: Either 1, 2, 3, 4, or 5.  
    #          1: Model with concomitant immunity starting at age 4. 
    #          2: Model with concomitant immunity gradually increasing at age 4.
    #          3: Model with concomitant immunity gradually increasing at 4 months and worm clearance.
    #          4: Model with concomitant immunity gradually increasing a fitted age and worm clearance 
    #          5: Same as model 4, but clearance probability is 0.

    dp = dim(params)

    if(model == 1){

        upper = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['max'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['max'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['max']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)
        lower = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['min'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['min'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['min']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)

    } else if(model == 2){

        upper = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['max'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['max'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['max'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['max']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)
        lower = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['min'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['min'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['min'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['min']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)

    } else if(model == 3){

        upper = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['max'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['max'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['max'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['max'],
                             CLEAR_PROB=param_priors[['CLEAR_PROB']]['max']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)
        lower = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['min'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['min'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['min'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['min'],
                             CLEAR_PROB=param_priors[['CLEAR_PROB']]['min']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)

    } else if(model == 4){

        upper = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['max'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['max'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['max'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['max'],
                             CLEAR_PROB=param_priors[['CLEAR_PROB']]['max'],
                             AGE_EGG_RESISTANCE=param_priors[['AGE_EGG_RESISTANCE']]['max']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)
        lower = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['min'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['min'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['min'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['min'],
                             CLEAR_PROB=param_priors[['CLEAR_PROB']]['min'],
                             AGE_EGG_RESISTANCE=param_priors[['AGE_EGG_RESISTANCE']]['min']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)

    } else if(model == 5){

        upper = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['max'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['max'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['max'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['max'],
                             AGE_EGG_RESISTANCE=param_priors[['AGE_EGG_RESISTANCE']]['max']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)
        lower = matrix(rep(c(ENCOUNTER_K=param_priors[['ENCOUNTER_K']]['min'],
                             EGG_CONTACT=param_priors[['EGG_CONTACT']]['min'],
                             RODENT_ENCOUNTER_PROB=param_priors[['RODENT_ENCOUNTER_PROB']]['min'],
                             AGE_SUSCEPTIBILITY=param_priors[['AGE_SUSCEPTIBILITY']]['min'],
                             AGE_EGG_RESISTANCE=param_priors[['AGE_EGG_RESISTANCE']]['min']),
                             dp[1]),
                             nrow=dp[1], ncol=dp[2], byrow=T)

    }

    gt = params_perturbed >= upper
    lt = params_perturbed <= lower 
    params_perturbed[gt] = params[gt]
    params_perturbed[lt] = params[lt]

    return(params_perturbed)
}

get_particles = function(num_particles, model){
    # Sample num_particles parameter vectors from the prior distributions of 
    # each of the parameters.
    # 
    # model: Either 1, 2, 3, 4, or 5.  
    #          1: Model with concomitant immunity starting at age 4. 
    #          2: Model with concomitant immunity gradually increasing at age 4.
    #          3: Model with concomitant immunity gradually increasing at 4 months and worm clearance.
    #          4: Model with concomitant immunity gradually increasing a fitted age and worm clearance 
    #          5: Same as model 4, but clearance probability is 0.
    #
    # Returns
    # -------
    # : Matrix of parameters where each row is a parameter vector

    # Prior distributions on parameters

    # ENCOUNTER_K is one the scale (1 / (1 + k))
    ENCOUNTER_K = runif(num_particles, min=param_priors[['ENCOUNTER_K']]['min'], 
                                       max=param_priors[['ENCOUNTER_K']]['max'])
    EGG_CONTACT = runif(num_particles, min=param_priors[['EGG_CONTACT']]['min'], 
                                       max=param_priors[['EGG_CONTACT']]['max'])
    RODENT_ENCOUNTER_PROB = runif(num_particles, min=param_priors[['RODENT_ENCOUNTER_PROB']]['min'], 
                                                 max=param_priors[['RODENT_ENCOUNTER_PROB']]['max'])
    AGE_SUSCEPTIBILITY = runif(num_particles, min=param_priors[['AGE_SUSCEPTIBILITY']]['min'], 
                                              max=param_priors[['AGE_SUSCEPTIBILITY']]['max'])
    CLEAR_PROB = runif(num_particles, min=param_priors[['CLEAR_PROB']]['min'], 
                                      max=param_priors[['CLEAR_PROB']]['max'])

    AGE_EGG_RESISTANCE = runif(num_particles, min=param_priors[['AGE_EGG_RESISTANCE']]['min'], 
                                      max=param_priors[['AGE_EGG_RESISTANCE']]['max'])

    # Model 1 is default
    params = cbind(ENCOUNTER_K, EGG_CONTACT, RODENT_ENCOUNTER_PROB)

    if(model == 1)
        params = params
    else if(model == 2)
        params = cbind(params, AGE_SUSCEPTIBILITY)
    else if(model == 3)
        params = cbind(params, AGE_SUSCEPTIBILITY, CLEAR_PROB)
    else if(model == 4)
        params = cbind(params, AGE_SUSCEPTIBILITY, CLEAR_PROB, AGE_EGG_RESISTANCE)
    else if(model == 5)
        params = cbind(params, AGE_SUSCEPTIBILITY, AGE_EGG_RESISTANCE)

    return(params)

}


calculate_weights = function(prev_weights, prev_params, 
                    current_params, perturb_sds){
    # Calculate the weights of the current parameters using the perturbation
    # kernel and the weights of the past parameters.  This is the weighting 
    # step in the ABC-SMC algorithm.
    #
    # Parameters
    # ----------
    # prev_weights : vector of previous weights on previous parameters
    # prev_params : Matrix of previous steps parameter vector
    # current_params : The current parameter matrix on which we want to compute
    #                   weights
    # perturb_sds : The previous perturbed standard deviations used to generate
    #                  the current perturbed parameters
    #
    # Returns
    # -------
    # : Normalized weights on current parameters

    # Loop through each current parameter vector
    weights = array(NA, dim=nrow(current_params))

    for(i in 1:nrow(current_params)){

        current = current_params[i, ]
        marginal = array(NA, dim=nrow(prev_params))

        # Loop through each previous parameter vector
        for(j in 1:nrow(prev_params)){

            previous = prev_params[j, ]
            pw = prev_weights[j]

            # Uniform perturbation
            vals = dunif((current - previous) / perturb_sds, min=-1, max=1)

            # Normal perturbation
            #vals = dnorm(current, mean=previous, sd=perturb_sds)
            marginal[j] = pw * prod(vals)
        }

        weights[i] = 1 / sum(marginal) # 1 in numerator due to uniform dist
    }

    return(weights / sum(weights))
}

build_simulated_datasets = function(model_params, num_sets=1, TIME_STEPS=100,
                                datasource="../data/formatted/raccoon_age_intensity_full.csv"){
    # Simulate the raccoon IBM and then build a mock dataset
    #
    # Parameters
    # ----------
    # model_params : list, should include RODENT_ENCOUNTER_PROB, ENCOUNTER_K, 
    #        EGG_CONTACT, ENCOUNTER_MEAN, AGE_SUSCEPTIBILITY, 
    #        CLEAR_PROB, AGE_EGG_RESISTANCE
    # TIME_STEPS: Number of time_steps to run the model. 

    params = do.call(get_simulation_parameters, c(list(TIME_STEPS=TIME_STEPS), 
                                                    model_params))

    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(params, init_arrays, 
                                  cull_params=NULL, 
                                  birth_control_params=NULL,
                                  worm_control_params=NULL,
                                  latrine_cleanup_params=NULL, 
                                  management_time=1000,
                                  print_it = TRUE)

    # Get the simulated data
    obs_dat = read.csv(datasource)
    mock_datasets = list()
    for(i in 1:num_sets){

        pred = compare_to_data(all_res, TIME_STEPS, datasource=datasource)
        abund = pred[1:8]
        prev = pred[9:16]
        iqr = pred[17:24]
        age = obs_dat$age_lower
        sim_num = rep(i, length(age))
        tres = cbind(sim_num, abund, prev, iqr, age)

        mock_datasets[[i]] = tres
    }

    return(mock_datasets)
}

get_sds = function(models, new_params, new_model_sample, old_perturb_sds){
    # Get the standard deviations necessary for pertubring a kernel

    new_perturb_sds = list()
    for(x in models){

        psds = apply(do.call(rbind, new_params[new_model_sample == x]), 2, sd)

        # If these is only one kernel (i.e. SD doesn't exist) use old sds
        if(all(is.na(psds)))
            psds = old_perturb_sds[[x]]

        new_perturb_sds[[x]]  = psds
    }

    return(new_perturb_sds)
}


run_abc = function(steps, num_particles, models,  stat_set="all", method="euclidean",
            datasource="../data/formatted/raccoon_age_intensity_full.csv", 
            percent_rj=0.05, cores=4){
    # Run the abc algorithm for some steps and some particles to fit raccoon
    # model based on some age-intensity profile in datasource.
    # 
    # Parameters
    # ----------
    # steps: int, number of iterations in the sequential monte carlo
    # num_particles: int, number of particles (parameter vectors) to sample
    # stat_set : str, see `compare_to_data` function
    # method : str, method for calculating distance
    # datasource: str, path to the data to which the model is being fit
    # models : Any combination of 1, 2, or 3.  
    #          1: Model with concomitant immunity starting at some age. 
    #          2: Model with concomitant immunity gradually increasing
    #          3: Model with concomitant immunity gradually increasing and worm clearance.
    #          4: Model with concomitant immunity gradually increasing a fitted age and worm clearance 
    #
    # Returns
    # -------
    # : list contains (past_params, weight_array)
    #       past_params is of length steps and each item contains an array
    #       of the sampled particles after than step. Weights contains the 
    #       weights of each parameter set.


    # Sample models with equal weights
    models = models[order(models)]
    current_model_sample = unlist(sample(as.list(models), num_particles, replace=TRUE))
    current_model_sample = current_model_sample[order(current_model_sample)]

    current_params_by_model = lapply(models, function(x) get_particles(sum(current_model_sample == x), x))

    # Unpack parameters into a list
    current_params = do.call(c, lapply(1:length(models), function(x) lapply(1:nrow(current_params_by_model[[x]]), 
                                                    function(i) current_params_by_model[[x]][i, ])))

    perturb_sds = lapply(models, function(x) apply(do.call(rbind, current_params[current_model_sample == x]), 2, sd))

    weight_array = list() # Holds all parameter weights
    past_params = list()  # Holds all parameter sets
    past_models = list()

    for(t in 1:steps){

        print(paste("Iteration", t))

        # Parallelizing over all parameter sets for all models
        saveRDS(current_params, "current_params.rds")
        res = mclapply(1:num_particles, simulate_and_compare, current_params, 
                                datasource=datasource, stat_set=stat_set,
                                mc.cores=cores)

        # Extract the results
        print(res)
        stats = do.call(rbind, lapply(res, function(x) x$stats))

        # Extract the parameter values and put them in a list
        params = lapply(res, function(x) x$params)

        # Standardize and convert the predicted statistics into PCA space
        dists = distances(stats, method=method, stat_set=stat_set, 
                                                datasource=datasource)

        # Extract the alpha = X% smallest distances
        lower = quantile(dists, c(percent_rj))
        lowest_inds = dists < lower

        # Get the new params
        new_params = params[lowest_inds]
        past_params[[t]] = new_params

        # Update and save model indices
        new_model_sample = current_model_sample[lowest_inds]
        past_models[[t]] = new_model_sample

        # Calculate weights
        if(t == 1){
            # Compute models separately as lists and recombine
            weights = do.call(c, lapply(models, function(x) rep(1  / sum(new_model_sample == x), sum(new_model_sample == x))))
        } else{

            weights = do.call(c, lapply(models, function(x) calculate_weights(weight_array[[t - 1]][past_models[[t - 1]] == x], 
                                                  toarray(prev_params[past_models[[t - 1]] == x]), 
                                                  toarray(new_params[new_model_sample == x]), 
                                                  perturb_sds[[x]])))
        }

        weight_array[[t]] = weights

        # Perturb using uniform perturbation
        perturb_sds = get_sds(models, new_params, new_model_sample, perturb_sds)

        # Sample new model indexes from model prior
        current_model_sample = unlist(sample(as.list(models), num_particles, replace=TRUE))
        current_model_sample = current_model_sample[order(current_model_sample)]
        current_params = do.call(c, lapply(models, function(x) tolist(resample_and_perturb(toarray(new_params[new_model_sample == x]), 
                                                 weights[new_model_sample == x], 
                                                 sum(current_model_sample == x), 
                                                 perturb_sds[[x]], x))))

        prev_params = new_params
    }

    return(list(past_params=past_params, weight_array=weight_array, past_models=past_models))
}

toarray = function(x){
    if(length(x) == 0)
        return(numeric())
    else
        return(do.call(rbind, x))
}

tolist = function(x){
    if(length(x) == 0)
        return(list())
    else
        return(lapply(1:nrow(x), function(i) x[i, ]))
}

