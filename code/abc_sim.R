## ABC-SMC algorithm for estimating two/three parameters of RacoonWormModel.
## T

library(parallel)
library(ggplot2)
source("raccoon_fxns.R")
datasource = "../data/raccoon_age_intensity.csv"


compare_to_data = function(all_res, time_steps, lag=24){
    # Given a simulation result, extract the simulation data that matches with
    # the empirical data so they can be compared.
    #
    # Parameters
    # ----------
    # all_res : results from the IBM model simulation
    # time_steps : number of time steps the model was run for
    # lag : the number of time steps at the end of the model from which to
    #       extract raccoons to compare to data
    #
    # Returns
    # ------
    # : means, prevs
    #    means is the age-abundance profile of raccoons for the same age bins
    #    as the and prevs is the age-prevalence profile.

    ages = all_res$age_array[(time_steps + 1 - lag):(time_steps + 1), ]
    dim(ages) = NULL
    worms = all_res$raccoon_worm_array[(time_steps + 1 - lag):(time_steps + 1), ]
    dim(worms) = NULL

    age_dat = data.frame(age=ages, worm=worms)
    age_dat = age_dat[complete.cases(age_dat), ]

    dat = read.csv(datasource)
    means = array(NA, dim=nrow(dat))
    prevs = array(NA, dim=nrow(dat))

    for(i in 1:nrow(dat)){

        lower = dat$age_lower[i] 
        upper = dat$age_upper[i] 
        trun_age = age_dat[age_dat$age >= lower & 
                                    age_dat$age < upper, ]

        # Sample
        inds = sample(1:nrow(trun_age), dat$sample_size[i], replace=FALSE)
        means[i] = mean(trun_age$worm[inds])
        prevs[i] = mean(trun_age$worm[inds] > 0)

    }

    return(c(means, prevs))
}

simulate_and_compare = function(i, abc_params, time_steps=100){
    # Simulate the raccoon model and compare it to the age intensity/prevalence
    # data.
    #
    # Parameters
    # ----------
    # abc_params : matrix of N rows and 4 columns where the columns correspond
    #               to the four parameters of interest drawn from their priors
    # time_steps : How many time steps to run the model.
    #
    # Returns
    # -------
    # : list
    #       `summary_stats`: Summary statistics for each simulation
    #       `params`: Parameter values for each simulation

    print(paste("Working on sim", i))

    abc_params_vect = abc_params[i, ]
    em = abc_params_vect["em"]
    ex = abc_params_vect["ex"]
    ec = abc_params_vect["ec"]
    rp = abc_params_vect["rp"]
    as = abc_params_vect["as"]
    ek = (1 - ex) / ex # Convert ex to ek

    params = get_simulation_parameters(TIME_STEPS=time_steps, 
                        RODENT_ENCOUNTER_PROB=rp, ENCOUNTER_K=ek, 
                        EGG_CONTACT=ec, ENCOUNTER_MEAN=em,
                        AGE_SUSCEPTIBILITY=as)

    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(params, init_arrays, 
                                  cull_params=NULL, 
                                  birth_control_params=NULL,
                                  worm_control_params=NULL,
                                  latrine_cleanup_params=NULL, 
                                  management_time=1000,
                                  print_it=FALSE)

    summary_stats = compare_to_data(all_res, time_steps)

    return(list(stats=summary_stats, params=c(em, ex, ec, rp, as)))

}

distances = function(stats){
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

    obs_stat = c(obs_dat$abundance, obs_dat$prev / 100)[non_zero_inds]
    scaled_obs = t(matrix((obs_stat - attr(scale(stats), "scaled:center")) / 
                        attr(scale(stats), "scaled:scale")))
    pca_obs = scaled_obs

    # loadings = eigen(cov(scale(stats)))$vectors
    # pca_obs = scaled_obs %*% loadings

    # Distances between observed and predicted data
    dists = as.matrix(dist(rbind(pca_stats, pca_obs), method="euclidean"))[sims + 1, 1:(sims)]

    return(dists)
}

resample_and_perturb = function(params, weights, num_samps, perturb_sds){
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
    #
    # Returns
    # -------
    # : a matrix of perturbed, resampled parameters

    inds = 1:nrow(params)
    new_param_inds = sample(inds, num_samps, replace=TRUE, prob=weights)
    new_params = params[new_param_inds, ]

    # Uniform perturbation
    new_params_perturbed = t(apply(new_params, 1, 
            function(x) x + perturb_sds*runif(length(x), min=-1, max=1)))

    # Normal perturbation
    # new_params_perturbed = t(apply(new_params, 1, 
    #         function(x) rnorm(length(x), mean=x, sd=perturb_sds)))

    # Ensure the these are giving the right results
    new_params_perturbed = check_params(new_params_perturbed, new_params)
    return(new_params_perturbed)

}

check_params = function(params_perturbed, params){
    # Makes sure each parameter is within its proper bounds
    # If a pertubation pushed it out set it to old parameters

    dp = dim(params)

    upper = matrix(rep(c(em=1000, ex=1, ec=5, rp=1, as=10), dp[1]), 
                                nrow=dp[1], ncol=dp[2], byrow=T)
    lower = matrix(rep(c(em=10, ex=0, ec=0.001, rp=0, as=0), dp[1]),
                                nrow=dp[1], ncol=dp[2], byrow=T)


    gt = params_perturbed >= upper
    lt = params_perturbed <= lower 
    params_perturbed[gt] = params[gt]
    params_perturbed[lt] = params[lt]

    return(params_perturbed)
}

get_particles = function(num_particles){
    # Sample num_particles parameter vectors from the prior distributions of 
    # each of the parameters.
    #
    # Returns
    # -------
    # : Matrix of parameters where each row is a parameter vector

    # Prior distributions on parameters
    em = runif(num_particles, min=10, max=1000)
    ex = runif(num_particles, min=0, max=1)
    ec = runif(num_particles, min=0.001, max=5)
    rp = runif(num_particles, min=0, max=1)
    as = runif(num_particles, min=0, max=10)

    params = cbind(em, ex, ec, rp, as)
    return(params)

}

particle_likelihood = function(particle){
    # Calculate the likelihood of a particle (vector of parameters) given the
    # prior distributions
    #
    # Returns
    # -------
    # : Matrix of parameters where each row is a parameter vector

    # Prior distributions on parameters
    em = dunif(particle['em'], min=10, max=1000)
    ex = dunif(particle['ex'], min=0, max=1)
    ec = dunif(particle['ec'], min=0.001, max=5)
    rp = dunif(particle['rp'], min=0, max=1)
    as = dunif(particle['as'], min=0, max=10)

    part_like = prod(c(em, ex, ec, rp, as))
    return(part_like)

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

### MAIN ###

STEPS = 4
num_particles = 48 # Number of parameter vectors to consider
current_params = get_particles(num_particles)
weight_array = list()
past_params = list()

# Perturbation standard deviations for each parameter
#perturb_sds = c(em=2, ex=0.02, ec=1, rp=0.02)

for(t in 1:STEPS){

    print(paste("Iteration", t))

    res = mclapply(1:num_particles, simulate_and_compare, current_params, mc.cores=4)

    # Extract the results
    stats = do.call(rbind, lapply(res, function(x) x$stats))
    params = do.call(rbind, lapply(res, function(x) x$params))

    # Standardize and convert the predicted statistics into PCA space
    dists = distances(stats)

    # Extract the alpha = 20% smallest distances
    lower = quantile(dists, c(.20))
    lowest_inds = dists < lower

    # Get the new params
    new_params = params[lowest_inds, ]
    past_params[[t]] = new_params

    # Calculate weights
    if(t == 1){
        weights = rep(1 / nrow(new_params), nrow(new_params))
    } else{
        weights = calculate_weights(weight_array[[t - 1]], prev_params, 
                                            new_params, perturb_sds)
    }
    weight_array[[t]] = weights

    # Perturb using uniform perturbation
    perturb_sds = apply(new_params, 2, sd) #(apply(new_params, 2, max) - apply(new_params, 2, min)) * 0.5 # From Filippi 2012
    current_params = resample_and_perturb(new_params, weights, num_particles, 
                            perturb_sds)

    prev_params = new_params
}

saveRDS(past_params, "past_params_concomitant.rds")


## Plot trajectories from the mean parameters ##
# best_params = data.table(readRDS("past_params_concomitant.rds")[[5]])
# best_params_m = data.frame(melt(best_params))

# median_params = apply(best_params, 2, quantile, 0.5)

# params = get_simulation_parameters(TIME_STEPS=100, 
#                         RODENT_ENCOUNTER_PROB=median_params['rp'], 
#                         ENCOUNTER_K=0.1,#(1 - median_params['ex']) / median_params['ex'], 
#                         EGG_CONTACT=median_params['ec'], 
#                         ENCOUNTER_MEAN=400,
#                         AGE_SUSCEPTIBILITY=10,
#                         AGE_EGG_RESISTANCE=4)#median_params['em'])

# init_arrays = get_init_arrays(params) # Load in init arrays
# all_res = full_simulation(params, init_arrays, 
#                                   cull_params=NULL, 
#                                   birth_control_params=NULL,
#                                   worm_control_params=NULL,
#                                   latrine_cleanup_params=NULL, 
#                                   management_time=1000,
#                                   print_it=TRUE)

# NUM_DRAWS = 100
# obs = read.csv(datasource) # Load observed data
# draw_results = list()

# for(i in 1:NUM_DRAWS){

#     pred = compare_to_data(all_res, 100)
#     prev = pred[9:16]
#     abund = pred[1:8]
#     sim_num = rep(i, 8)
#     age = obs$age_lower
#     tres = cbind(sim_num, abund, prev, age)

#     draw_results[[i]] = tres
# }

# draw_results = data.frame(do.call(rbind, draw_results))
# ggplot(draw_results, aes(x=age, y=abund, group=sim_num)) + geom_point() + 
#             geom_point(data=obs, aes(x=age_lower, y=abundance, color="red")) +
#             geom_errorbar(data=obs, aes(x=age_lower, ymin=abund_lower, ymax=abund_upper, color="red"))


# ggplot(draw_results, aes(x=age, y=prev, group=sim_num)) + geom_point() + 
#             geom_point(data=obs, aes(x=age_lower, y=prev/100, color="red")) + 
#             geom_errorbar(data=obs, aes(x=age_lower, ymin=prev_lower/100, ymax=prev_upper/100, color="red")) +
#             ylim(0, 1)









