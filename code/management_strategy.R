## IBM Model for Raccoon Parasites
##
## Description
## -----------
## The following script simulates the raccoon IBM model.  This script
## can three types of simulations
##
## 1. A single simulation of the IBM.  This is achieved by setting the
## simulation_types = c("single sim").  All the simulation results are then stored
## in an object called all_res. This is helpful for testing out different
## strategies to see how they are working.
##
## 2. Bulk simulations of single management strategies.  This is achieved by
## setting simulation_types = c("one management").
##
## 3. Bulk simulations of dual management strategies. This is achieved by 
## setting simulation_types = c("two management")
##
## Any combination of these options can also be used.
##
## Authors: Mark Wilber and Sara Weinstein 

library(parallel)
# library(profvis)
source("raccoon_fxns.R") # Load in the helper functions
#set.seed(1)


# A vector that can contain 1 or more of the three analyses
# "single sim", "one management", "two_management"
simulation_types = c("single sim")#c("one management", "two management") #"single sim" # Either single sim, one management bulk, or two management bulk 


results_dir = "../test_dir"
fitted_path = "fit_abc_results.rds"
cores = 8
SIMS = 40
management_time = 100
time_steps = 300
use_fitted = FALSE


# Load the ABC-SMC fitted transmission parameters
if(use_fitted){
    fitted_params = get_fitted_params(fitted_path)
} else{
    fitted_params = list()
}


col_names = c("min_rac_pop", "mean_rac_pop", "max_rac_pop", 
         "min_worm_pop", "mean_worm_pop", "max_worm_pop",
         "min_human_risk", "mean_human_risk", "max_human_risk",
         "min_prev", "mean_prev", "max_prev",
         "min_intensity", "mean_intensity", "max_intensity")


if(any(simulation_types == "single sim")) { # Run a single simulation

    # Simulation parameters
    cull_params = list(strategy="human", quota=500, overlap_threshold=0, age=12)
    birth_control_params = list(strategy="random", quota=1000000, overlap_threshold=0.5)
    worm_control_params = list(strategy="random", quota=10000, overlap_threshold=0.8)
    latrine_cleanup_params = list(strategy="human", overlap_threshold=0.3, quota=0.5)

    management_time = 100
    time_steps = 80
    fitted_params[['TIME_STEPS']] = time_steps
    fitted_params$ZONES = 1
    params = do.call(get_simulation_parameters, fitted_params)
    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(params, init_arrays, 
                              cull_params=NULL, 
                              birth_control_params=NULL,
                              worm_control_params=NULL,
                              latrine_cleanup_params=NULL, 
                              management_time=management_time,
                              print_it=TRUE)

} 

if(any(simulation_types == "one management")) {
    
    SIMS = SIMS # Number of simulations for combination
    management_time = management_time # Time at which management starts 
    time_steps = time_steps # Length of a single simulation
    
    management_scenarios = 
                list(cull_only=list(
                               cull_params=
                                    list("human"=c(0.1, 0.5, 0.9), 
                                         "age"=c(3, 12, 25), 
                                         "random"=c(1)),
                               birth_control_params=NULL,
                               worm_control_params=NULL,
                               latrine_cleanup_params=NULL),

                     birth_control_only=list(
                                cull_params=NULL,
                                birth_control_params=
                                     list("human"=c(0.1, 0.5, 0.9),
                                          "random"=c(1)), 
                                worm_control_params=NULL,
                                latrine_cleanup_params=NULL),

                     worm_control_only=list(
                                cull_params=NULL,
                                birth_control_params=NULL,
                                worm_control_params=
                                     list("human"=c(0.1, 0.5, 0.9),
                                          "random"=c(1)),
                                latrine_cleanup_params=NULL),

                    latrine_cleanup_only=list(
                                cull_params=NULL,
                                birth_control_params=NULL,
                                worm_control_params=NULL,
                                latrine_cleanup_params=
                                    list("human"=c(0.1, 0.5, 0.9))))


    for(scenario_nm in names(management_scenarios)) {

        scenario = management_scenarios[[scenario_nm]]

        # Assign scenario list to environment: assigns cull_params, birth_control_params, worm_control_params
        list2env(scenario, envir=environment())

        # Figure out which action is happening: cull, birth control, or worm control
        for(action in names(scenario)) {

            if(length(scenario[[action]]) > 0){
                current_action = action
            }
        }

        # Set quotas for the correct scenario
        if(scenario_nm == "worm_control_only"){
            quotas = c(0, 1, 10, 100, 1000, 10000, 100000)
        } else if(scenario_nm == "latrine_cleanup_only"){
            quotas = c(0, 0.1, 0.5, 1) # Quotas are cleaning efficiencies of latrines (proportions between 0 and 1)
        } else{
           quotas = c(0, 10, 20, 50, 100, 500, 1000, 5000) 
        }

        # Loop through different management strategies within a current action
        for(strategy in names(scenario[[current_action]])) {

            print(paste("Beginning", strategy, "for", current_action))

            controls = scenario[[current_action]][[strategy]]

            for(control in controls) { # Loops controls on the strategy

                eval(parse(text=paste(current_action, "$strategy=strategy", sep="")))
                eval(parse(text=paste(current_action, "$overlap_threshold=control", sep="")))
                eval(parse(text=paste(current_action, "$age=control", sep="")))

                # Arrays to hold sim results
                sim_mean_results = list()
                sim_var_results = list()

                for(j in 1:length(quotas)) { # Looping through quotas
                    
                    # Parallelize within each quota

                    eval(parse(text=paste(current_action, "$quota=", quotas[j], sep="")))
                    sim_vals = mclapply(1:SIMS, run_and_extract_results, 
                                                        quotas[j], management_time, 
                                                        cull_params, birth_control_params,
                                                        worm_control_params, latrine_cleanup_params, 
                                                        time_steps, fitted_params, 
                                                        mc.cores=cores)
                    
                    # Cull is a stand in for any management...not just culling
                    cull_matrix = do.call(rbind, sim_vals)

                    cull_means = colMeans(cull_matrix, na.rm=TRUE)
                    cull_vars = apply(cull_matrix, 2, sd, na.rm=TRUE)
                    names(cull_means) = col_names
                    names(cull_vars) = col_names
                    sim_mean_results[[j]] = cull_means
                    sim_var_results[[j]] = cull_vars

                }

                save_prefix = strsplit(current_action, "params")[[1]]

                folder = file.path(results_dir, paste(save_prefix, strategy, sep=""))

                dir.create(folder, showWarnings=FALSE)
                filenm_mean = paste("sim_mean_results_", 
                                   strategy, control, ".rds", sep="")
                filenm_var = paste("sim_var_results_", 
                                   strategy, control, ".rds", sep="")
                saveRDS(sim_mean_results, file.path(folder, filenm_mean))

                saveRDS(sim_var_results, file.path(folder, filenm_mean))

                print(paste("Analysis complete for", current_action, strategy, control))
            }

        }
    } 
}

if(any(simulation_types == "two management")) {

    # Just focus on high quota, for three different human overlap

    SIMS = SIMS # Number of simulations for combination
    management_time = management_time # Time at which management starts 
    time_steps = time_steps # Length of a single simulation

    human_overlaps = c(0.1, 0.5, 0.9)

    # List of different dual managements
    two_management_strats = list(

        dual_latrine_cleanup_and_bait = list(
            birth_control_params=NULL,
            cull_params = NULL,
            worm_control_params=list(strategy="human", quota=100000, 
                                        overlap_threshold=0.1),
            latrine_cleanup_params=list(strategy="human", overlap_threshold=0.1, 
                                        cleanup_efficiency=1)),

        dual_cull_and_bait = list(
            birth_control_params=NULL,
            cull_params = list(strategy="human", quota=5000, 
                                        overlap_threshold=0.1),
            worm_control_params=list(strategy="human", quota=100000, 
                                        overlap_threshold=0.1),
            latrine_cleanup_params=NULL),

        dual_birth_control_and_bait = list(
            birth_control_params=list(strategy="human", quota=5000, 
                                        overlap_threshold=0.1),
            cull_params = NULL,
            worm_control_params=list(strategy="human", quota=100000, 
                                        overlap_threshold=0.1),
            latrine_cleanup_params=NULL)
        )
    
    # Loop through different double strategies
    for(strategy_nm in names(two_management_strats)){

        print(paste("Beginning analysis for", strategy_nm))

        strategies = two_management_strats[[strategy_nm]]

        # Arrays to hold sim results
        sim_mean_results = list()
        sim_var_results = list()

        # Loop through different human overlaps
        for(j in 1:length(human_overlaps)){

            human_overlap = human_overlaps[j]

            # Assign the human overlap in each of the strategies
            for(strat_type_nm in names(strategies)){

                if(length(strategies[[strat_type_nm]]) > 0) {

                    strategies[[strat_type_nm]][['overlap_threshold']] = human_overlap
                }
            }

            # Assign all updated strategies to the environment
            list2env(strategies, envir=environment())

            sim_vals = mclapply(1:SIMS, run_and_extract_results, 
                                    "HIGH", management_time, 
                                    cull_params, birth_control_params,
                                    worm_control_params, latrine_cleanup_params, 
                                    time_steps, 
                                    mc.cores=cores)

            cull_matrix = do.call(rbind, sim_vals)

            cull_means = colMeans(cull_matrix, na.rm=TRUE)
            cull_vars = apply(cull_matrix, 2, sd, na.rm=TRUE)
            names(cull_means) = col_names
            names(cull_vars) = col_names
            sim_mean_results[[paste("overlap=", human_overlap, sep="")]] = cull_means
            sim_var_results[[paste("overlap=", human_overlap, sep="")]] = cull_vars

        } # End human overlap for loop

        # Save the results for the particular strategy
        res_dir = file.path(results_dir, strategy_nm)
        dir.create(res_dir, showWarnings=FALSE)

        saveRDS(sim_mean_results, file.path(res_dir, 
                        "sim_mean_results_high_quota.rds"))
        saveRDS(sim_var_results, file.path(res_dir, 
                        "sim_var_results_high_quota.rds"))
        print(paste("Completed analysis for", strategy_nm))

    } # End double strategy for loop

} # End if/else
