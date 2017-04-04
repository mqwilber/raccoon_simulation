## IBM Model Draft for Raccoon Parasites

## Authors: SW and MW

## Notes
##
## We are assuming an all female raccoon population and they are
## only reproducing females. 

library(parallel)
# library(profvis)
source("raccoon_fxns.R") # Load in the helper functions
set.seed(1)

run_and_extract_results = function(i, quota, management_time, 
                                    cull_params, birth_control_params, 
                                    worm_control_params, time_steps){
    # Runs simulations and extracts results

    print(paste("Simulation", i, "for quota", quota))

    params = get_simulation_parameters(TIME_STEPS=time_steps) # Load in simulation parameters
    init_arrays = get_init_arrays(params) # Load in init arrays

    all_res = full_simulation(params, init_arrays,
                                cull_params=cull_params, 
                                birth_control_params=birth_control_params, 
                                worm_control_params=worm_control_params, 
                                management_time=management_time)

    # Extract min, max rac pop sizes
    decrement = 24 # Use the first two years of the simulation
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
    human_risk = get_human_risk_metric(all_res$eggproduction_array, params$EGG_DECAY)
    max_human_risk = max(human_risk[(tot_steps - decrement):tot_steps])
    min_human_risk = min(human_risk[(tot_steps - decrement):tot_steps])
    mean_human_risk = mean(human_risk[(tot_steps - decrement):tot_steps])

    # Prevalence measures
    prev_traj = get_prevalence(all_res$raccoon_worm_array)
    max_prev = max(prev_traj[(tot_steps - decrement):tot_steps])
    min_prev = min(prev_traj[(tot_steps - decrement):tot_steps])
    mean_prev = mean(prev_traj[(tot_steps - decrement):tot_steps])

    # Prevalence measures
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


## RUNNING SIMULATION ###

# Simulation parameters
cull_params = list(strategy="age", quota=20, overlap_threshold=0.8, age=12)
birth_control_params = list(strategy="random", quota=10, overlap_threshold=0.7)
worm_control_params = list(strategy="random", quota=10000, overlap_threshold=0.8)

# Worm control, this one might not be quota based.  However, for comparison
# purposes it might make sense to implement this one as a quota as well...

SIMS = 50
management_time = 100
time_steps = 400
col_names = c("min_rac_pop", "mean_rac_pop", "max_rac_pop", 
             "min_worm_pop", "mean_worm_pop", "max_worm_pop",
             "min_human_risk", "mean_human_risk", "max_human_risk",
	         "min_prev", "mean_prev", "max_prev",
	         "min_intensity", "mean_intensity", "max_intensity")
single_sim = TRUE # IF TRUE JUST RUNS A SINGLE SIMULATION 


if(single_sim){ # Run a single simulation

    management_time = 50
    time_steps = 180
    params = get_simulation_parameters(TIME_STEPS=time_steps)
    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(params, init_arrays, 
                              cull_params=cull_params, 
                              birth_control_params=NULL,
                              worm_control_params=NULL, 
                              management_time=management_time,
                              print_it=TRUE)

} else { # Run a full simulation
    
    management_scenarios = 
                list(cull_only=list(
                               cull_params=
                                    list("human"=c(0.1, 0.5, 0.9), 
                                         "age"=c(3, 12, 25), 
                                         "random"=c(1)),
                               birth_control_params=NULL,
                               worm_control_params=NULL),

                     birth_control_only=list(
                                cull_params=NULL,
                                birth_control_params=
                                     list("human"=c(0.1, 0.5, 0.9),
                                          "random"=c(1)), 
                                worm_control_params=NULL),

                     worm_control_only=list(
                                cull_params=NULL,
                                birth_control_params=NULL,
                                worm_control_params=
                                     list("human"=c(0.1, 0.5, 0.9),
                                          "random"=c(1))))

                # list(

                #      worm_control_only=list(
                #                 cull_params=NULL,
                #                 birth_control_params=NULL,
                #                 worm_control_params=
                #                      list("human"=seq(0.1, 0.9, length=9),
                #                           "random"=c(1)))) 


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
        } else{
            quotas = c(0:10, 20, 50, 100, 200)
        }

        # Loop through different management strategies within a current action
        for(strategy in names(scenario[[current_action]])) {

            print(paste("Beginning", strategy, "for", current_action))

            controls = scenario[[current_action]][[strategy]]



            for(control in controls){ # Loops controls on the strategy

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
                                                        worm_control_params, time_steps, 
                                                        mc.cores=4)

                    cull_matrix = do.call(rbind, sim_vals)

                    cull_means = colMeans(cull_matrix, na.rm=TRUE)
                    cull_vars = apply(cull_matrix, 2, sd, na.rm=TRUE)
                    names(cull_means) = col_names
                    names(cull_vars) = col_names
                    sim_mean_results[[j]] = cull_means
                    sim_var_results[[j]] = cull_vars

                }

                save_prefix = strsplit(current_action, "params")[[1]]
                folder = "" #"../results/"
                saveRDS(sim_mean_results, 
                            paste(folder, save_prefix, strategy, "/sim_mean_results_", 
                                   strategy, control, ".rds", sep=""))
                saveRDS(sim_var_results, 
                            paste(folder, save_prefix, strategy, "/sim_var_results_", 
                                   strategy, control, ".rds", sep=""))

                print(paste("Analysis complete for", current_action, strategy, control))
            }

        }
    }
} # End else
