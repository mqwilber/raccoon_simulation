## IBM Model Draft for Raccoon Parasites

## Authors: SW and MW

## Notes
##
## We are assuming an all female raccoon population and they are
## only reproducing females. 

library(parallel)
library(profvis)
source("raccoon_fxns.R") # Load in the helper functions
set.seed(1)

run_and_extract_results = function(i, quota, management_time, 
                                    cull_params, birth_control_params, 
                                    worm_control_params, time_steps){
    # Runs simulations and extracts results

    print(paste("Simulation", i, "for quota", quota))

    params = get_simulation_parameters(TIME_STEPS=time_steps) # Load in simulation parameters
    init_arrays = get_init_arrays(params) # Load in init arrays

    cull_params$quota = quota
    all_res = full_simulation(params, init_arrays,
                                cull_params=cull_params, 
                                birth_control_params=birth_control_params, 
                                worm_control_params=worm_control_params, 
                                management_time=management_time)

    # Extract min, max rac pop sizes
    pop_traj = rowSums(all_res$raccoon_dead_alive_array, na.rm=T)
    tot_steps = length(pop_traj)
    max_rac_pop = max(pop_traj[(tot_steps - 12):tot_steps])
    min_rac_pop = min(pop_traj[(tot_steps - 12):tot_steps])
    mean_rac_pop = mean(pop_traj[(tot_steps - 12):tot_steps])

    # Extract min, max, worm pop sizes
    worm_pop_traj = rowSums(all_res$raccoon_worm_array, na.rm=T)
    max_worm_pop = max(worm_pop_traj[(tot_steps - 12):tot_steps])
    min_worm_pop = min(worm_pop_traj[(tot_steps - 12):tot_steps])
    mean_worm_pop = mean(worm_pop_traj[(tot_steps - 12):tot_steps])

    # Human risk # UPDATE THIS
    human_risk = get_human_risk_metric(all_res$eggproduction_array, params$EGG_DECAY)
    max_human_risk = max(human_risk[(tot_steps - 12):tot_steps])
    min_human_risk = min(human_risk[(tot_steps - 12):tot_steps])
    mean_human_risk = mean(human_risk[(tot_steps - 12):tot_steps])


    return(c(min_rac_pop, mean_rac_pop, max_rac_pop, 
             min_worm_pop, mean_worm_pop, max_worm_pop,
             min_human_risk, mean_human_risk, max_human_risk))

}


## RUNNING SIMULATION ###

# Simulation parameters
cull_params = list(strategy="random", quota=9, overlap_threshold=0.8)
birth_control_params = list(strategy="random", quota=20, overlap_threshold=0.7)

# Worm control, this one might not be quota based.  However, for comparison
# purposes it might make sense to implement this one as a quota as well...
worm_control_params = NULL #list(strategy="random", distribution=0.5)


quotas = 0:10#0:5
SIMS = 50
management_time = 10
time_steps = 50
col_names = c("min_rac_pop", "mean_rac_pop", "max_rac_pop", 
             "min_worm_pop", "mean_worm_pop", "max_worm_pop",
             "min_human_risk", "mean_human_risk", "max_human_risk")
single_sim = TRUE # IF TRUE JUST RUNS A SINGLE SIMULATION 


if(single_sim){ # Run a single simulation

    params = get_simulation_parameters(TIME_STEPS=time_steps)
    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(params, init_arrays, 
                              cull_params=NULL, 
                              birth_control_params=birth_control_params, 
                              worm_control_params=NULL, 
                              management_time=management_time,
                              print_it=TRUE)

} else{ # Run a full simulation

    # Different strategies
    strategies = list("human"=c(0.5), "age"=c(12), "random"=c(1))

    for(strategy in names(strategies)){ # Loop through different strategies

        print(paste("Beginning", strategy))

        controls = strategies[[strategy]]

        for(control in controls){ # Loops through different ages of overlaps

            cull_params$strategy = strategy
            cull_params$overlap_threshold = control
            cull_params$age = control

            # Arrays to hold sim results
            sim_mean_results = list()
            sim_var_results = list()

            for(j in 1:length(quotas)){ # Looping through quotas
                
                # Parallelize within each quota
                sim_vals = mclapply(1:SIMS, run_and_extract_results, 
                                                    quotas[j], management_time, 
                                                    cull_params, birth_control_params,
                                                    worm_control_params, time_steps, 
                                                    mc.cores=4)

                cull_matrix = do.call(rbind, sim_vals)

                cull_means = colMeans(cull_matrix)
                cull_vars = apply(cull_matrix, 2, sd)
                names(cull_means) = col_names
                names(cull_vars) = col_names
                sim_mean_results[[j]] = cull_means
                sim_var_results[[j]] = cull_vars

            }

            saveRDS(sim_mean_results, 
                        paste("../results/cull_", strategy, "/sim_mean_results_", 
                               strategy, control, ".rds", sep=""))
            saveRDS(sim_var_results, 
                        paste("../results/cull_", strategy, "/sim_var_results_", 
                               strategy, control, ".rds", sep=""))

            print(paste("Analysis complete for", strategy, control))
        }

    }
}