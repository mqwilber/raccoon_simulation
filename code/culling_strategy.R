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

full_simulation = function(cull_params, birth_control_params, 
            worm_control_params, management_time, 
            prms, init_arrays, print_it=FALSE){
    # Runs the full individual based simulation
    #
    # Parameters
    # ----------
    # cull_params : list or NULL, parameters for culling strategies. See cull_strategy
    #    i.e. list(strategy="age", quota=5, overlap_threshold=0.5)
    # management_time : int, time in simulation where management begins
    # prms : list, simulation parameters
    # init_arrays : initial arrays for holding results
    # print_it : bool, if TRUE prints progress, otherwise False.


    # Save the init arrays to the current environment
    list2env(init_arrays, envir=environment())

    # Loop through time steps
    for(time in 2:(prms$TIME_STEPS + 1)){

        if(print_it){
            print(paste("Beginning step", time - 1, "of", prms$TIME_STEPS))
        }


        new_babies = 0
        babies_at_this_time_vect = array(NA, dim=dim(raccoon_worm_array)[2])

        # Calculate eggs remaining in each zone for this time step
        eggs_remaining = get_cumulative_egg_load(time - 1, eggproduction_array, 
                                            prms$EGG_DECAY)

        # If time is >= management time begin management
        if(time >= management_time){

            # TODO: Set this up so you don't comment
            cull_indices = integer(0)#get_cull_indices(cull_params,
                            #             raccoon_dead_alive_array[time - 1, ],
                            #             age_array[time - 1, ], 
                            #             human_overlap_through_time[[time - 1]])

            # Pick individuals for birth control
            birth_control_indices = get_birth_control_indices(birth_control_params,
                                                repro_able_vect)

            # If picked for birth control, can't have babies
            repro_able_vect[birth_control_indices] = 0 

            tworm_control_params = NULL #worm_control_params 
        } else{
            birth_control_indices = integer(0) 
            cull_indices = integer(0)
            tworm_control_params = NULL 
        }

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

                tot_racs = sum(raccoon_dead_alive_array[time, ], na.rm=T)
                new_babies_now = give_birth(age_now, time, tot_racs, 
                                                        repro_able_vect[rac],
                                                        prms$MONTH_AT_REPRO,
                                                        prms$FIRST_REPRO_AGE,
                                                        prms$LITTER_SIZE, prms$BETA)

                # Add new babies to array to assign human contacts later
                babies_at_this_time_vect[rac] = new_babies_now
                new_babies = new_babies + new_babies_now

                # 2. Kill old worms for each possible source of worms
                got_bait = picked_up_bait(overlap_now, tworm_control_params)


                previous_cohorts = infra_worm_array[[rac]][time - 1, 1:(time - 1)]
                new_cohort = kill_raccoon_worms(previous_cohorts,
                                                    prms$WORM_SURV_TRESH,
                                                    prms$WORM_SURV_SLOPE,
                                                    got_bait=got_bait)

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
                    worms_acquired = pick_up_rodents(prms$MOUSE_WORM_MEAN,
                                                     prms$RODENT_ENCOUNTER_PROB,
                                                     prms$LARVAL_WORM_INFECTIVITY,
                                                     eggs_remaining[zone_now],
                                                     prms$EGG_CONTACT)

                }

                infra_worm_array[[rac]][time, time] = worms_acquired
                raccoon_worm_array[time, rac] = sum(infra_worm_array[[rac]][time, ], na.rm=T)

                # Sum over all worm sources and add to total raccoon array
                # raccoon_worm_array[time, rac] = sum(all_worms_infra_array[[1]][[rac]][time, ], na.rm=T) + 
                #                                 sum(all_worms_infra_array[[2]][[rac]][time, ], na.rm=T)

                # 4. Disperse if raccoon is 6
                if(age_now == prms$DISPERSAL_AGE){
                    human_vect[rac] = assign_human_contacts(1)
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

        # Save the new human risk array
        human_overlap_through_time[[time]] = human_vect
        eggproduction_array[time, ] = assign_egg_production(
                                            raccoon_worm_array[time, ],
                                            human_vect,
                                            prms$ZONES)

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
                eggproduction_array=eggproduction_array))
} # End full simulation


run_and_extract_results = function(i, quota, management_time, 
                                    cull_params, birth_control_params, 
                                    worm_control_params, time_steps){

    print(paste("Simulation", i, "for quota", quota))

    params = get_simulation_parameters(TIME_STEPS=time_steps) # Load in simulation parameters
    init_arrays = get_init_arrays(params) # Load in init arrays

    cull_params$quota = quota
    all_res = full_simulation(cull_params, birth_control_params, 
                                worm_control_params, management_time,
                                params, init_arrays)

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
cull_params = NULL # list(strategy="random", quota=9, overlap_threshold=0.8)
birth_control_params = list(strategy="random", quota=20)
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

    params = get_simulation_parameters(TIME_STEPS=time_steps)#, K_CAPACITY=10, 
                                              #INIT_NUM_RACCOONS=20) # Load in simulation parameters
    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(cull_params, birth_control_params, 
                                    worm_control_params, management_time,
                                    params, init_arrays, print_it=TRUE)

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