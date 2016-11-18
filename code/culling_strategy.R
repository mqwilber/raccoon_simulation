## IBM Model Draft for Raccoon Parasites

## Date: 2 / 1 / 2016

## Authors: SW and MW

## Notes
##
## We are assuming (right now) an all female raccoon population and they are
## only reproducing females.  We will modify this later (2-26-2016)
##
##

source("raccoon_fxns.R") # Load in the helper functions
source("raccoon_parameters.R") # Load in the parameters
source("raccoon_init_arrays.R") # Load in the initial arrays
set.seed(1)

cull_params = NULL #list(strategy="age", cull_prob=0.9, overlap_threshold=0.5)
birth_control_params = NULL #list(strategy="random", distribution=0.9)
worm_control_params = NULL #list(strategy="random", distribution=0.5)
management_time = 10


# Loop through time steps
for(time in 2:(TIME_STEPS + 1)){

    print(paste("Beginning step", time - 1, "of", TIME_STEPS))

    if(time == management_time){
        birth_control_params = NULL #list(strategy="random", distribution=0.2)
        cull_params = list(strategy="random", cull_prob=0.05, overlap_threshold=0.5)
        worm_control_params = NULL #list(strategy="human", distribution=1, 
                                            #overlap_threshold=0.99) #list(strategy="random", distribution=0)
    }

    new_babies = 0
    babies_at_this_time_vect = array(NA, dim=dim(raccoon_worm_array)[2])
    previous_prevalence = get_prevalence(raccoon_worm_array)

    # Loop through raccoons
    for(rac in 1:dim(raccoon_worm_array)[2]){

        alive_then = raccoon_dead_alive_array[time - 1, rac]

        if(alive_then == 1){

            # 1. Raccoon dies (Killing at the beginning of the time interval)

            age_now = initial_age_vector[rac] +
                        sum(raccoon_dead_alive_array[, rac], na.rm=T)

            overlap_now = human_risk_through_time[[time - 1]][rac]

            alive_now = kill_my_raccoon(raccoon_worm_array[time - 1, rac],
                                            age_now, overlap_now, 
                                            DEATH_THRESHOLD, PATHOGENICITY,
                                            BABY_DEATH,
                                            INTRINSIC_DEATH_RATE,
                                            RANDOM_DEATH_PROB, OLD_DEATH,
                                            cull_params=cull_params)

            raccoon_dead_alive_array[time, rac] = alive_now


            if(alive_now){

                # 1. Give birth if appropriate

                age_array[time, rac] = age_now

                tot_racs = sum(raccoon_dead_alive_array[time, ], na.rm=T)
                new_babies_now = give_birth(age_now, time, tot_racs,
                                                        MONTH_AT_REPRO,
                                                        FIRST_REPRO_AGE,
                                                        LITTER_SIZE, BETA,
                                    birth_control_params=birth_control_params)

                # Add new babies to array to assign human contacts later
                babies_at_this_time_vect[rac] = new_babies_now
                new_babies = new_babies + new_babies_now

                # 2. Deposit Eggs (ignoring for now)

                # 3. Kill old worms for each possible source of worms
                got_bait = picked_up_bait(overlap_now, worm_control_params)

                for(tw in 1:length(all_worms_infra_array)){

                    previous_cohorts = all_worms_infra_array[[tw]][[rac]][time - 1, 1:(time - 1)]
                    new_cohort = kill_raccoon_worms(previous_cohorts,
                                                    WORM_SURV_TRESH,
                                                    WORM_SURV_SLOPE,
                                                    got_bait=got_bait)

                    # Assign previous cohort to current cohort
                    all_worms_infra_array[[tw]][[rac]][time, 1:(time - 1)] = new_cohort

                }

                # 4. Pick up eggs or pick up rodents depending on age

                if(age_now <= AGE_EGG_RESISTANCE){ # Only pick up worms from eggs

                    worms_acquired = pick_up_eggs(ENCOUNTER_MEAN,
                                            ENCOUNTER_K,
                                            INFECTIVITY, RESISTANCE,
                                            previous_prevalence[1:(time - 1)],
                                            raccoon_worm_array[time - 1, rac],
                                            EGG_DECAY, ENCOUNTER_PARAMS)

                    all_worms_infra_array[[1]][[rac]][time, time] = worms_acquired
                    all_worms_infra_array[[2]][[rac]][time, time] = 0

                } else{ # Only pick up worms from rodent
                    worms_acquired = pick_up_rodents(MOUSE_WORM_MEAN,
                                                     MOUSE_WORM_AGG,
                                                     RODENT_ENCOUNTER_PROB,
                                                     LARVAL_WORM_INFECTIVITY)

                    all_worms_infra_array[[1]][[rac]][time, time] = 0
                    all_worms_infra_array[[2]][[rac]][time, time] = worms_acquired
                }

                # Sum over all worm sources and add to total raccoon array
                raccoon_worm_array[time, rac] = sum(all_worms_infra_array[[1]][[rac]][time, ], na.rm=T) + 
                                                sum(all_worms_infra_array[[2]][[rac]][time, ], na.rm=T)

                # sum(unlist(lapply(all_worms_infra_array, 
                #             function(wa) sum(wa[[rac]][time, ], na.rm=T))))


                # 5. Disperse if raccoon is 6
                if(age_now == DISPERSAL_AGE){
                    human_array[rac] = assign_human_contacts(1)
                }

            } else { # Raccoon dead and its worms die
                raccoon_worm_array[time, rac] = NA
                human_array[rac] = NA
            }


        } else{ # This keeps deads dead

            raccoon_worm_array[time, rac] = NA
            raccoon_dead_alive_array[time, rac] = 0

        }
    }

    # Update the various vectors if babies were born
    if(new_babies > 0){

        updated_arrays = update_arrays(time, new_babies, new_babies_vect,
                                           initial_age_vector,
                                           raccoon_dead_alive_array,
                                           raccoon_worm_array,
                                           age_array, all_worms_infra_array,
                                           human_array,
                                           babies_at_this_time_vect,
                                           TIME_STEPS)

        new_babies_vect = updated_arrays$new_babies_vect
        initial_age_vector = updated_arrays$initial_age_vector
        age_array = updated_arrays$age_array
        raccoon_dead_alive_array = updated_arrays$raccoon_dead_alive_array
        raccoon_worm_array = updated_arrays$raccoon_worm_array
        all_worms_infra_array = updated_arrays$all_worms_infra_array
        human_array = updated_arrays$human_array

    }

    # Save the new human risk array
    human_risk_through_time[[time]] = human_array

}

