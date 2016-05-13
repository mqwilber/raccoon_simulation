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

# Loop through time steps
for(time in 2:(TIME_STEPS + 1)){

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

            alive_now = kill_my_raccoon(raccoon_worm_array[time - 1, rac],
                                            age_now,
                                            DEATH_THRESHOLD, PATHOGENICITY,
                                            BABY_DEATH,
                                            INTRINSIC_DEATH_RATE,
                                            RANDOM_DEATH_PROB, OLD_DEATH)
            raccoon_dead_alive_array[time, rac] = alive_now


            if(alive_now){

                # 1. Give birth if appropriate
                # Calculate age. Minus 1 to account for alive dead decision

                age_array[time, rac] = age_now

                tot_racs = sum(raccoon_dead_alive_array[time, ], na.rm=T)
                new_babies_now = give_birth(age_now, time, tot_racs,
                                                        MONTH_AT_REPRO,
                                                        FIRST_REPRO_AGE,
                                                        LITTER_SIZE, BETA)

                # Add to new babies to array to assign human contacts later
                babies_at_this_time_vect[rac] = new_babies_now
                new_babies = new_babies + new_babies_now

                # 2. Deposit Eggs (ignoring for now)

                # 3. Kill old worms.
                previous_cohorts = infra_worm_array[[rac]][time - 1, 1:(time - 1)]
                new_cohort = kill_raccoon_worms(previous_cohorts,
                                                WORM_SURV_TRESH,
                                                WORM_SURV_SLOPE)
                # Assign previous cohort to current cohort
                infra_worm_array[[rac]][time, 1:(time - 1)] = new_cohort

                # 4. Pick eggs
                worms_acquired = pick_up_eggs(ENCOUNTER_MEAN,
                                            ENCOUNTER_K,
                                            INFECTIVITY, RESISTANCE,
                                            previous_prevalence[1:(time - 1)],
                                            raccoon_worm_array[time - 1, rac],
                                            EGG_DECAY, ENCOUNTER_PARAMS)

                infra_worm_array[[rac]][time, time] = worms_acquired
                raccoon_worm_array[time, rac] = sum(infra_worm_array[[rac]][time, ], na.rm=T)#raccoon_worm_array[time - 1, rac] + worms_acquired

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
                                           age_array, infra_worm_array,
                                           human_array,
                                           babies_at_this_time_vect,
                                           TIME_STEPS)

        new_babies_vect = updated_arrays$new_babies_vect
        initial_age_vector = updated_arrays$initial_age_vector
        age_array = updated_arrays$age_array
        raccoon_dead_alive_array = updated_arrays$raccoon_dead_alive_array
        raccoon_worm_array = updated_arrays$raccoon_worm_array
        infra_worm_array = updated_arrays$infra_worm_array
        human_array = updated_arrays$human_array

    }

    # Save the new human risk array
    human_risk_through_time[[time]] = human_array


    ## TESTING CODE ##
    # if(time == 10){
    #     ind = !is.na(raccoon_worm_array[time, ])
    #     raccoon_worm_array[time, ind] = 0
    # }

}



