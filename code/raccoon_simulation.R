## IBM Model Draft for Raccoon Parasites

## Date: 2 / 1 / 2016

## Authors: SW and MW

## Notes
##
## We are assuming (right now) an all female raccoon population and they are
## only reproducing females.  We will modify this later (2-26-2016)
##
##


## TODO:
## 2. Carrying capacity on raccoon population. Looking at literature
## - Thought are to kill babies based on density
## 3. Look at actual population data to get parameter values and functional
## forms.
## 4. Document to discuss what we are doing with each function and why we made
## each decision.

source("raccoon_fxns.R") # Load in the helper functions
source("raccoon_parameters.R") # Load in the parameters
source("raccoon_init_arrays.R") # Load in the initial arrays

# Loop through time steps
for(time in 2:(TIME_STEPS + 1)){

    new_babies = 0

    # Loop through raccoons
    for(rac in 1:dim(raccoon_worm_array)[2]){

        alive_then = raccoon_dead_alive_array[time - 1, rac]


        if(alive_then == 1){

            # 1. Raccoon dies (Killing at the end of the time interval)

            age_now = initial_age_vector[rac] +
                        sum(raccoon_dead_alive_array[, rac], na.rm=T)

            alive_now = kill_my_raccoon(raccoon_worm_array[time - 1, rac],
                                            age_now,
                                            DEATH_THRESHOLD, PATHOGENICITY,
                                            BABY_DEATH,
                                            INTRINSIC_DEATH_RATE,
                                            RANDOM_DEATH_PROB)
            raccoon_dead_alive_array[time, rac] = alive_now


            if(alive_now){

                # 1. Give birth if appropriate
                # Calculate age. Minus 1 to account for alive dead decision

                age_array[time, rac] = age_now
                new_babies = new_babies + give_birth(age_now, time,
                                                        MONTH_AT_REPRO,
                                                        FIRST_REPRO_AGE,
                                                        LITTER_SIZE)

                # 2. Deposit Eggs (ignoring for now)

                # 3. Kill old worms.
                previous_cohorts = infra_worm_array[[rac]][time - 1, 1:(time - 1)]
                new_cohort = kill_raccoon_worms(previous_cohorts,
                                                WORM_SURV_TRESH,
                                                WORM_SURV_SLOPE)
                # Assign previous cohort to current cohort
                infra_worm_array[[rac]][time, 1:(time - 1)] = new_cohort

                # 4. Pick eggs
                worms_acquired = pick_up_eggs(ENCOUNTER_PROB, ENCOUNTER_MEAN,
                                              INFECTIVITY, RESISTANCE,
                                             raccoon_worm_array[time - 1, rac])

                infra_worm_array[[rac]][time, time] = worms_acquired
                raccoon_worm_array[time, rac] = sum(infra_worm_array[[rac]][time, ], na.rm=T)#raccoon_worm_array[time - 1, rac] + worms_acquired

            } else { # Raccoon dead and its worms die
                raccoon_worm_array[time, rac] = NA
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
                                           TIME_STEPS)

        new_babies_vect = updated_arrays$new_babies_vect
        initial_age_vector = updated_arrays$initial_age_vector
        age_array = updated_arrays$age_array
        raccoon_dead_alive_array = updated_arrays$raccoon_dead_alive_array
        raccoon_worm_array = updated_arrays$raccoon_worm_array
        infra_worm_array = updated_arrays$infra_worm_array

    }

}

