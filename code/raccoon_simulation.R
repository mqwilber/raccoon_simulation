## IBM Model Draft for Raccoon Parasites

## Date: 2 / 1 / 2016

## Authors: SW and MW


## TODOS: 1. Worm death 2. Raccoon birth and age structure 3. Infectivity with age
source("raccoon_fxns.R")

## Parameters

## Raccoon parameters
INIT_NUM_RACCOONS = 100
DEATH_PROB = 0
DEATH_THRESHOLD = 300 # Worms at which Raccoon is dead
PATHOGENICITY = (1 - DEATH_PROB) / DEATH_THRESHOLD
K = 50 # Karrying capacity --> exponential growth with threshold


# PARASITE PARAMETERS
ENVIRONMENTAL_POOL = 0
EGG_PRODUCTION_PER_WORM = 100 # eggs per month (note this is very wrong)
ENCOUNTER_PROB = 1 # 100% probability that
ENCOUNTER_MEAN = 100 # Probability of a raccoon and egg coming into contact
INFECTIVITY = 0.5 # Probability of infectivity
EGG_LOSS_PROB =  .14 # Probability that an egg survives of 7 months

# Time parameters: Each time step is a month
TIME_STEPS = 12

# Set up the raccoon array
raccoon_worm_array = array(10, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))
raccoon_dead_alive_array = array(10, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))

# Initialize all arrays
raccoon_worm_array[1, ] = 0 # Initialize all raccoons with 10 worms
raccoon_dead_alive_array[1, ] = 1 # All raccoons are alive

# Loop through time steps
for(time in 2:(TIME_STEPS + 1)){

    # Loop through raccoons
    for(rac in 1:dim(raccoon_worm_array)[2]){

        alive_then = raccoon_dead_alive_array[time - 1, rac]

        if(alive_then == 1){

            # 1. Raccoon dies (Killing at the end of the time interval)
            alive_now = kill_my_raccoon(raccoon_worm_array[time - 1, rac],
                                                    DEATH_PROB, PATHOGENICITY)
            raccoon_dead_alive_array[time, rac] = alive_now

            if(alive_now){
                # 1. Deposit Eggs (ignoring for now)

                # 2. Pick eggs
                worms_acquired = pick_up_eggs(ENCOUNTER_PROB, ENCOUNTER_MEAN, INFECTIVITY)
                raccoon_worm_array[time, rac] = raccoon_worm_array[time - 1, rac] + worms_acquired

            } else { # Raccoon dead and its worms die
                raccoon_worm_array[time, rac] = NA
            }


        } else{ # This keeps deads dead

            raccoon_worm_array[time, rac] = NA
            raccoon_dead_alive_array[time, rac] = 0

        }
    }
}

