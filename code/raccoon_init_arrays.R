## Initial arrays for raccoon simulation ###

# Set up the raccoon arrays
raccoon_worm_array = array(NA, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))
raccoon_dead_alive_array = array(NA, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))
initial_age_vector = rep(1, INIT_NUM_RACCOONS)
age_array = array(NA, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))
age_array[1, ] = initial_age_vector
new_babies_vect = array(NA, dim=TIME_STEPS + 1)
new_babies_vect[1] = 0

# Initialize all arrays
raccoon_worm_array[1, ] = 0 # Initialize all raccoons with 0 worms
raccoon_dead_alive_array[1, ] = 1 # All raccoons are alive

# Set up worm arrays.  This array hold the cohort of worms for each raccoon
# so that we can track the age-dependent death of the worms in the raccoons
infra_worm_array = lapply(1:INIT_NUM_RACCOONS,
                function(x) array(NA, dim=c(TIME_STEPS + 1, TIME_STEPS + 1)))

# Initialize worm arrays
for(i in 1:length(infra_worm_array)){
    infra_worm_array[[i]][1, 1] = raccoon_worm_array[1, i]

}