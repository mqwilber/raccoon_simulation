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