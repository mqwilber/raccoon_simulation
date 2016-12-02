## Initial arrays for raccoon simulation ###

# Set up the raccoon arrays

raccoon_dead_alive_array = array(NA, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))
initial_age_vector = rep(10, INIT_NUM_RACCOONS)
age_array = array(NA, dim=c(TIME_STEPS + 1, INIT_NUM_RACCOONS))
age_array[1, ] = initial_age_vector
human_array = assign_human_contacts(INIT_NUM_RACCOONS)
human_risk_through_time = list()
human_risk_through_time[[1]] = human_array


new_babies_vect = array(NA, dim=TIME_STEPS + 1)
new_babies_vect[1] = 0

# Initialize all arrays
# Seed a raccoon with some worms
# raccoon_worm_array[1, ] = 10 # Seeding worms
raccoon_dead_alive_array[1, ] = 1 # All raccoons are alive

# Set up worm arrays.  This array hold the cohort of worms for each raccoon
# so that we can track the age-dependent death of the worms in the raccoons
infra_mouse_worm_array = lapply(1:INIT_NUM_RACCOONS,
                function(x) array(NA, dim=c(TIME_STEPS + 1, TIME_STEPS + 1)))

infra_nonmouse_worm_array = lapply(1:INIT_NUM_RACCOONS,
                function(x) array(NA, dim=c(TIME_STEPS + 1, TIME_STEPS + 1)))


# # Initialize worm arrays
# for(i in 1:length(infra_worm_array)){
#     infra_worm_array[[i]][1, 1] = raccoon_worm_array[1, i]

# }

# Initialize non-rodent worm array
for(i in 1:length(infra_nonmouse_worm_array)){
    infra_nonmouse_worm_array[[i]][1, 1] = 10

}

# Initialize non-rodent worm array
for(i in 1:length(infra_nonmouse_worm_array)){
    infra_mouse_worm_array[[i]][1, 1] = 0

}

all_worms_infra_array = list(infra_nonmouse_worm_array, infra_mouse_worm_array)

# Make total worm arrays
# infra_worm_array = combine_worm_arrays(all_worms_infra_array, raccoon_dead_alive_array)
raccoon_worm_array = get_tot_worm_array_from_infra(infra_nonmouse_worm_array, raccoon_dead_alive_array) 


init_arrays = list(raccoon_dead_alive_array=raccoon_dead_alive_array,
                   initial_age_vector=initial_age_vector,
                   age_array=age_array,
                   human_array=human_array,
                   human_risk_through_time=human_risk_through_time,
                   new_babies_vect=new_babies_vect,
                   infra_mouse_worm_array=infra_mouse_worm_array,
                   infra_nonmouse_worm_array=infra_nonmouse_worm_array,
                   all_worms_infra_array=all_worms_infra_array,
                   raccoon_worm_array=raccoon_worm_array)

return(init_arrays)