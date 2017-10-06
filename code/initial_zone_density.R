# Calculate mean initial density of model

ZONES = 10
sims = 50
init_dens = array(NA, dim=c(sims, ZONES))

for(i in 1:sims){
    print(paste("Starting sim", i))
    source("management_strategy.R")
    init_dens[i, ] = get_zone_density(all_res$human_overlap_through_time[[1]], 10)
}

mean_init = apply(init_dens, 2, mean)