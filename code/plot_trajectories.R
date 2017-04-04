## Plot different metric trajectories 
library(ggplot2)
library(gridExtra)
setwd("~/Repos/raccoon_simulation/code")
source("raccoon_fxns.R") # Load in the helper functions
set.seed(1)

# Set the possible management strategies
cull_params = list(strategy="random", quota=10, overlap_threshold=0.4)
birth_control_params = list(strategy="random", quota=8, overlap_threshold=0.7)
worm_control_params = list(strategy="random", quota=10000, overlap_threshold=0.8)

strats = list(cull=list(cull_params=cull_params, birth_control_params=NULL, worm_control_params=NULL),
     birth=list(cull_params=NULL, birth_control_params=birth_control_params, worm_control_params=NULL),
     worm=list(cull_params=NULL, birth_control_params=NULL, worm_control_params=worm_control_params))

# Loop through different strategies and make some example plots
for(stratnm in names(strats)){

    strat = strats[[stratnm]]

    # Run simulation
    management_time = 60
    time_steps = 200
    params = get_simulation_parameters(TIME_STEPS=time_steps)
    init_arrays = get_init_arrays(params) # Load in init arrays
    all_res = full_simulation(params, init_arrays, 
                              cull_params=strat[["cull_params"]], 
                              birth_control_params=strat[["birth_control_params"]],
                              worm_control_params=strat[["worm_control_params"]], 
                              management_time=management_time,
                              print_it=TRUE)

    # Make some plots
    pop_temp = rowSums(all_res$raccoon_dead_alive_array, na.rm=T)
    time = 1:length(pop_temp) / 12
    pop_dat = data.frame(list(population=pop_temp, time=time))
    pop_plot = ggplot(data=pop_dat, aes(x=time, y=population)) + geom_line() + 
                        xlab("Years") + ylab("Raccoon density") +
                        geom_vline(xintercept=management_time / 12, color="black", linetype="dashed") +
                        geom_rect(aes(xmin=time[length(time)] - 2, 
                                      xmax=time[length(time)], 
                                      ymax=600, ymin=0), fill="light gray", alpha=0.02) +
                        ylim(0, 600) +
                        theme_bw()

    worm_pop_traj = rowSums(all_res$raccoon_worm_array, na.rm=T)
    worm_dat = data.frame(worm=worm_pop_traj, time=time)
    worm_plot = ggplot(data=worm_dat, aes(x=time, y=log10(worm + 1))) + geom_line(color="#ED1C24") + 
                        xlab("Years") + ylab("log(worm population)") +
                        geom_vline(xintercept=management_time / 12, color="black", linetype="dashed") +
                        geom_rect(aes(xmin=time[length(time)] - 2, 
                                      xmax=time[length(time)], 
                                      ymax=3.75, ymin=0), fill="light gray", alpha=0.02) + 
                        ylim(0, 3.75) + 
                        theme_bw()
    worm_plot

    human_risk = get_human_risk_metric(all_res$eggproduction_array, params$EGG_DECAY)
    human_dat = data.frame(human=human_risk, time=time)
    human_plot = ggplot(data=human_dat, aes(x=time, y=human)) + geom_line(color="#27AAE1") + 
                        xlab("Years") + ylab("Human risk") +
                        geom_vline(xintercept=management_time / 12, color="black", linetype="dashed") + 
                        geom_rect(aes(xmin=time[length(time)] - 2, 
                                      xmax=time[length(time)], 
                                      ymax=800, ymin=0), fill="light gray", alpha=0.02) +
                        ylim(0, 800) +
                        theme_bw()
    human_plot

    prev = get_prevalence(all_res$raccoon_worm_array)
    prev_dat = data.frame(prev=prev, time=time)
    prev_plot = ggplot(data=prev_dat, aes(x=time, y=prev)) + geom_line(color="blue") + 
                    xlab("Years") + ylab("Prevalence") +
                    geom_vline(xintercept=management_time / 12, color="black", linetype="dashed") + 
                    geom_rect(aes(xmin=time[length(time)] - 2, 
                                  xmax=time[length(time)], 
                                  ymax=1, ymin=0), fill="light gray", alpha=0.02) +
                    ylim(0, 1) +
                    theme_bw()

    intensity = get_mean_intensity(all_res$raccoon_worm_array)
    intensity_dat = data.frame(intensity=intensity, time=time)
    intensity_plot = ggplot(data=prev_dat, aes(x=time, y=intensity)) + geom_line(color="purple") + 
                    xlab("Years") + ylab("Mean intensity") +
                    geom_vline(xintercept=management_time / 12, color="black", linetype="dashed") + 
                    geom_rect(aes(xmin=time[length(time)] - 2, 
                                  xmax=time[length(time)], 
                                  ymax=30, ymin=0), fill="light gray", alpha=0.02) +
                    ylim(0, 30) +
                    theme_bw()

    # Save the plots all together
    grid.arrange(pop_plot, worm_plot, human_plot, prev_plot, intensity_plot, nrow=3, ncol=2)
    g = arrangeGrob(pop_plot, worm_plot, human_plot, prev_plot, intensity_plot, nrow=3, ncol=2)
    ggsave(paste("../results/plots/example_trajectories_", stratnm,  ".jpg", sep=""), plot=g, width=10, height=6)
}


# Run simulation and make some plots
management_time = 1000
time_steps = 48
params = get_simulation_parameters(TIME_STEPS=time_steps)
init_arrays = get_init_arrays(params) # Load in init arrays
all_res = full_simulation(params, init_arrays, 
                          cull_params=NULL, 
                          birth_control_params=NULL,
                          worm_control_params=NULL, 
                          management_time=management_time,
                          print_it=TRUE)

pop_temp = rowSums(all_res$raccoon_dead_alive_array, na.rm=T)
worm_pop_traj = rowSums(all_res$raccoon_worm_array, na.rm=T)
mean_intens = get_mean_intensity(all_res$raccoon_worm_array)
prev = get_prevalence(all_res$raccoon_worm_array)
time = 1:length(pop_temp) / 12
pop_dat = data.frame(list(raccoons=log10(pop_temp), time=time,
                          worms=log10(worm_pop_traj),
                          mean_intensity=log10(mean_intens),
                          prevalance=prev))

melt_dat = melt(pop_dat, id.vars=c('time'))
traj_plot = ggplot(data=melt_dat, aes(x=time, y=value, color=variable)) + 
                geom_line() + xlab("Years") +
                ylim(0, 4) + ylab("log(raccoons), log(worms), log(mean_intensity), prevalence") +
                theme_bw() + scale_color_manual(values=c("black", "#ED1C24", "#56B4E9", "blue"))
traj_plot
ggsave("../results/plots/example_output.jpg", width=8, height=6)


# Age intensity plot

time_steps = 200
params = get_simulation_parameters(TIME_STEPS=time_steps)
init_arrays = get_init_arrays(params) # Load in init arrays
all_res = full_simulation(params, init_arrays, 
                          cull_params=NULL, 
                          birth_control_params=NULL,
                          worm_control_params=NULL, 
                          management_time=management_time,
                          print_it=TRUE)

range = 100:200
flat_age = as.vector(all_res$age_array[range, ])
flat_worms = as.vector(all_res$raccoon_worm_array[range, ])

full_dat = as.data.table(data.frame(list(age=flat_age, worms=flat_worms)))

prev_fxn = function(x){
    no_na_x = x[!is.na(x)]
    return(sum(no_na_x != 0) / length(no_na_x))
}
means = full_dat[, list(mean_worms=mean(worms, na.rm=T), 
                        se_worms=sd(worms) / sqrt(length(worms)),
                        prev=prev_fxn(worms)), 
                        by="age"]
means = means[!is.na(means$age), ]
means$age = means$age / 12
means[, c("lower", "upper") := list(mean_worms - se_worms, mean_worms + se_worms)]


age_plot = ggplot(means, aes(x=age, y=mean_worms)) +
          geom_line() + 
          geom_ribbon(aes(x=age, ymin=lower, ymax=upper), alpha=0.2) + 
          theme_bw() + ylim(0, 30) + xlab("Raccoon age (years)") +
          ylab("Mean intensity (+/- SE)")

age_plot

ggsave("../results/plots/mean_intensity.jpg", width=8, height=6)

prev_plot = ggplot(means, aes(x=age, y=prev)) +
          geom_line() + 
          theme_bw() + xlab("Raccoon age (years)") +
          ylab("Average prevalence")

prev_plot
ggsave("../results/plots/prevalence.jpg", width=8, height=6)


