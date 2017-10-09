## Script simulates a dataset from raccoon IBM and then uses the results
## to test whether the ABC method can recover the known parameter values
## from different attributes of the age-intensity/prevalence profile
source("abc_sim.R")
set.seed(5)


############# BUILD THE DATASET TO USE ########################

# Set default model parameters that we will try to estimate
# AGE_SUSCEPTIBILITY and CLEAR_PROB will not be estimated so they are set to 
# default values.  AGE_EGG_RESISTANCE is fixed at 4. 
model_params = list(RODENT_ENCOUNTER_PROB=0.8,
                        ENCOUNTER_K=0.6,
                        EGG_CONTACT=3.696582, 
                        ENCOUNTER_MEAN=459.174, 
                        AGE_SUSCEPTIBILITY=10000, 
                        CLEAR_PROB=0, 
                        AGE_EGG_RESISTANCE=4)

sim_data = build_simulated_datasets(model_params, num_sets=1, TIME_STEPS=100,
                datasource="../data/formatted/raccoon_age_intensity_full.csv")[[1]]

# Convert sim data to correct format and save
sim_data = as.data.frame(sim_data)
colnames(sim_data) = c("sim_num", "abundance", "prev", "iqr", "age_lower")
sim_data$prev = sim_data$prev * 100
sim_data$age_upper = c(3, 4, 5, 7, 12, 26, 48, 77)
sim_data$sample_size = c(17, 26, 37, 27, 22, 26, 21, 13)
write.csv(sim_data, "sim_data.csv", row.names=FALSE)

###############################################################

############# LOOP THROUGH DIFFERENT ABC STRATEGIES ###########

methods = c("euclidean", "manhattan")
stat_sets = c("means", "prevs", "means_prevs", "all")
 

for(method in methods){
    for(stat_set in stat_sets){

        print(paste(method, stat_set))

        # Use the ABC method to estimate the parameters of the model given the simulated
        # data
        steps = 5
        particles = 10000
        percent_rj = 0.05
        cores = 10
        abc_fit = run_abc(4, 24, method=method, stat_set=stat_set, 
                            datasource="sim_data.csv", percent_rj=.50, cores=4)

        # Save results
        sim_results = list(abc_fit, model_params)
        saveRDS(sim_results, paste("test_abc_results", "_", method, "_", 
                                                stat_set, ".rds", sep=""))
    }
}

###############################################################


# Make some plots of the ABC results
# library(ggplot2)
# library(data.table)
# library(gridExtra)

# # Examine the posterior distributions after ABC
# plots = list()

# for(i in 1:5){

#     tdat = melt(abc_fit[[i]])
#     colnames(tdat) = c("index", "param", "value")
#     tdat$param = plyr::revalue(tdat$param, replace=c("ex"="encounter agg.", "em"="encounter mean",
#                                   "ec"="egg contact decay", 
#                                   "rp"="rodent encounter prob."))

#     plots[[i]] = ggplot(data=tdat, aes(x=value)) + geom_histogram() + 
#                     facet_wrap(~param, scales="free") +
#                     ggtitle(paste("Iteration", i))
# }


# do.call(grid.arrange, plots)
