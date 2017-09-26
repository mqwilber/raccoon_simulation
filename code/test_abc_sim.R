## Script simulates a dataset from raccoon IBM and then uses the results
## to test whether the ABC method can recover the known parameter values.

source("abc_sim.R")

# Set default model parameters that we will try to estimate
# AGE_SUSCEPTIBILITY and CLEAR_PROB will not be estimated so they are set to 
# default values
model_params = list(RODENT_ENCOUNTER_PROB=0.8,
                        ENCOUNTER_K=0.6,
                        EGG_CONTACT=3.696582, 
                        ENCOUNTER_MEAN=459.174, 
                        AGE_SUSCEPTIBILITY=10000, 
                        CLEAR_PROB=0, 
                        AGE_EGG_RESISTANCE=4)

sim_data = build_simulated_datasets(model_params, num_sets=1, TIME_STEPS=100)[[1]]

# Convert sim data to correct format and save
sim_data = as.data.frame(sim_data)
colnames(sim_data) = c("sim_num", "abundance", "prev", "iqr", "age_lower")
sim_data$prev = sim_data$prev * 100
sim_data$age_upper = c(3, 4, 5, 7, 12, 26, 48, 77)
sim_data$sample_size = c(17, 26, 37, 27, 22, 26, 21, 13)
write.csv(sim_data, "sim_data.csv", row.names=FALSE)

# Use the ABC method to estimate the parameters of the model given the simulated
# data
abc_fit = run_abc(4, 48, datasource="sim_data.csv", percent_rj=.40, cores=4)

# Save results

sim_results = list(abc_fit, model_params)
saveRDS(sim_results, "test_abc_results.rds")
