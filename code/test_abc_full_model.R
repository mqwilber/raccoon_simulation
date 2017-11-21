## Script simulates a dataset from raccoon IBM and then uses the results
## to test whether the ABC method can recover the known parameter values
## from different attributes of the age-intensity/prevalence profile
source("abc_sim.R")
set.seed(5)


############# BUILD THE DATASET TO USE ########################

# Set default model parameters that we will try to estimate  
# AGE_EGG_RESISTANCE is fixed at 4. 
model_params = list(RODENT_ENCOUNTER_PROB=0.8,
                        ENCOUNTER_K=0.6,
                        EGG_CONTACT=5e-6, 
                        AGE_SUSCEPTIBILITY=1000, 
                        CLEAR_PROB=0, 
                        AGE_EGG_RESISTANCE=4)

run_sim = TRUE # IF you want to run the ABC simulation
if(run_sim){

    # sim_data = build_simulated_datasets(model_params, num_sets=1, TIME_STEPS=100,
    #                 datasource="../data/formatted/raccoon_age_intensity_full.csv")[[1]]

    # # Convert sim data to correct format and save
    # sim_data = as.data.frame(sim_data)
    # colnames(sim_data) = c("sim_num", "abundance", "prev", "iqr", "age_lower")
    # sim_data$prev = sim_data$prev * 100
    # sim_data$age_upper = c(3, 4, 5, 7, 12, 26, 48, 77)
    # sim_data$sample_size = c(17, 26, 37, 27, 22, 26, 21, 13)
    # write.csv(sim_data, "sim_data_test.csv", row.names=FALSE)

    ###############################################################

    ############# LOOP THROUGH DIFFERENT ABC STRATEGIES ###########

    methods = c("euclidean") #c("euclidean", "manhattan")
    stat_sets = c("all") #c("means", "prevs", "means_prevs", "all")
     

    for(method in methods){
        for(stat_set in stat_sets){

            print(paste(method, stat_set))

            # Use the ABC method to estimate the parameters of the model given the simulated
            # data
            steps = 3
            particles = 16
            percent_rj = 0.5
            cores = 4
            abc_fit = run_abc(steps, particles, models=c(1, 2), method=method, 
                                stat_set=stat_set, 
                                datasource="sim_data_test.csv", 
                                percent_rj=percent_rj, 
                                cores=cores)

            # Save results
            sim_results = list(abc_fit, model_params)
            # saveRDS(sim_results, paste("test_abc_results_multimodel", "_", method, "_", 
            #                                          stat_set, ".rds", sep=""))
        }
    }

}

###############################################################

# library(ggplot2)
# library(data.table)
# library(gridExtra)

# res_files = Sys.glob("sim_results/test_abc_results_*lowtol*.rds")

# prior_params = list(RODENT_ENCOUNTER_PROB=c(lower=0, upper=1),
#                 ENCOUNTER_K=c(lower=0, upper=1),
#                 EGG_CONTACT=c(lower=0.001, upper=20), 
#                 ENCOUNTER_MEAN=c(lower=10, upper=1000),
#                 AGE_SUSCEPTIBILITY=c(lower=0.001, upper=5),
#                 CLEAR_PROB=c(lower=0, upper=1))

# for(res_file in res_files){

#     tabc_fit = readRDS(res_file)[[1]][[1]]
#     # Examine the posterior distributions after ABC
#     plots = list()

#     for(i in 1:5){

#         tdat = melt(do.call(rbind, tabc_fit[[i]]))
#         colnames(tdat) = c("index", "param", "value")

#         tdat$true_val = -1
#         tdat$prior_samp = -1
#         for(nm in names(prior_params)){
#             tdat$true_val[tdat$param == nm] = model_params[[nm]]
#             tdat$prior_samp[tdat$param == nm] = runif(sum(tdat$param == nm), 
#                                             min=prior_params[[nm]]['lower'],
#                                             max=prior_params[[nm]]['upper'])
#         } 


#         plots[[i]] = ggplotGrob(ggplot(data=tdat, aes(x=value)) + geom_density(fill="gray", alpha=0.5) + 
#                         geom_density(aes(x=prior_samp), fill="blue", alpha=0.1) +
#                         geom_vline(aes(xintercept=true_val), color="red") +
#                         facet_wrap(~param, scales="free") + theme_bw() +
#                         xlab("Parameter value") + ylab("Density") +
#                         ggtitle(paste("Iteration", i)))
#     }

#     (paste("../results/plots/", strsplit(strsplit(res_file, "[.]")[[1]][1], "[/]")[[1]][2], ".pdf", sep=""))
#     do.call(grid.arrange, plots)
#     dev.off()

#     ggsave(paste("../results/plots/", strsplit(strsplit(res_file, "[.]")[[1]][1], "[/]")[[1]][2], "_iteration5.pdf", sep=""),
#                 plots[[5]])

#     # Plot the model predictions against model data

#     rand_samp = 30
#     rand_draws = list()
#     for(samp in 1:rand_samp){
#         print(paste("Sample", samp))

#         tmodel_params = as.list(tabc_fit[[steps]][sample(1:length(tabc_fit[[steps]]), 1)])

#         draw_results = build_simulated_datasets(tmodel_params, num_sets=1, 
#                                 TIME_STEPS=100,
#                                 datasource="sim_data_full_model.csv")[[1]]
#         draw_results[, 'sim_num'] = samp
#         rand_draws[[samp]] = draw_results
#     }

#     obs = read.csv("sim_data_full_model.csv")
#     draw_resultsdf = data.frame(do.call(rbind, rand_draws))
#     colnames(draw_resultsdf) = c("sim_num", "abundance", "prev", "iqr", "age")
#     dt = as.data.table(draw_resultsdf)
#     summ_dt = dt[, list(abundance=quantile(abundance, 0.5), 
#                         abund_lower=quantile(abundance, 0.025),
#                         abund_upper=quantile(abundance, 0.975),
#                         prev=quantile(prev, 0.5), 
#                         prev_lower=quantile(prev, 0.025),
#                         prev_upper=quantile(prev, 0.975)), by=c("age")]
    
#     # Compare observed ad predicted
#     abund_plot = ggplot() + geom_point(data=summ_dt, aes(x=age, y=abundance))  +
#                     geom_errorbar(data=summ_dt, aes(x=age, ymin=abund_lower, ymax=abund_upper)) +
#                     geom_point(data=obs, aes(x=age_lower + 0.5, y=abundance, color="red")) 


#     prev_plot = ggplot(summ_dt, aes(x=age, y=prev)) + geom_point() + 
#                 geom_errorbar(aes(x=age, ymin=prev_lower, ymax=prev_upper)) +
#                 geom_point(data=obs, aes(x=age_lower, y=prev/100, color="red")) + 
#                 #geom_errorbar(data=obs, aes(x=age_lower, ymin=prev_lower/100, ymax=prev_upper/100, color="red")) +
#                 ylim(0, 1)
# }

# # ## Using the the 

