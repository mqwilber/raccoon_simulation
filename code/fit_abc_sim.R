## Script fits the raccoon IBM to the data from Weinstein (2016) using ABC-SMC
source("abc_sim.R")


# Use the ABC method to estimate the parameters of the model given the simulated
# data
steps = 5
particles = 10000
percent_rj = 0.05
cores = 10
abc_fit = run_abc(steps, particles, datasource="../data/formatted/raccoon_age_intensity_full.csv", 
                    percent_rj=percent_rj, cores=cores)

# Save results
saveRDS(abc_fit, "fit_abc_results.rds")
abc_fit = readRDS("fit_abc_results.rds")

knot = FALSE

if(!knot){
    # Make some plots of the ABC results
    library(ggplot2)
    library(data.table)
    library(gridExtra)

    # Examine the posterior distributions after ABC
    plots = list()

    for(i in 1:steps){

        tdat = melt(abc_fit[[1]][[i]])
        colnames(tdat) = c("index", "param", "value")
        tdat$param = plyr::revalue(tdat$param, replace=c("ex"="encounter agg.", "em"="encounter mean",
                                      "ec"="egg contact decay", 
                                      "rp"="rodent encounter prob."))

        plots[[i]] = ggplot(data=tdat, aes(x=value)) + geom_histogram() + 
                        facet_wrap(~param, scales="free") +
                        ggtitle(paste("Iteration", i))
    }


    gridplot = do.call(grid.arrange, plots)

    # Plot the model predictions against model data

    rand_samp = 100
    rand_draws = list()
    for(samp in 1:rand_samp){
        print(paste("Sample", samp))

        param_samp = abc_fit[[1]][[steps]][sample(1:nrow(abc_fit[[1]][[steps]]), 1), ]
        tmodel_params = list(RODENT_ENCOUNTER_PROB=param_samp['rp'], 
                            ENCOUNTER_K=(1 - param_samp['ex']) / param_samp['ex'], 
                            EGG_CONTACT=param_samp['ec'], 
                            ENCOUNTER_MEAN=param_samp['em'])

        draw_results = build_simulated_datasets(tmodel_params, num_sets=1, 
                                TIME_STEPS=100,
                                datasource="../data/formatted/raccoon_age_intensity_full.csv")[[1]]
        draw_results[, 'sim_num'] = samp
        rand_draws[[samp]] = draw_results
    }

    obs = read.csv("../data/formatted/raccoon_age_intensity_full.csv")
    draw_resultsdf = data.frame(do.call(rbind, rand_draws))
    colnames(draw_resultsdf) = c("sim_num", "abundance", "prev", "iqr", "age")
    dt = as.data.table(draw_resultsdf)
    summ_dt = dt[, list(abundance=quantile(abundance, 0.5), 
                        abund_lower=quantile(abundance, 0.025),
                        abund_upper=quantile(abundance, 0.975),
                        prev=quantile(prev, 0.5), 
                        prev_lower=quantile(prev, 0.025),
                        prev_upper=quantile(prev, 0.975)), by=c("age")]
    
    # Compare observed ad predicted
    abund_plot = ggplot(data=summ_dt, aes(x=age, y=abundance)) + geom_point()  +
                    geom_errorbar(aes(x=age, ymin=abund_lower, ymax=abund_upper)) +
                    geom_point(data=obs, aes(x=age_lower + 0.5, y=abundance, color="red")) + 
                    geom_errorbar(data=obs, aes(x=age_lower + 0.5, ymin=abund_lower, ymax=abund_upper, color="red"))


    prev_plot = ggplot(summ_dt, aes(x=age, y=prev)) + geom_point() + 
                geom_errorbar(aes(x=age, ymin=prev_lower, ymax=prev_upper)) +
                geom_point(data=obs, aes(x=age_lower, y=prev/100, color="red")) + 
                geom_errorbar(data=obs, aes(x=age_lower, ymin=prev_lower/100, ymax=prev_upper/100, color="red")) +
                ylim(0, 1)


}
