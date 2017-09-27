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

knot = TRUE

if(!knot){
    # Make some plots of the ABC results
    library(ggplot2)
    library(data.table)
    library(gridExtra)

    # Examine the posterior distributions after ABC
    plots = list()

    for(i in 1:4){

        tdat = melt(abc_fit[[1]][[i]])
        colnames(tdat) = c("index", "param", "value")
        tdat$param = plyr::revalue(tdat$param, replace=c("ex"="encounter agg.", "em"="encounter mean",
                                      "ec"="egg contact decay", 
                                      "rp"="rodent encounter prob."))

        plots[[i]] = ggplot(data=tdat, aes(x=value)) + geom_histogram() + 
                        facet_wrap(~param, scales="free") +
                        ggtitle(paste("Iteration", i))
    }


    do.call(grid.arrange, plots)
}
