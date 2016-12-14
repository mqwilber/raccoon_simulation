# Plotting some of the simulation results
library(ggplot2)
library(data.table)

# Get results folders
folder_names = Sys.glob("../results/*")

# In each result folder, extract all mean and var results
mean_arrays = list()
var_arrays = list()
for(fname in folder_names){

    mean_files = Sys.glob(file.path(fname, "sim*mean*.rds"))

    for(mfile in mean_files){
        mean_dat = readRDS(mfile)
        var_dat = readRDS(gsub("mean", "var", mfile))
        splt = strsplit(mfile, "_")[[1]]
        mean_arrays[[splt[length(splt)]]] = do.call(rbind, mean_dat)
        var_arrays[[splt[length(splt)]]] = do.call(rbind, var_dat)

    }

}

# Extract and format means for plotting
vars = c("mean_rac_pop", "mean_worm_pop", "mean_human_risk")
plotting_df = list()

for(varname in vars){

    tmelt_df = melt(sapply(mean_arrays, function(x) x[, varname]))
    tmeltvar_df = melt(sapply(var_arrays, function(x) x[, varname]))

    # Rename and adjust variables for easy ggplotting later
    colnames(tmelt_df) = c("quota", "strategy", "mean_value")
    tmelt_df$sd_value = tmeltvar_df$value
    tmelt_df$quota = tmelt_df$quota - 1
    tmelt_df$metric = varname
    plotting_df[[varname]] = tmelt_df

}

# Plot the simulation results

plotting_df = do.call(rbind, plotting_df)
gp = ggplot(plotting_df, aes(x=quota, y=log10(mean_value), color=strategy)) + 
                    geom_line() +
                    geom_ribbon(aes(ymin=log10(mean_value - sd_value), 
                                    ymax=log10(mean_value + sd_value),
                                    color=strategy), alpha=0.3) +
                    facet_wrap(~metric) + theme_bw() + xlab("Quota") + 
                    ylab("log(population) or log(metric)") +
                    scale_x_continuous(breaks=0:10)


gp