# Plotting some of the simulation results
library(ggplot2)
library(data.table)

# Get results folders
folder_names = Sys.glob("../results/cull_age/")

#c(Sys.glob("../results/birth_control_human/"), 
               #  Sys.glob("../results/cull_random/"))

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

quota_vals = 0:10
#quota_vals = log10(c(1, 10, 100, 1000, 10000)) # Use for worm control
plotting_df = do.call(rbind, plotting_df)
plotting_df = plotting_df[plotting_df$metric != "mean_worm_pop", ]
gp = ggplot(plotting_df, aes(x=quota, y=log10(mean_value + 1), color=strategy)) + 
                    geom_line() +
                    # geom_ribbon(aes(ymin=log10(mean_value + 1 - sd_value), 
                    #                 ymax=log10(mean_value  + 1 + sd_value),
                    #                 color=strategy), alpha=0.3) +
                    facet_wrap(~metric) + theme_bw() + xlab("Quota") + 
                    ylab("log(population) or log(metric)") + ylim(c(0, 3))# +
                    #scale_x_continuous(breaks=quota_vals)

gp2 = ggplot(plotting_df, aes(x=quota, y=mean_value, color=strategy)) + 
                    geom_line() +
                    # geom_ribbon(aes(ymin=mean_value - sd_value, 
                    #                 ymax=mean_value + sd_value,
                    #                 color=strategy), alpha=0.3) +
                    facet_wrap(~metric) + theme_bw() + xlab("Quota") + 
                    ylab("Population or metric") + ylim(c(0, 550))# +
                    #scale_x_continuous(breaks=quota_vals)

