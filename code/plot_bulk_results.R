# Plotting some of the simulation results
library(ggplot2)
library(data.table)

# Get results folders
folder_combos = list(human=Sys.glob("../results/latrine_cleanup_human/"))

# list(human=c(Sys.glob("../results/birth_control_human/"),
#                         Sys.glob("../results/cull_human/")),
#                      random=c(Sys.glob("../results/birth_control_random/"),
#                         Sys.glob("../results/cull_random/")),
#                      worm=Sys.glob("../results/worm_control_*/"),
#                      age=c(Sys.glob("../results/cull_age/"), 
#                            Sys.glob("../results/cull_random/")))
                
                     
                 # Sys.glob("../results/cull_random/"),
                 # Sys.glob("../results/birth_control_random/"))
                # Sys.glob("../results/birth_control_random/"))

col_names = c("min_rac_pop", "mean_rac_pop", "max_rac_pop", 
             "min_worm_pop", "mean_worm_pop", "max_worm_pop",
             "min_human_risk", "mean_human_risk", "max_human_risk",
             "min_prev", "mean_prev", "max_prev",
             "min_intensity", "mean_intensity", "max_intensity")



#c(Sys.glob("../results/birth_control_human/"), 
               #  Sys.glob("../results/cull_random/"))

for(combo in names(folder_combos)){

    folder_names = folder_combos[[combo]]

    # In each result folder, extract all mean and var results
    mean_arrays = list()
    var_arrays = list()

    for(fname in folder_names){

        mean_files = Sys.glob(file.path(fname, "sim*mean*.rds"))
        fsplit = strsplit(fname, "/")[[1]][[3]]

        if(length(grep("worm_control", fsplit)) == 1) {
            quotas = c(0, log10(c(1, 10, 100, 1000, 10000, 100000)) + 1)
        } else if((length(grep("latrine_cleanup", fsplit)) == 1)){
            quotas = c(1)
        } else{
            quotas = log10(c(0:10, 20, 50, 100, 200)) # CHANGE
        }

        for(mfile in mean_files){

            mean_dat = readRDS(mfile)
            var_dat = readRDS(gsub("mean", "var", mfile))
            splt = strsplit(mfile, "_")[[1]]

            mean_arrays[[paste(fsplit, splt[length(splt)], sep="_")]] = do.call(rbind, mean_dat)
            var_arrays[[paste(fsplit, splt[length(splt)], sep="_")]] = do.call(rbind, var_dat)

            rownames(mean_arrays[[paste(fsplit, splt[length(splt)], sep="_")]]) = quotas
            rownames(var_arrays[[paste(fsplit, splt[length(splt)], sep="_")]]) = quotas

            colnames(mean_arrays[[paste(fsplit, splt[length(splt)], sep="_")]]) = col_names
            colnames(var_arrays[[paste(fsplit, splt[length(splt)], sep="_")]]) = col_names

            # mean_arrays[[paste(fsplit, splt[length(splt)], sep="_")]][, "mean_intensity"] = 
            #             mean_arrays[[paste(fsplit, splt[length(splt)], sep="_")]][, 'mean_worm_pop'] / mean_arrays[[paste(fsplit, splt[length(splt)], sep="_")]][, 'mean_rac_pop']
        }

    }

    # Extract and format means for plotting
    vars = c("mean_rac_pop", "mean_worm_pop", "mean_human_risk", "mean_intensity") 
             #   "mean_prev", "mean_intensity")
    plotting_df = list()

    for(varname in vars){

        tmelt_df = melt(sapply(mean_arrays, function(x) x[, varname]))
        tmeltvar_df = melt(sapply(var_arrays, function(x) x[, varname]))

        # Rename and adjust variables for easy ggplotting later
        colnames(tmelt_df) = c("quota", "strategy", "mean_value")
        tmelt_df$sd_value = tmeltvar_df$value
        tmelt_df$quota = tmelt_df$quota
        tmelt_df$metric = varname
        plotting_df[[varname]] = tmelt_df

    }

    # Plot the simulation results

    #quota_vals = log10(c(1, 10, 100, 1000, 10000)) # Use for worm control
    plotting_df = do.call(rbind, plotting_df)

    # Make extra column to make comparison easier
    specific_type = sapply(strsplit(as.character(plotting_df$strategy), "_"), function(x) x[length(x)])
    generic_type = sapply(strsplit(as.character(plotting_df$strategy), "_"), function(x) do.call(paste, as.list(x[1:(length(x) - 1)])))
    plotting_df$specific_type = specific_type
    plotting_df$generic_type = generic_type


    #plotting_df = plotting_df[plotting_df$metric != "mean_worm_pop", ]
    gp = ggplot(plotting_df, aes(x=quota, y=log10(mean_value + 1), linetype=generic_type, color=specific_type)) + 
                        geom_line() +
                        # geom_ribbon(aes(ymin=log10(mean_value + 1 - sd_value), 
                        #                 ymax=log10(mean_value  + 1 + sd_value),
                        #                 color=strategy), linetype=0, alpha=0.3) +
                        facet_wrap(~metric) + theme_bw() + xlab("log(quota)") + ylim(0,3.5) +
                        ylab("log(population) or log(metric)")
                        #scale_x_continuous(breaks=quota_vals)


    gp2 = ggplot(plotting_df, aes(x=quota, y=mean_value, linetype=generic_type, color=specific_type)) + 
                        geom_line() +
                        # geom_ribbon(aes(ymin=mean_value - sd_value, 
                        #                 ymax=mean_value + sd_value,
                        #                 color=strategy), alpha=0.3) +
                        facet_wrap(~metric, scales="free") + theme_bw() + xlab("log(quota)") + 
                        ylab("Population or metric") #+ ylim(c(0, 550))# +

                        #scale_x_continuous(breaks=quota_vals)

    ggsave(paste("../results/plots/bulk_results_", combo, ".jpg", sep=""), plot=gp, width=8, height=6)

}