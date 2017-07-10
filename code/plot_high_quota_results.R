# Plotting some of the simulation results
library(ggplot2)
library(data.table)
library(plyr)
library(RColorBrewer)


# Get results folders
folder_combos =
    list(human=c(Sys.glob("../results/birth_control_human/"),
                        Sys.glob("../results/cull_human/")),
        worm=Sys.glob("../results/worm_control_human/"),
        latrine=Sys.glob("../results/latrine_cleanup_human/"))

col_names = c("min_rac_pop", "mean_rac_pop", "max_rac_pop", 
             "min_worm_pop", "mean_worm_pop", "max_worm_pop",
             "min_human_risk", "mean_human_risk", "max_human_risk",
             "min_prev", "mean_prev", "max_prev",
             "min_intensity", "mean_intensity", "max_intensity")

high_effort = list()

all_data = list()

for(combo in names(folder_combos)){

    folder_names = folder_combos[[combo]]

    # In each result folder, extract all mean and var results
    mean_arrays = list()
    var_arrays = list()

    for(fname in folder_names){

        mean_files = Sys.glob(file.path(fname, "sim*mean*.rds"))
        fsplit = strsplit(fname, "/")[[1]][[3]]

        if(length(grep("worm_control", fsplit)) == 1) {
            quotas = c(0, 1, 10, 100, 1000, 10000, 100000)

        } else if((length(grep("latrine_cleanup", fsplit)) == 1)){
            quotas = c(1)
        } else{
            quotas = c(0, 10, 20, 50, 100, 500, 1000, 5000) # CHANGE
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
    vars = c("mean_rac_pop", "mean_worm_pop", "mean_human_risk", "mean_intensity", "mean_prev") 
             #   "mean_prev", "mean_intensity")
    plotting_df = list()

    for(varname in vars){

        tmelt_df = melt(sapply(mean_arrays, function(x) x[, varname]))
        tmeltvar_df = melt(sapply(var_arrays, function(x) x[, varname]))

        # Rename and adjust variables for easy ggplotting later

        if(combo != "latrine"){

            colnames(tmelt_df) = c("quota", "strategy", "mean_value")
            tmelt_df$sd_value = tmeltvar_df$value
            tmelt_df$quota = tmelt_df$quota
            tmelt_df$metric = varname
            plotting_df[[varname]] = tmelt_df
        } else{
            colnames(tmelt_df) = c("mean_value")
            tmelt_df$sd_value = tmeltvar_df$value
            tmelt_df$metric = varname
            tmelt_df$strategy = rownames(tmelt_df)
            plotting_df[[varname]] = tmelt_df

        }

    }

    # Plot the simulation results

    #quota_vals = log10(c(1, 10, 100, 1000, 10000)) # Use for worm control
    plotting_df = do.call(rbind, plotting_df)

    # Make extra column to make comparison easier
    specific_type = sapply(strsplit(as.character(plotting_df$strategy), "_"), function(x) x[length(x)])
    generic_type = sapply(strsplit(as.character(plotting_df$strategy), "_"), function(x) do.call(paste, as.list(x[1:(length(x) - 1)])))
    plotting_df$specific_type = specific_type
    plotting_df$generic_type = generic_type


    all_data[[combo]] = plotting_df

}

# Join the data and set it up for plotting

# Set the particular quotas to high
all_data[['human']]['high_low'] = ifelse(all_data[['human']]$quota == 5000, "High", "Low")
all_data[['worm']]['high_low'] = ifelse(all_data[['worm']]$quota == 100000, "High", "Low")
all_data[['latrine']]['quota'] = 1
all_data[['latrine']]['high_low'] = ifelse(all_data[['latrine']]$quota == 1, "High", "Low")

for(nm in names(all_data)){

    colnms = colnames(all_data[[nm]])
    all_data[[nm]] = all_data[[nm]][, colnms[order(colnms)]]
}

full_df = do.call(rbind, all_data)
full_dt  = as.data.table(full_df)

################# Tile plots ######################

tile_dt = full_dt[(specific_type == "human0.1.rds" | specific_type == "human0.5.rds" | specific_type == "human0.9.rds"), ]
tile_dt$mean_value[is.na(tile_dt$mean_value)] = 0

# Merge dataframe to ca
base_line_tile = full_dt[quota == 0 & generic_type == "cull human" & specific_type == "human0.1.rds", list(mean_value, metric)]
colnames(base_line_tile) = c('mean_value_baseline', 'metric')

setkey(tile_dt, metric)
setkey(base_line_tile, metric)

tile_dt2 = merge(tile_dt, base_line_tile, all.x=T)
tile_dt2$management_effect = tile_dt2$mean_value / tile_dt2$mean_value_baseline

tile_dt2$specific_type = revalue(tile_dt2$specific_type, 
                            c("human0.1.rds"="overlap > 0.1",
                              "human0.5.rds"="overlap > 0.5",
                              "human0.9.rds"="overlap > 0.9"))

tile_dt2$generic_type = revalue(tile_dt2$generic_type,
                               c("birth control human"="Birth control",
                                 "cull human"="Cull",
                                 "latrine cleanup human"="Latrine cleanup",
                                 "worm control human"="Worm control"))

tile_dt2$generic_type = ordered(tile_dt2$generic_type, c("Birth control", "Cull", "Worm control", "Latrine cleanup"))

tile_dt2$metric = revalue(tile_dt2$metric,
                        c("mean_human_risk" = "Human risk",
                          "mean_intensity" = "Mean infection intensity",
                          "mean_prev" = "Mean prevalence",
                          "mean_rac_pop"="Raccoon population",
                          "mean_worm_pop"="Worm population"))

# Difference from baseline

myPalette <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")

tile_plot = ggplot(as.data.frame(tile_dt2), aes(x=as.factor(quota), y=specific_type, fill=management_effect)) +
      geom_tile() + 
      scale_fill_gradientn(colours = myPalette(8), name="Management effect", limits=c(0, 2)) +
      facet_grid(metric~generic_type, scales="free") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("Human overlap") + xlab("Management effort")

ggsave("../results/plots/tile_plot_of_management_results.pdf", width=12, height=8)

################# Bar plots ######################

base_line = full_dt[quota == 0 & generic_type == "cull human" & specific_type == "human0.1.rds", ]

# Only include the relevant human overlaps
high_dt = full_dt[(high_low == "High") & (specific_type == "human0.1.rds" | specific_type == "human0.5.rds" | specific_type == "human0.9.rds"), ]

# Replace values for prettier plotting
high_dt$specific_type = revalue(high_dt$specific_type, 
                            c("human0.1.rds"="overlap > 0.1",
                              "human0.5.rds"="overlap > 0.5",
                              "human0.9.rds"="overlap > 0.9"))

high_dt$generic_type = revalue(high_dt$generic_type,
                               c("birth control human"="birth control",
                                 "cull human"="cull",
                                 "latrine cleanup human"="latrine cleanup",
                                 "worm control human"="worm control"))

high_dt$metric = revalue(high_dt$metric,
                        c("mean_human_risk" = "Human risk",
                          "mean_intensity" = "Mean infection intensity",
                          "mean_prev" = "Mean prevalence",
                          "mean_rac_pop"="Raccoon population",
                          "mean_worm_pop"="Worm population"))

base_line$metric = revalue(base_line$metric,
                        c("mean_human_risk" = "Human risk",
                          "mean_intensity" = "Mean infection intensity",
                          "mean_prev" = "Mean prevalence",
                          "mean_rac_pop"="Raccoon population",
                          "mean_worm_pop"="Worm population"))

high_dt$mean_value[is.na(high_dt$mean_value)] = 0
high_dt$sd_value[is.na(high_dt$sd_value)] = 0

high_plot = ggplot(as.data.frame(high_dt), aes(x=generic_type, y=mean_value, fill=specific_type)) + 
            geom_bar(position=position_dodge(), stat="identity") +
            geom_errorbar(aes(ymin=mean_value - sd_value, ymax=mean_value + sd_value),
                  width=.3, 
                  size=0.3,                   # Width of the error bars
                  position=position_dodge(0.9)) +
            scale_fill_manual(name="Human overlap", values=c("#deebf7", "#9ecae1", "#3182bd")) +
            geom_hline(data=as.data.frame(base_line), aes(yintercept=mean_value), linetype="dashed") +
            facet_wrap(~ metric, scales="free") + 
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                        legend.position = c(0.8, 0.2)) +
            xlab("Management strategy") + ylab("Mean value +/- SE") 
ggsave("../results/plots/high_effort_management_results.pdf", width=10, height=7.5)