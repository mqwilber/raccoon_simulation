## Read in Sara's raccoon data and extract the IQR for each age class
library(data.table)

raw_dat = as.data.table(read.csv("../archival/raccoon_worms.csv"))
intens_dat = read.csv("../archival/raccoon_age_intensity.csv")

age_classes = unique(raw_dat$Age_class)
age_means = array(NA, dim=length(age_classes))
age_iqrs = array(NA, dim=length(age_classes))

# Get the IQR range
iqr_dat = raw_dat[, list(mean_worms=mean(ad_worms), 
              iqr=diff(quantile(ad_worms, c(0.25, 0.75)))),
               by=c("Age_class")]

intens_dat$iqr = iqr_dat$iqr

write.csv(intens_dat, "raccoon_age_intensity_full.csv", row.names=FALSE)

