# Loading in the trapping data
dat = read.csv("trapping.csv")

# TODO: SHOW PLOT

# Don't include the really good trappers
dat_sub = dat[dat$gehrt == "N", ]

# Convert the trapping prob. to a log scale so we can use
# a linear model. Type I functional response
prob_trap = dat_sub$capture.trapnight
prob_no_trap = log(1 - prob_trap)

# Divide density by 2 to only get females
density = dat_sub$rac.km / 2

# Fit the linear model
fit = lm(prob_no_trap ~ density - 1)

# Plot the results of the model
plot(density, prob_no_trap, ylim=c(-.3, 0))
lines(density, predict(fit), type="l")

summary(fit)

# Plot the results
plot(density, prob_trap, ylim=c(0, .4))
points(density, 1 - exp(predict(fit)), col="red")
