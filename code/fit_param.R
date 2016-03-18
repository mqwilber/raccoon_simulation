## Estimating age-survival curve
library(nlme)

# Fitting data from parameters document
age = c(0, 3, 10, 22, 34) # Months
death_prob_data = c(1 - (.52^(1/7))^4, 1 - .88^(1/3),
            1 - .92^(1/5), 1 - 0.72^(1 / 12), 1 -  0.84^(1 / 12))

y = death_prob_data + 0.1  # Random death rate is 0.1
fit = lm(log(y) ~ log(age + 1))
plot(age, exp(predict(fit)), type='l', ylim=c(0, 0.6))
points(age, y)
summary(fit)

# Check this with nlme functions to see if we can get a better looking fit...
# Something other than exponential.