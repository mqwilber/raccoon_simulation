## Estimating age-survival curve
library(nlme)

# Fitting data from parameters document
age = c(0, 3, 10, 22, 34) # Months
death_prob_data = c(1 - (.52^(1/7))^4, 1 - .88^(1/3),
            1 - .92^(1/5), 1 - 0.72^(1 / 12), 1 -  0.84^(1 / 12))

y = death_prob_data + 0.01  # Random death rate is 0.01
fit = lm(log(y) ~ log(age + 1))
plot(age, exp(predict(fit)), type='l', ylim=c(0, 0.6))
points(age, y)
summary(fit)

# Check this with nlme functions to see if we can get a better looking fit...
# Something other than exponential.

#the fit looks a lot better if the first point is turned into two (identical y value) points, 
#that seems really sketchy
#age = c(0,1, 3, 10, 22, 34) # Months
#death_prob_data = c(1 - (.52^(1/7))^4, .3118 , 1 - .88^(1/3),
           # 1 - .92^(1/5), 1 - 0.72^(1 / 12), 1 -  0.84^(1 / 12))
