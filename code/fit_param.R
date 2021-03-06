## Estimating raccoon age-survival curve
library(nlme)

# Fitting data from parameters document
age = c(0, 3, 10, 22, 34) # Months
death_prob_data = c(1 - (.52^(1/7))^4, 1 - .88^(1/3),
            1 - .92^(1/5), 1 - 0.72^(1 / 12), 1 -  0.84^(1 / 12))

y = death_prob_data  # Random death prob is 0.01
fit = lm(log(y) ~ age, weight=c(82, 24, 14, 18, 70)) # Weight by sample size!
plot(age, exp(predict(fit)), type='l', ylim=c(0, 0.6))
points(age, y)
summary(fit)

## estimating parasite age-survival curve
#data from Olsen 1958
weeks <- c(1,15,16,17,20,21,22,40,51)
wormpop_alive <- c(13,11,10,9,7,4,3,1,0)
wormpop_dead <- 13 - wormpop_alive
months <- weeks / 4

# Fit a Binomial GLM
fit1 = glm(cbind(wormpop_alive, wormpop_dead) ~ months, family=binomial)
summary(fit1)

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   4.7104     1.2854   3.665 0.000248 ***
# months       -0.9446     0.2693  -3.507 0.000452 ***

# Check how our model is doing against the data
plot(months, wormpop_alive / 13)
newdata=data.frame(months=min(months):max(months))
lines(newdata$months, predict(fit1, newdata=newdata, type="response"))
