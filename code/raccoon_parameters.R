## Parameters for raccoon simulation ##

## Raccoon parameters
INIT_NUM_RACCOONS = 7
#DEATH_PROB = 22 # Natural death
DEATH_THRESHOLD = 22 #22 # beta in death_probability fxn
PATHOGENICITY = -4.2 # alpha in death_probability fxn
BABY_DEATH = 1 - (.52^(1/7))^4 # Probability of dying as baby
INTRINSIC_DEATH_RATE = 0.33 # Age related death rate
RANDOM_DEATH_PROB = 0.01 # Lower bound to death prob


## Age parameters
FIRST_REPRO_AGE = 10 # months
LITTER_SIZE = 2 # Maximum litter size
MONTH_AT_REPRO = 12
DISPERSAL_AGE = 6 # months

## Ricker function for density-dependent recruitment of new babies
K_CAPACITY = 3 # "Carrying" capacity for raccoons. Need to figure out what
                # what is determining carrying capacity for the deterministic
                # version of this model
DI_NEW_BABY_DEATH_RATE = 0 # Density independent new baby death
BIRTH_RATE = log(LITTER_SIZE) # Gives birth rate. Little r in ricker function
BETA = BIRTH_RATE / K_CAPACITY


# PARASITE PARAMETERS
# ENVIRONMENTAL_POOL = 0
# EGG_PRODUCTION_PER_WORM = 100 # eggs per month (note this is very wrong)
ENCOUNTER_PROB = 1 # 100% probability that
ENCOUNTER_MEAN = 100 # Probability of a raccoon and egg coming into contact
ENCOUNTER_K = 1 # Aggregation parameter of the NBD
INFECTIVITY = 0.02 # Probability of infectivity
RESISTANCE = 0.1 # How quickly a raccoon gains resistance based on previous load

# See fit_param.R for how we got these values
WORM_SURV_TRESH = 4.7104 # Threshold parameter of worm survival probability
WORM_SURV_SLOPE = -0.9446 # Slope of worm surv probability

# Time parameters: Each time step is a month
TIME_STEPS = 23