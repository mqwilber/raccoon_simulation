## Parameters for raccoon simulation ##

## Raccoon parameters
INIT_NUM_RACCOONS = 5
#DEATH_PROB = 22 # Natural death
DEATH_THRESHOLD = 22 #22 # beta in death_probability fxn
PATHOGENICITY = -4.2 # alpha in death_probability fxn

## Age parameters
FIRST_REPRO_AGE = 3 # months
LITTER_SIZE = 2 # Number of raccoons per litter
MONTH_AT_REPRO = 12 # Reproduce once a year

# PARASITE PARAMETERS
ENVIRONMENTAL_POOL = 0
EGG_PRODUCTION_PER_WORM = 100 # eggs per month (note this is very wrong)
ENCOUNTER_PROB = 1 # 100% probability that
ENCOUNTER_MEAN = 100 # Probability of a raccoon and egg coming into contact
INFECTIVITY = 0.1 # Probability of infectivity
RESISTANCE = 0.01 # How quickly a raccoon gains resistance based on previous load
EGG_LOSS_PROB =  .14 # Probability that an egg survives of 7 months

WORM_SURV_TRESH = 7 # Threshold parameter of worm survival probability
WORM_SURV_SLOPE = -1 # Slope of worm surv probability
# WORM_SURV_TRESH / WORM_SURV_SLOPE = 7 -> 0.5 prob of dying after 7 months

# Time parameters: Each time step is a month
TIME_STEPS = 10