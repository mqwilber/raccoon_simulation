## Parameters for raccoon simulation ##

## Raccoon parameters
INIT_NUM_RACCOONS = 5
DEATH_PROB = 0
DEATH_THRESHOLD = 300 # Worms at which Raccoon is dead
PATHOGENICITY = (1 - DEATH_PROB) / DEATH_THRESHOLD

## Age parameters
FIRST_REPRO_AGE = 3 # months
LITTER_SIZE = 2 # Number of raccoons per litter
MONTH_AT_REPRO = 4 # Reproduce once a year

# PARASITE PARAMETERS
ENVIRONMENTAL_POOL = 0
EGG_PRODUCTION_PER_WORM = 100 # eggs per month (note this is very wrong)
ENCOUNTER_PROB = 1 # 100% probability that
ENCOUNTER_MEAN = 100 # Probability of a raccoon and egg coming into contact
INFECTIVITY = 0.5 # Probability of infectivity
EGG_LOSS_PROB =  .14 # Probability that an egg survives of 7 months

# Time parameters: Each time step is a month
TIME_STEPS = 8