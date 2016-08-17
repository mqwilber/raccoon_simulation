## Parameters for raccoon simulation ##

## Raccoon parameters
INIT_NUM_RACCOONS = 2000
#DEATH_PROB = 22 # Natural death
DEATH_THRESHOLD = 22 #22 # beta in death_probability fxn
PATHOGENICITY = -4.2 # alpha in death_probability fxn
BABY_DEATH = 1 - (.52^(1/7))^4 # Probability of dying as baby
INTRINSIC_DEATH_RATE = 0.33 # Age related death rate
RANDOM_DEATH_PROB = 0.01 # Lower bound to death prob
OLD_DEATH = (1 / (20 * 12)^2) # Above 20 years old the raccoon dies
AGE_EGG_RESISTANCE = 4 # Age above which raccoons no longer pick up eggs
RODENT_ENCOUNTER_PROB = 0.5 # Monthly probability of encountering a rodent

## Rodent parameters
MOUSE_WORM_MEAN = 3.49 # Abundance of worms in peromyscus estimated from Sara's data
MOUSE_WORM_AGG = 0.22 # Aggregation of worms (k parameter) in peromyscus
LARVAL_WORM_INFECTIVITY = 0.25 # Probability of larval worm establishing


## Age parameters
FIRST_REPRO_AGE = 10 # months
LITTER_SIZE = 4 # Maximum litter size
MONTH_AT_REPRO = 12
DISPERSAL_AGE = 12 # months

## Ricker function for density-dependent recruitment of new babies
K_CAPACITY = 400 # "Carrying" capacity for raccoons. Need to figure out what
                # what is determining carrying capacity for the deterministic
                # version of this model
DI_NEW_BABY_DEATH_RATE = 0 # Density independent new baby death
BIRTH_RATE = log(LITTER_SIZE) # Gives birth rate. Little r in ricker function
BETA = BIRTH_RATE / K_CAPACITY
# TODO: CHECK WHY I DID THIS?


# PARASITE PARAMETERS
# ENVIRONMENTAL_POOL = 0
# EGG_PRODUCTION_PER_WORM = 100 # eggs per month (note this is very wrong)
ENCOUNTER_MEAN = 500 # Mean number of encountered eggs
ENCOUNTER_K = 1 # Aggregation parameter of the NBD

# A set of parameters that determines the
# probability of encounter given a weighted history of prevalence. The first
# parameter dictates at what level of the log metric eprob is close to 1 and
# the second parameter dictates how quickly eprob goes to 1.
# TODO: WE WILL NEED TO FIDDLE WITH THESE
ENCOUNTER_PARAMS = c(2, 5)

INFECTIVITY = 0.02 # Probability of infectivity
RESISTANCE = 0.03 # How quickly a raccoon gains resistance based on previous load
EGG_DECAY = 0.3 # Rate of egg decay such that 3% chance of survival after year

# See fit_param.R for how we got these values
WORM_SURV_TRESH = 4.7104 #/ 2 # Threshold parameter of worm survival probability
WORM_SURV_SLOPE = -0.9446 #/ 2 # Slope of worm surv probability

# Time parameters: Each time step is a month
TIME_STEPS = 100