## ABC-SMC algorithm for estimating two/three parameters of RacoonWormModel

# Step 1: Draw parameters from prior distributions






# Step 2: Simulate model with drawn parameters
# Step 3: Calculate distance between observed and simulated data
#   - Do this in such a way that we can combine the two measures
# Step 4: Repeat for X number of parameter draws
# Step 5: Keep only the parameter combinations with 10% lowest distances
# Step 6: Use importance resampling, a bootstrap filter and some perturbations
# to smooth the particles and redraw from the weighted sample 
# Repeat until the distribution does not change. \
