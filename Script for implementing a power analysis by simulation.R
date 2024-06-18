# ================================
#   
#   Power analysis presentation
# 
# ================================

# The purpose of this file is to present an R code which can be used to conduct power
# analyses by simulation in a regression context.
# The basic idea is to repeatedly generate a dataset given certain assumed parameters
# such as the sample size, N, and the effect size, b1. For each dataset we estimate a
# regression and get the p-values. The power of the test is equal to the number of p-values
# below the critical threshold (0.05, conventionally) divided by the total number of simulations.



#___________________ Simulating data _________________________ #

# First we're going to simulate data and run a regression on the data, just to get
# an idea of how that process works. We generate data featuring two variables: 
# individual height, our dependent variable, and individual milk consumption, our independent variable.

# Generate data
b_1 = 0.1 # We assume that a one unit increase in milk consumption increases height by 0.1
N = 750 # set sample size
df <- data.frame( # generate data
  milk <- rnorm(n = N, mean = 70, sd = 10),
  height <- 110 + b_1*milk + rnorm(n = 500, mean = 0, sd = 8)
); names(df) <- c("milk","height")

# Now we run a regression and get a (simple) summary
coef(summary(lm(height ~ milk, data = df))) # the estimated coefficeint on milk, b_1, 
# should be approx. 0.1



#___________________ Power analysis on simulated data _________________________ #

# There are two functions we must specify. 1) a data- and statistic-generating function
# which will generate a simulated data set according to given parameters (N, effect
# size etc) and give a p-value. 2) a simulation function which will run the
# function specified in 1) many times to generate a distribution of p-values. The power is
# the number of p-values below 0.05 divided by the number of simulations.


## specify data- and statistic-generating function
# this function generates simulted data according to given parameters (N and b1),
# runs a regression, and returns the p-value for b1 from that regression.

# define function
power_function <- function(N, b_1) {
  
  # generate data
  df <- data.frame( 
    milk <- rnorm(N,70,10),
    height <- 100 + b_1*milk + rnorm(N,0,8)
  ); names(df) <- c("milk","height")
  
  # run regression and return p-value
  pvalue_b1 <- coef(summary(lm(height ~ milk, data = df)))[2,4] 
  return(pvalue_b1)
}

# test function
power_function(N, b_1) # the returned value is the pvalue for b_1


## specify simulation function
# this function, power_sim, executes the power_function repeatedly (the total number of exections
# is determined by the variable nsims). For each execution, it extracts the pvalue 
# and saves it the vector pvals. When this is finished, it estimates the power (pow) by
# generating the proportion of pvalues in the pvals vectors which are less than 0.05 (the 
# convetional threshold)

# define function
power_sim <- function(N, b_1, nsims) {
  
  pvals <- replicate(
    n = nsims,
    power_function(N, b_1)
  )
  pow <- sum(pvals < .05) / nsims
  return(pow)
}

# test function (note: the higher nsims is, the closer the estimated power will be to the
# true value, but the longer the simulation process takes)
set.seed(130593)
power_sim(N = 500, b_1 = 0.1, nsims = 1000) # should get about 0.8
power_sim(N = 500, b_1 = 0, nsims = 1000) # When the true effect, b_1, is zero, 
# the power should be about 0.05, the probability of a false positive when the null is true


## visualise pvalues
# we might also wish to visualise the pvalues in a histogram, in which case we generate
# the vector of pvalues and visualise
pvals <- replicate(
  n = 1000,
  power_function(N = 500, b_1 = 0.1)
)
hist(pvals, breaks = 100, main = "distribution of p-values")

# the greater the power, the more right-skewed the function. If the true effect is
# zero, the distribution of pvalues should be approximately uniform
pvals <- replicate(
  n = 1000,
  power_function(N = 500, b_1 = 0)
)
hist(pvals, breaks = 100, main = "distribution of p-values")




### Multiple regression
# the above example was for the case of bivariate regression. It is straightforward to
# add control variables to the simulation

# specify data- and statistic-generating function
power_function <- function(N, b_1) {
  
  df <- data.frame( # generate data
    milk <- rnorm(N,70,10),
    exercise <- rnorm(N,4,1) + 0.1*milk, # we add this regressor, correlated with milk
    height <- 100 + b_1*milk + exercise + rnorm(N,0,8)
  ); names(df) <- c("milk","exercise","height")
  
  pvalue_b1 <- coef(summary(lm(height ~ milk + exercise, data = df)))[2,4] # run regression and get (simple) summary
  return(pvalue_b1)
}


# specify simulation function
power_sim <- function(N, b_1, nsims) {
  pvals <- replicate(
    n = nsims,
    power_function(N, b_1)
  )
  pow <- sum(pvals < .05) / nsims
  return(pow)
}

# run simulations
set.seed(130593)
power_sim(N = 500, b_1 = 0.1, nsims = 1000) 

# visualise pvalues
pvals <- replicate(
  n = 1000,
  power_function(N = 500, b_1 = 0.1)
)
hist(pvals, breaks = 100)




### Running a power analysis over a parameter space
# we would generally want to know the power for different values of N (the sample
# size), and b_1 (the effect size).

# Specify different parameter values and create empty vectors to record the value of the 
# parameter for a given simulation
n_values <- c(300, 500, 700) # values of N you want to estimate power for
b1_values <- c(0.3, 0.7, 0.1) # effect sizes you want to estimate power for
n_vector <- c() # vector to record values of N for each sim
b1_vector <- c() # vector to record values of effect size for each sim
power_vector <- c() # vector to register values of power for each sim

# run simulation for each combination of N and b1
for (N in n_values) {
  for (b_1 in b1_values) {
    
    # the following two lines are just useful to receord the simulations progress
    iteration_progress <- round(length(power_vector)/(length(n_values)*length(b1_values)), 2)
    print(paste("Iteration progress:", iteration_progress," Time:", Sys.time()))
    
    # run simulation and record parameter values
    power_it <- power_sim(N = N, b_1 = b_1, nsims = 1000) # this generates an estimate of the power for the given parameters
    power_vector <- append(power_vector, power_it) # this adds the power estimate to the power vector
    n_vector <- append(n_vector, N)
    b1_vector <- append(b1_vector, b_1)
  }
}

# generate a data frame to record the parameter values for each simulation
powerdf <- data.frame(power_vector, n_vector, b1_vector)

