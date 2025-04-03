######################################.
## FORWARD ALGORITHM - BETA-BINOMAL ##
######################################.


## BETA-BINOMIAL DISTRIBUTION NIMBLE FUNCTION ####

dbetabinomial = nimbleFunction(
  
  run = function(x = double(0),
                 nu = double(0),
                 phi = double(0),
                 size = double(0),
                 log = integer(0)){  
    
    returnType(double(0))
    
    # UPPER LIMIT TO ENSURE STABIITY IN THE DISTRIBUTION
    phi <- min(phi,1e+04)    
    
    if(x >= 0 & x <= size){  
      
      return(lgamma(size + 1) + lgamma(x + nu * phi) + lgamma(size - x + (1 - nu) * phi) + lgamma(phi) -                
               lgamma(size + phi) - lgamma(nu * phi) - lgamma((1 - nu) * phi) - lgamma(size - x + 1) - lgamma(x + 1))     
      
    } else {       
      return(-Inf)     
    }   
  })     



## COLUMN SUM NIMBLE FUNCTION ####

column_sum = nimbleFunction(
  
  run = function(x = double(2)){
    
    y = numeric(dim(x)[2])
    
    for(j in 1:dim(x)[2]){
      
      y[j] = sum(x[,j])
      
    }
    
    returnType(double(1))
    return(y)
    
  })


## FORWARD ALGORITHM NIMBLE FUNCTION ####

forward_alg = nimbleFunction(
  
  run = function(p0 = double(1),
                 p = double(2),
                 n = double(0),
                 N_states = double(0),
                 dens = double(2)){
    
    c = numeric(n)
    
    c[1] = sum(dens[1, ] * p0)
    
    alpha = (dens[1, ] * p0) / c[1]
    
    for(t in 2:n){
      
      delta = dens[t, ] * column_sum(matrix(rep(alpha, N_states), ncol = N_states) * p)
      c[t] = sum(delta)
      alpha = delta / c[t]
      
    }
    
    returnType(double(0))
    return(sum(log(c)))
    
  })


## HMM DISTRIBUTION NIMBLE FUNCTION ####


### dhmm is the distribution you will use in nimble. It's a distribution for the whole 
### time series of the response variable, i.e. y[1:n] ~ dhmm(...)

# x is the response variable for the whole time series (can't have any missing values!!!)
# P_zero is the initial state probabilities
# P is the transition matrix
# N is the length of the time series
# mean is the mean of the Normal for each state
# sd is the standard deviation of the Normal for each state

dhmm_betabinomial = nimbleFunction(
  
  run = function(x = double(1),
                 y = double(1),
                 p0 = double(1),
                 p = double(2),
                 N = double(0),
                 N_states = double(0),
                 nu = double(1),
                 phi = double(1),
                 log = integer(0)) {
    
    dens = matrix(nrow = N, ncol = N_states)
    
    # Loop over states
    for(j in 1:N_states){
      
      for (w in 1:N) {
        
        dens[w, j] = exp(dbetabinomial(x[w], nu[j], phi[j], y[w], log = TRUE))
      
      }
    }
    
    returnType(double(0))
    return(forward_alg(p0, p, N, N_states, dens))
    
  }
)



rhmm_betabinomial = nimbleFunction(
  
  run = function(n = integer(0),
                 y = double(1),
                 p0 = double(1),
                 p = double(2),
                 N = double(0),
                 N_states = double(0),
                 nu = double(1),
                 phi = double(1)) {
    
    returnType(double(1))
    
    return(rep(0, n))
  }
)


## REGISTER DISTRIBUTION ####

# registerDistributions(list(
#   dhmm = list(
#     BUGSdist = 'dhmm(p0, p, N, N_states, mean, sd)',
#     types = c('value = double(1)','p0 = double(1)','p = double(2)', 'N = double(0)', 'N_states = double(0)','mean = double(1)','sd = double(1)'),
#     mixedSizes = TRUE
#   )
# ))

#######################.