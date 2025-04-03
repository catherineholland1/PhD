################################.
## BETA-BINOMIAL DISTRIBUTION ##
################################.


## BETA-BINOMIAL DISTRIBUTION NIMBLE FUNCTION ####


dbetabinomial = nimbleFunction(
  
  run = function(x = double(0),
                 nu = double(0),
                 phi = double(0),
                 size = double(0),
                 log = integer(0)){  
    
    returnType(double(0))
    
    # UPPER LIMIT TO ENSURE STABIITY IN THE DISTRIBUTION
    #phi <- min(phi,1e+04)    
    
    if(x >= 0 & x <= size){  
      
      return(lgamma(size + 1) + lgamma(x + nu * phi) + lgamma(size - x + (1 - nu) * phi) + lgamma(phi) -                
               lgamma(size + phi) - lgamma(nu * phi) - lgamma((1 - nu) * phi) - lgamma(size - x + 1) - lgamma(x + 1))     
      
    }else{       
      return(-Inf)     
    }   
  },
  buildDerivs = TRUE)     


rbetabinomial = nimbleFunction(
  
  run = function(n = integer(0),
                 nu = double(0),
                 phi = double(0),
                 size = double(0)){   
    
    returnType(double(0)) 
    
    # UPPER LIMIT TO ENSURE STABIITY IN THE DISTRIBUTION
    phi <- min(phi, 1e+04) 
    pi = rbeta(1, nu * phi, (1 - nu) * phi)     
    
    return(rbinom(1, size, pi))   
    
  })      


## REGISTER DISTRIBUTION ####
registerDistributions(
  list(dbetabinomial = list(BUGSdist='dbetabinomial(nu, phi, size)', discrete=TRUE))
)


#######################.