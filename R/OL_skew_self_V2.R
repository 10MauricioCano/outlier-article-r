##############################################
#  Skewness projection directions for outlier detection
#  Outputs adapted to kurtosis routines
##############################################s

################### Function Sk.Child

Sk.Child = function(X, Proy){
  #
  # 'Sk.Child' computes random and specific directions from the maximum skewness direction 
  # of a multivariate dataset, making use of the maximum projection
  # method.
  #
  # Required input arguments:
  #     X : Matrix of size 'n*p' with the multivariate data; observations
  #     by rows and variables by columns ('X' must be standardized).
  #     Proy: Maximum skewness proyection
  #
  # Outputs:
  #     chil.proy : Skewness-based random directions children
  #
  # Author: Santiago Ortiz (santiagoortiz00@usc.edu.co - sortiza2@eafit.edu.co - saortizar@unal.edu.co - santy_ortiz@hotmail.com)
  # Date: 07/2023
  p = ncol(X)
  perc = as.vector(quantile(Proy, probs = c(0.25, 0.75)))
  rand.low = sample(which(Proy <= perc[1]), 2*p)
  rand.upp = sample(which(Proy >= perc[2]), 2*p)
  sem.low = X[rand.low,]
  sem.upp = X[rand.upp,]
  chil.proy = c()#matrix(0 ,4*(p^2), p)
  for (i in 1:(2*p)) {
    ch1 = matrix(rep(sem.low[1,], 2*p), 2*p, p, byrow = T) - sem.upp
    ch.norm = matrix(rep(sqrt(apply(ch1^2, 1, sum)), 2*p),  2*p, p, byrow = F)
    ch1.unit = ch1 / ch.norm
    chil.proy = rbind(chil.proy, ch1.unit)
  }
  return(t(chil.proy))
}

################### Function adj.outly

adj.outly = function(X){
  #
  # 'adj.outly' computes the adjusted outlyingness measure 
  # of a univariate dataset.
  #
  # Required input arguments:
  #     X : Vector of size 'n' with the univariate data.
  #
  # Outputs:
  #     Adj.AO : Outlyingness measure for each observation
  #
  # Author: Santiago Ortiz (santiagoortiz00@usc.edu.co - sortiza2@eafit.edu.co - saortizar@unal.edu.co - santy_ortiz@hotmail.com)
  # Date: 07/2023
  mu.rob = median(X)
  mc = medcouple(X)
  iqr = IQR(X)
  n = length(X)
  Adj.AO = rep(0, n)
  if (mc >= 0){
    w1 = quantile(X, 0.25) - (1.5*exp(-4*mc)*iqr)
    w2 = quantile(X, 0.75) + (1.5*exp(3*mc)*iqr)
  } else {
    w1 = quantile(X, 0.25) - (1.5*exp(-3*mc)*iqr)
    w2 = quantile(X, 0.75) + (1.5*exp(4*mc)*iqr)
  }
  for (i in 1:n) {
    if (X[i] > mu.rob){
      Adj.AO[i] = (X[i]- mu.rob) / (w2 - mu.rob)
    } else {
      Adj.AO[i] = (mu.rob - X[i]) / (mu.rob - w1)
    }
  }
  return(Adj.AO)
}

################### Function Adj.Outlier

Adj.Outlier = function(X1, type){
  #
  # 'Adj.Outlier' computes the univariate outlier detection based on
  # the Skewness Adjusted Outlyingness of a univariate dataset.
  #
  # Required input arguments:
  #     X1 : Vector of size 'n' with the univariate data.
  #     type: scalar, "1" for Medcouple calibration, otherwise traditional upper whisker
  #
  # Outputs:
  #     lab.out : Binary, "1" Outlier, "0" Non-outlier
  #
  # Author: Santiago Ortiz (santiagoortiz00@usc.edu.co - sortiza2@eafit.edu.co - saortizar@unal.edu.co - santy_ortiz@hotmail.com)
  # Date: 07/2023
  X = X1[which(X1 <= quantile(X1, 0.6))]
  n = length(X1)
  lab.out = rep(0, n)
  mc = medcouple(X)
  iqr = IQR(X)
  p75 = quantile(X, 0.75)
  if (type == 1){
    cutoff = p75 + (1.5*exp(3*mc)*iqr)
  } else {
    cutoff = p75 + (1.5*iqr)
  }
  outliers = which(X1 >= cutoff)
  lab.out[outliers] = 1
  return(lab.out)
}

################### Function maxskew

max_skew <- function(X) {
   
   #
   # max_skew computes the maximum skewness projection of a multivariate
   # dataset. This function requires the function 'val_skew.m' and the
   # R library 'mrfDepth' (function 'medcouple')
   #
   # Required input arguments:
   #     X : Matrix of size 'n*p' with the multivariate data. observations
   #     by rows and variables by columns.
   #
   # Outputs:
   #     dir_vec : Vector of size 'p*1' that maximize the skewnees coefficient
   #     of the data in x_in.
   #
   #     skew_val : Value of the third standardized central moment of the
   #     projection of x_in in dir_vec.
   #
   # Author: Santiago Ortiz (santiagoortiz00@usc.edu.co - sortiza2@eafit.edu.co - saortizar@unal.edu.co - santy_ortiz@hotmail.com)
   # Date: 02/2020

   ## Initial verification
     
   if (missing(X)) stop("Input argument X is undefined")
   
   tol = 1.0e-10

   X.dm = dim(X)
   n.x = X.dm[1]
   p.x = X.dm[2]
   
   x = scale(X)                   # Standardized Data
   x1 = t(x)                      # Transposed Data

   ## Projection vector's Initializer - Taken from:
   ##    Peña & Prieto (2001), Multivariate Outlier Detection and Robust Covariance
   ##    Matrix Estimation

   # Normalizing data by columns
   uv = apply(x*x,1,sum)
   uw = 1.0/(tol + sqrt(uv))
   uu = x * uw

   # Computing the Principal Components of normalized data in uu
   Sue = eigen(cov(uu))
   V = Sue$vectors

   # Checking which eigenvalue has the biggest absolute third centered moment
   r = matrix(0,1,p.x)
   for (i in 1:p.x) {
      d.i = as.matrix(V[,i])
      r[i] = val_skew(x,d.i)
   }
   
   # Eigenvector with the biggest third moment, this is the INITIALIZER
   ik = order(abs(r))
   d0 = as.matrix(V[,ik[p.x]])
   if (val_skew(x,d0) < 0) d0 = -d0

   ## Search Algorithm
   d1 = matrix(0,p.x,1)                # Dummy vector to initialize
   cont = 0                            # Counter to stop loop if the skewness does not converge
   # Initialize direction before the loop starts, using the vector d0 as the best projection vector
   best_dir = d0
   tolerance = 1e-4                    # Tolerance level to stop the loop
   
   ## We need this for the corrected procedure
   S.sqr = pracma::sqrtm(cov(x))
   S.sqri = S.sqr$Binv

   while (norm(d0 - d1) > tolerance) {
      cont = cont + 1

      # Skewness matrix
      w.vec = t(x1) %*% d0
      aux.vec = c(w.vec) * t(x1)
      M_D = x1 %*% aux.vec
      ## Corrected procedure
      M_D = S.sqri %*% M_D %*% S.sqri

      eig.v = eigen(M_D)                      # Eigenvalues and Eigenvectors of matrix M_D
      ind = order(abs(eig.v$values))          # Sort the eigenvalues in ascending order
      eig.vec.sorted = eig.v$vectors[,ind]    # Ordered eigenvectors, greatest eigenvector in the last column

      d1 = d0
      d0 = S.sqri %*% as.matrix(eig.vec.sorted[,p.x])    # New projection vector in the loop
      d0 = d0/norm(d0,type="F")
      
      if (cont == 1) best_min = d0
      if ((cont >= 2) && (abs(medcouple(x %*% d0)) < abs(medcouple(x %*% best_min)))) best_min = d0
      # Saving the vector which has the biggest third central moment
      if (abs(medcouple(x %*% d0)) > abs(medcouple(x %*% best_dir))) best_dir = d0
      if (cont == 200) break         # If in 200 iterations there's not convergence, the loop breaks
   }

   # Selecting the best projection vector
   dir_saved = d0

   x.bd = X %*% best_dir
   x.bm = X %*% best_min
   x.ds = X %*% dir_saved
   am.bd = abs(medcouple(x.bd))
   am.bm = abs(medcouple(x.bm))
   am.ds = abs(medcouple(x.ds))
   as.bd = abs(skewness(x.bd))
   as.bm = abs(skewness(x.bm))
   as.ds = abs(skewness(x.ds))
   
   if ((am.bd > am.bm && am.bd > am.ds) || (as.bd > as.bm && as.bd > as.ds)) {
      best_of_all = best_dir
      x.ba = x.bd
   }
   else if (am.bm > am.bd && am.bm > am.ds) {
      best_of_all = best_min
      x.ba = x.bm
   }
   else {
      best_of_all = dir_saved
      x.ba = x.ds
   }

   if (skewness(x.ba) > 0) dir_vec = best_of_all
   else dir_vec = -best_of_all

   # Computing the third standardized central moment of the projection of 'X' in 'dir_vec'
   skew_v = val_skew(X,dir_vec)
   
   # Return values
   rval = list(dv = dir_vec, sv = skew_v)
   return(rval)
}

################### Function val_skew

val_skew <- function(x,d,km = 3) {
   
   #
   #     mc = val_skew(x,dir,k)
   #
   # Evaluate the moment coefficient of order k
   # for the univariate projection of multivariate data
   #
   # Inputs:  x, observations (by rows)
   #          dir, projection direction
   #          k, order of the moment (k=3 by default)
   # Output:  mc, value of the moment coefficient for the
   #          univariate data
   #
   #
   # Daniel Pena/Francisco J Prieto 23/5/00

   vs.xd = dim(x)
   vs.n = vs.xd[1]
   vs.p = vs.xd[2]
   
   vs.dd = dim(d)
   vs.p1 = vs.dd[1]
   
   if (vs.p != vs.p1) stop("Data dimensions are not correct")

   vs.t = x %*% d
   vs.tm = mean(vs.t)
   vs.tt = abs(vs.t - vs.tm)
   vs.vr = sum(vs.tt^2)/(vs.n-1)
   vs.kr = sum(vs.tt^km)/vs.n
   vs.mc = vs.kr/vs.vr^(km/2)

   return(vs.mc)
}

################### Function gen_rcorr

gen_rcorr <- function(cond.S,p.x) {
  
  #
  # Generate random correlation matrix.
  #   R = random_correlation_generator(p, iteration_times) is a function
  #   that generate a random correlation estructure for a p-dimensional
  #   multivariate random variable. The random correlation structure
  #   is described in [1].
  #
  #   [1] Agostinelli, C., Leung, A., Yohai, V.J., Zamar, R.H., 2015. Robust
  #        estimation of multivariate location and scatter in the presence of
  #        cellwise and casewise contamination. TEST 24(3), 441-461.
  #   INPUTS:
  #       p: Is the dimension of the normal multivariate random variable.
  #       iteration_times (Not required): Is the maximum iteration times for
  #       the convergence of the algorithm to estimate R, default 99 times.
  #   OUTPUTS:
  #       R: Is the correlation matrix that contain the random correlation
  #       estructure generated.
  #   Autor: Henry G. Velasco (hgvelascov@eafit.edu.co)
  
  maxits = 100
  tol = 1.0e-5
  
  if (p.x < 3) stop("p must be larger than 2")
  
  lambda = sort(c(1,runif(p.x-2,1,cond.S),cond.S))
  x = matrix(rnorm(p.x*p.x),p.x,p.x)
  Sigma = x %*% t(x)
  
  S.eig = eigen(Sigma)
  Q.S = S.eig$vectors
  ratio = 0
  iter = 0
  
  while ((abs(ratio - cond.S) > tol) && (iter < maxits)) {
    iter = iter + 1
    Sigma = Q.S %*% diag(lambda) %*% t(Q.S)
    Sig.sr = 1.0/sqrt(c(diag(Sigma)))
    Sigma = t(Sig.sr * Sigma) * Sig.sr
    S.eig = eigen(Sigma)
    Q.S = S.eig$vectors
    lambda = S.eig$values
    ratio = lambda[1]/lambda[p.x]
    lambda[p.x] = lambda[1]/cond.S
  }
  
  return(Sigma)
}

################### Function GenAtip

GenAtip <- function(n.x,p.x,par.lst,sim.mode = 1) {
  
  #
  # Generate contaminated observations controlling the parameters
  #
  #     GenAtip(n.x,p.x,par.lst,sim.mode)
  #
  # Inputs:  n.x,  number of observations
  #          p.x,  dimension of each observation
  #          par.lst, parameters to use in the data generation model
  #          sim.mode, = 1 , normal observations and contamination
  # Output:  x, observations
  #          lbl, labels of the observations (lbl = 0 for an uncontaminated observation
  #                                           lbl = 1 for an outlier)
  #
  
  #  Daniel Pena / Francisco J Prieto 06/03/2020
  
  ## Labels for outliers and central observations
  
  lbl = matrix(1,n.x,1)
  
  ## Generate outliers according to the selected pattern
  
  if (sim.mode == 1) {
    
    ## Usual contamination model (1 group of outliers)
    
    ## Parameters
    
    alpha = par.lst[1]
    n.1 = floor(n.x*alpha)
    n.0 = n.x - n.1
    
    lbl[1:n.0] = 0
    
    ## Centers of clusters
    
    delta = par.lst[2]            # Default value = 1
    m.0 = matrix(0,1,p.x)
    m.1 = delta*matrix(2,1,p.x)
    
    ## Data for each cluster
    
    lambda = par.lst[3]
    x = matrix(rnorm(n.x*p.x),n.x,p.x)
    x[1:n.0,] = scale(x[1:n.0,])
    x[(n.0+1):n.x,] = scale(x[(n.0+1):n.x,])
    x[(n.0+1):n.x,] = sqrt(lambda)*x[(n.0+1):n.x,] + matrix(1,n.1,1) %*% m.1
    
  } else if (sim.mode == 2) {
    
    ## Contamination type A from HL and SO
    
    cn.mx = 100
    cn.mn = 1
    dst.mx = 10
    dst.mn = 5
    
    n.gr = 2
    
    alpha = par.lst[1]
    n.cnt.g = round(n.x*alpha/n.gr)
    n.cln = n.x - n.gr*n.cnt.g
    
    delta = par.lst[2]        # Default value = 1
    lambda = par.lst[3]       # Default value = 1
    
    lbl[1:n.cln] = 0
    
    ## Observations in each cluster
    
    Cov.cln = gen_rcorr(cn.mx,p.x)
    mn.cln = matrix(0,p.x,1)
    x.cln = mvrnorm(n.cln,mn.cln,Cov.cln)
    
    Cov.cnt = lambda*diag(p.x)
    v.1 = runif(p.x-1,dst.mn,dst.mx)
    v.2 = -runif(p.x-1,dst.mn,dst.mx)
    mn.cnt.1 = delta*matrix(c(v.1,0),p.x,1)
    mn.cnt.2 = delta*matrix(c(v.2,0),p.x,1)
    x.cnt.1 = mvrnorm(n.cnt.g,mn.cnt.1,Cov.cnt)
    x.cnt.2 = mvrnorm(n.cnt.g,mn.cnt.2,Cov.cnt)
    
    ## Simulated data
    
    x = rbind(x.cln,x.cnt.1,x.cnt.2)
    
  } else if (sim.mode == 3) {
    
    ## Contamination type B from HL and SO
    
    cn.mx = 100
    cn.mn = 1
    dst.mx = 10
    dst.mn = 5
    
    n.gr = 4
    
    alpha = par.lst[1]
    n.cnt.g = round(n.x*alpha/n.gr)
    n.cln = n.x - n.gr*n.cnt.g
    
    delta = par.lst[2]         # Default value = 1
    lambda = par.lst[3]        # Default value = 1
    
    lbl[1:n.cln] = 0
    
    ## Observations in each cluster
    
    Cov.cln = gen_rcorr(cn.mx,p.x)
    mn.cln = matrix(0,p.x,1)
    x.cln = mvrnorm(n.cln,mn.cln,Cov.cln)
    
    Cov.cnt = lambda*diag(p.x)
    v.1 = runif(p.x-1,dst.mn,dst.mx)
    v.2 = -runif(p.x-1,dst.mn,dst.mx)
    v.3 = c(matrix(c(1,-1),p.x-1,1))*runif(p.x-1,dst.mn,dst.mx)
    v.4 = c(matrix(c(-1,1),p.x-1,1))*runif(p.x-1,dst.mn,dst.mx)
    mn.cnt.1 = delta*matrix(c(v.1,0),p.x,1)
    mn.cnt.2 = delta*matrix(c(v.2,0),p.x,1)
    mn.cnt.3 = delta*matrix(c(v.3,0),p.x,1)
    mn.cnt.4 = delta*matrix(c(v.4,0),p.x,1)
    x.cnt.1 = mvrnorm(n.cnt.g,mn.cnt.1,Cov.cnt)
    x.cnt.2 = mvrnorm(n.cnt.g,mn.cnt.2,Cov.cnt)
    x.cnt.3 = mvrnorm(n.cnt.g,mn.cnt.3,Cov.cnt)
    x.cnt.4 = mvrnorm(n.cnt.g,mn.cnt.4,Cov.cnt)
    
    ## Simulated data
    
    x = rbind(x.cln,x.cnt.1,x.cnt.2,x.cnt.3,x.cnt.4)
  }
  
  ## Return values
  cv = list(x = x, lbl = lbl)
  return(cv)
}


###############################################################
## Generador FICM  (α = probabilidad celda-a-celda)----
###############################################################
GenAtip_FICM <- function(n.x, p.x, par.lst, sim.mode = 1) {
  # par.lst = c(alpha_cell, delta, lambda)
  alpha_cell <- par.lst[1]
  delta      <- par.lst[2]
  lambda     <- par.lst[3]
  
  ## 1) datos limpios
  Y <- matrix(rnorm(n.x * p.x), nrow = n.x, ncol = p.x)
  
  ## 2) máscara de celdas contaminadas
  B <- matrix(rbinom(n.x * p.x, 1, alpha_cell), nrow = n.x, ncol = p.x)
  
  ## 3) medias: contaminación uniforme en dirección +delta (modo 1 únicamente)
  mu_vec <- rep(delta, p.x)
  
  ## 4) genera Z sólo donde B == 1
  idx <- col(B)[B == 1]                 # columnas de las celdas TRUE
  Z   <- rnorm(length(idx),             # mismo largo que sum(B)
               mean = mu_vec[idx],
               sd   = sqrt(lambda))
  
  ## 5) ensambla resultado y etiquetas
  X           <- Y
  X[B == 1]   <- Z
  lbl         <- as.integer(rowSums(B) > 0)   # 0 = fila limpia
  
  attr(X, "mask") <- B
  list(x = X, lbl = lbl)
}


######################################################################
################### Main code----
######################################################################

library(mrfDepth)             # Function medcouple
library(e1071)                # Function skewness
library(pracma)               # Function findpeaks
library(MASS)                 # Function mvrnorm

#########################
#  Simulation runs
#  Choice of contamination model
#########################

# Simulation parameters

mode.sim = 3    # Simulation mode
# = 1 one oultier group
# = 2 two outlier groups (first model from Henry and Santiago)
# = 3 four outlier groups (second model from Henry and Santiago)

alpha = 0.1
delta = 1            # Default value for methods 2 and 3
lambda = 0.1         # Only used for method 1
par.lst = c(alpha,delta,lambda)

n.x = 500
p.x = 3

x.val = GenAtip(n.x,p.x,par.lst,mode.sim)
#aa = out_skew(x.val$x)
