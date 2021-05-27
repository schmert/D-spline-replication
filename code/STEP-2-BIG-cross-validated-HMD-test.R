#---------------------------------------------------------------
# large Monte Carlo test of P-spline(MortalitySmooth), D1,D2,DLC
# estimators on simulated small-area data
# 
# for each D-spline method there is a separate set of penalty
# constants for each country, constructed from all the HMD
# schedules EXCEPT that country
# 
# for each combination of (schedule, sample size) there are
# nsize (probably=100) simulations. Thus 
# 222 schedules x 2 sample sizes x 100 simulations = 44,400 trials
# 
# in order to bypass GitHub and GitLab file size limits,
# the saved output is divided into three parts. For example,
# if the program starts saving output at 1:08pm on 1 Feb 2021,
# the output files will be called
#     BIG-cross-validated-HMD-test-2021-02-01-1308-part-1.Rdata
#     BIG-cross-validated-HMD-test-2021-02-01-1308-part-2.Rdata
#     BIG-cross-validated-HMD-test-2021-02-01-1308-part-3.Rdata
#
# This program may take a long time to run: on a standard
# desktop PC (circa 2020) the run time was nearly 1 hour
#---------------------------------------------------------------

library(MASS)
library(MortalitySmooth)
library(tidyverse)

rm(list = ls())

set.seed(6447100)


## load some utility functions
source('Dspline_fit_function.R')
source('vech-xpnd-functions.R')

## load the HMD data (female 1x10 lifetables, exposures, and deaths)
load('../data/HMD.RData')


# which decade corresponds to each HMD schedule?
endyear_vec = strsplit(colnames(HMD_N),"-") %>% 
  sapply(function(x) as.numeric(x[[2]]))


#' select a historical cutoff -- only schedules that start after a certain
#' year

earliest_year = 1970

keep = which(endyear_vec > earliest_year) 

# only keep schedules for which all years are after the cutoff

sel_logmx = HMD_logmx[,keep]
sel_N     = HMD_N[,keep]
sel_D     = HMD_D[,keep]

nsched = ncol(sel_logmx)

# which country corresponds to each HMD schedule,
# among those retained

country_vec = strsplit(colnames(sel_N),":") %>% 
  sapply(function(x) x[[1]]) 

## samples
ntrials          = 100            # trials per HMD schedule
SAVE.RESULTS     = TRUE
sample_size_vals = c(10000, 100000)
identical_bases  = TRUE           # use Mortality1DSmooth B matrix?


get_sel_logmx = function(ix) {
  sel_logmx[, ix]
}

get_sel_N = function(ix) {
  sel_N[, ix]
}

D1 = diff( diag(100), diff=1)

D2 = diff( diag(100), diff=2)

B1     = splines::bs(x=0:98, knots=c(seq(0,20,2),seq(30,100,10)), degree=2)
Proj1  = B1 %*% ginv(crossprod(B1)) %*% t(B1)

mean_first_diff = function(ix, smoothed=TRUE) {
  c1 = rowMeans(D1 %*% sel_logmx[,ix])
  if (smoothed) c1 = Proj1 %*% c1
  return(as.vector(c1))
}

# returns the lower-triangular part of the inverse variance
# matrix of residuals
inv_var_first_diff = function(ix,c1) {
  V = var( t( D1 %*% sel_logmx[,ix] - c1 ))
  return(ginv(V) %>% vech() )
}

B2     = splines::bs(x=0:97, knots=c(seq(0,20,2),seq(30,100,10)), degree=2)
Proj2  = B2 %*% ginv(crossprod(B2)) %*% t(B2)

mean_second_diff = function(ix, smoothed=TRUE) {
  c2 = rowMeans(D2 %*% sel_logmx[,ix])
  if (smoothed) c2 = Proj2 %*% c2
  return(as.vector(c2))
}

inv_var_second_diff = function(ix,c2) {
  V = var( t( D2 %*% sel_logmx[,ix] - c2 ))
  return(ginv(V) %>% vech() )
}


LC_A = function(ix) {
  L = sel_logmx[,ix]
  a = rowMeans(L)
  tmp = svd( sweep(L,1,a,'-') )
  b   = tmp$u[,1]
  M   = diag(length(b)) - tcrossprod(b)
  
  return(M %>% vech())
}

LC_c = function(ix) {
  L = sel_logmx[,ix]
  a = rowMeans(L)
  tmp = svd( sweep(L,1,a,'-') )
  b   = tmp$u[,1]
  M   = diag(length(b)) - tcrossprod(b)
  
  return(as.vector(M %*% a))
}


inv_var_LC_resid = function(ix,A,c) {
  V = var( t( xpnd(A) %*% sel_logmx[,ix] - c ))
  return(ginv(V) %>% vech())
}


make_sample_N = function(ix,sample_size) {
    tmp = as.integer( round( sample_size * prop.table(sel_N[,ix]) ))
    return(as.vector(tmp))
}


make_sample_D = function(ix,N) {
  repeat {
    mu     = exp( sel_logmx[,ix] )
    lambda = N * mu
    res    = rpois(length(lambda), lambda)
    
    if (sum(res) > 1) {break}  # ensure at least 2 deaths
  }
  return(as.vector(res))

}

## sample_structure holds the indices of the
## training schedules and the test schedules
## Dspline constants for each trial 
## will be developed from training subset of HMD


#' all estimators will maximize
#'    L(theta) -1/2 (A B theta - c)' Vinv (A B theta - c)
#' over theta   
#' for 1st diffs -- A = D1 , c = mean 1st diffs,
#' for 2nd diffs -- A = D2 , c = mean 2nd diffs 
#' for LC        -- A = Mb , c = Mb a

# calculate Dspline penalty vectors and matrices for each country
# by omitting that country and using data from the rest of the HMD
# schedules
# 
# the result, Dspline_constants, will be a data frame with
# one row per COUNTRY, columns that are LISTS
# for example, the first row will have country=AUS (Australia)
# and its variables will be 
#    c1 (a 99-vector of avg first differences in non-Australian HMD logmx schedules)
#    V1 (a 4950-vector containing the diagonal and lower
#        elements of the 99x99 covariance matrix of residuals
#        for first-differences in logmx)
#   etc. 
#   
#   MAIN POINT: variables in this tibble/data frame are LIST COLUMNS       
#
#  the temporary variable train_ix below is a vector
#  containing the index of all of the HMD_logmx columns
#  that do NOT correspond to the country in question 
#  
#     country c1         V1            c2         V2            A3            c3          V3           
#     <chr>   <list>     <list>        <list>     <list>        <list>        <list>      <list>       
#   1 AUS     <dbl [99]> <dbl [4,950]> <dbl [98]> <dbl [4,851]> <dbl [5,050]> <dbl [100]> <dbl [5,050]>
#   2 AUT     <dbl [99]> <dbl [4,950]> <dbl [98]> <dbl [4,851]> <dbl [5,050]> <dbl [100]> <dbl [5,050]>
#   3 BEL     <dbl [99]> <dbl [4,950]> <dbl [98]> <dbl [4,851]> <dbl [5,050]> <dbl [100]> <dbl [5,050]> 
#...
#

Dspline_constants = tibble(
    country  = unique(country_vec), 
    train_ix = map(country, function(cc) which(country_vec != cc))
  ) %>% 
  mutate(
    c1          = map (train_ix, mean_first_diff, smoothed=TRUE),
    V1          = map2(train_ix, c1, inv_var_first_diff),
    c2          = map (train_ix, mean_second_diff, smoothed=TRUE),
    V2          = map2(train_ix, c2, inv_var_second_diff),
    A3          = map (train_ix, LC_A),
    c3          = map (train_ix, LC_c),
    V3          = pmap(list(train_ix, A3, c3), inv_var_LC_resid )
  ) %>% 
  dplyr::select(-train_ix)


# create ntrials samples for each (schedule x sample size) pair
# ix represents the column number of the schedule in the HMD_logmx
# matrix that is save in HMD.RData
# 
# as above, the columns of sample are LIST COLUMNS
# 
# A tibble: 44,400 x 5
#        trial ix sample_size  N          D          
#       <int> <int>     <dbl>  <list>     <list>     
#  1     1     1       10000  <int [100]> <int [100]>
#  2     2     1       10000  <int [100]> <int [100]>
#  ...
#  
#  
samples = expand.grid(
            trial       = seq(ntrials), 
            ix          = seq(nsched),
            sample_size = sample_size_vals
              ) %>% 
          as_tibble() %>% 
          mutate( N = map2(ix,sample_size, make_sample_N),
                  D = map2(ix,N, make_sample_D))
 

## calculate the true values of logmu, e0, e60 for each schedule in HMD
## and arrange into a df
e = function(logmx, start.x=0) {
  lx = exp( -cumsum(c(0,exp(logmx))))  # ages 0,1,...,100
  age = 0:100
  keep.lx = lx[ age >= start.x ] / lx[age == start.x]
  return ( sum((head(keep.lx, -1) + tail(keep.lx, -1))) / 2)
}

ex.star = sapply( c(0,60,80),
                  function(this.x) apply(sel_logmx,2,e,start.x=this.x))
colnames(ex.star) = c(0,60,80)

true_values = tibble(
  ix    = 1:nrow(ex.star),
  logmx = map(ix, function(i) sel_logmx[,i]),
  e0    = ex.star[,'0'],
  e60   = ex.star[,'60'],
  e80   = ex.star[,'80']
)

fit = function(N,D,ix,method) {

  # get the constants for this country
  x = filter(Dspline_constants, country==country_vec[ix])

  if (method %in% c('D1','D2','LC')) {
  
    if (method == 'D1') {
      this_A = D1
      this_c = x$c1[[1]]
      this_V = xpnd( x$V1[[1]] )
  } else if (method == 'D2') {
      this_A = D2
      this_c = x$c2[[1]]
      this_V = xpnd( x$V2[[1]] )
  } else if (method == 'LC') {
      this_A = xpnd( x$A3[[1]])
      this_c = x$c3[[1]]
      this_V = xpnd( x$V3[[1]] )
  }
    
    tmp = try({ 
        Dspline_fit(N,D,this_A, this_c, this_V, 
                    use_MS_basis = identical_bases,
                    details=TRUE, 
                    max_iter = 100, theta_tol = .005)
      }, silent=TRUE)
    
    # if nonconvergent, try again without the default Dspline basis
    if (class(tmp) == 'try-error') {
      try({ 
        Dspline_fit(N,D,this_A, this_c, this_V, 
                    use_MS_basis = FALSE,
                    details=TRUE, 
                    max_iter = 100, theta_tol = .005)
      }, silent=TRUE)
    }
    
    if (class(tmp) == 'try-error') {
        result=NA } else {  
        result = list(logmx_hat = as.vector(tmp$lambda.hat),
                      df        = tmp$df,
                      bic       = tmp$bic)
        }
    
  } else if (method=='MS') {
    tmp = try( {
      w = 1*(N>0) # disregard ages with zero exposure
      z = Mort1Dsmooth( x      = 0:99, 
                        y      = D,
                        w      = w,
                        offset = log(N),
                        ndx    = 33,
                        method = 1,
                        control = list(RANGE=c(1e0,1e6)))
     }, silent=TRUE) 
    

    if (class(tmp) == 'try-error') {
         result = NA
      } else {  
result = list(logmx_hat = as.vector(z$logmortality),
                        df        = z$df,
                        bic       = z$bic)
        }          
  
  } ## if  MS

 return(result)
}
 
grab_component = function(this_list, this_name) {
  lapply( this_list, function(x) unlist(x[this_name]))  
}

# the next step takes the sample data frame and fits all 4 models to
# each of the 44,400 (D,N) samples. It then separates out the logmx estimates, 
# the estimated degrees of freedom, and the estimated e0 and e60 from each 
# simulated sample
# 
# result is another data frame/tibble with LIST COLUMNS, one per sample
# like...
#     trial    ix sample_size method logmx_hat      df   bic e0_hat e60_hat
#     <int> <int>       <dbl> <chr>  <list>      <dbl> <dbl>  <dbl>   <dbl>
# 1     1     1       10000     D1     <dbl [100]>  1.60  84.8   76.7    20.7
# 2     1     1       10000     D2     <dbl [100]>  3.27  92.5   76.7    20.9
# 3     1     1       10000     LC     <dbl [100]>  1.39  82.8   76.7    20.7
# 4     1     1       10000     MS     <dbl [100]>  5.82 104.    76.8    20.5
# 5     2     1       10000     D1     <dbl [100]>  1.71 102.    74.8    19.6
# 6     2     1       10000     D2     <dbl [100]>  3.37 105.    73.8    20.3
# 7     2     1       10000     LC     <dbl [100]>  1.48 110.    75.0    19.8
# 8     2     1       10000     MS     <dbl [100]>  3.33 122.    73.8    20.2
# ...


system.time( {
  result = samples %>% 
             mutate(method = map(trial, function(i) c('D1','D2','LC','MS'))) %>% 
             unnest(col=method) %>% 
             mutate(this_fit  = pmap(list(N,D,ix,method), fit),       # <<< MAIN FITTING
                    logmx_hat = grab_component(this_fit, 'logmx_hat'),
                    df        = unlist(grab_component(this_fit, 'df')),
                    bic       = unlist(grab_component(this_fit, 'bic')),  
                    e0_hat    = map_dbl(logmx_hat, e, start.x=0),
                    e60_hat   = map_dbl(logmx_hat, e, start.x=60)
             ) %>% 
           dplyr::select(-this_fit, -N,-D)
})


if (SAVE.RESULTS) {

    timestamp = format(Sys.time(), "%Y-%m-%d-%H%M")

# save in parts in order to stay under
# 100 MB Github limit  

  k = 0
  for (ssize in sample_size_vals) {
     k = k+1
     subresult = filter(result, sample_size == ssize)
     save( subresult, 
           file=paste0('../data/BIG-cross-validated-HMD-test-',
                      timestamp,'-part-',k,'.Rdata'))
  }
  
  # info about how to map from this sample back to orig HMD data
  crosswalk = tibble(ix              = seq(1:nsched),
                     orig_HMD_column = keep,
                     sched           = colnames(sel_D), 
                     country         = country_vec
                     )
  
  save( true_values, samples,  Dspline_constants,
        crosswalk,
        file=paste0('../data/BIG-cross-validated-HMD-test-',
                    timestamp,'-part-',k+1,'.Rdata'))

} # if SAVE.RESULTS

