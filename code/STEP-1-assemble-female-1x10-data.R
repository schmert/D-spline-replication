############################################################
# process a large # of 1x10 female tables from the HMD
# (already downloaded and unzipped)
#
# Calculate logmx by single year of age and decade, and
# then calculate residuals and covariances for each of
# four possible D-spline penalties
#
# Save results for use in other programs. These include
# the D-spline constants used in the penalties (D1,D2,DLC).
############################################################

rm(list=ls())
graphics.off()

library(tidyverse)

file.list = dir('../HMD-input-data/fltper_1x10')

HMD_logmx = matrix(NA,111,0)

for (this.file in file.list) {
  
  this.country = sub('.fltper_1x10.txt','',this.file)

  ## log mortality rates    
  tmp = read.table(paste0('../HMD-input-data/fltper_1x10/',this.file), 
                   stringsAsFactors = FALSE,
                  skip=2, header=TRUE)
  tmp$Age = 0:110
  
  # despite using library(tidyverse) above, it's safer to be redunant
  # in the next command and use dplyr::select 
  
  tmp = tmp %>%
         mutate(logmx=log(mx)) %>%
         dplyr::select(Year,Age,logmx)
  
  HMD_logmx = cbind( HMD_logmx, 
             matrix(tmp$logmx, nrow=111, 
                    dimnames=list(0:110,
                          paste(this.country,unique(tmp$Year),sep=":")))
             )
  
}



# quick and dirty: keep only columns without zeroes, and 
# only ages < 100

ok = apply(HMD_logmx, 2, function(x) all(is.finite(x)))

HMD_logmx = HMD_logmx[paste(0:99),ok]

##################################################################
#  get exposure and death data for the selected schedules
##################################################################

HMD_N = NA * HMD_logmx
HMD_D = NA * HMD_logmx

for (j in colnames(HMD_logmx)) {
  
  this.schedule = j
  info          = unlist(strsplit(this.schedule,":"))
  this.pop      = info[1]
  this.period   = info[2]
  
  logmx.star = HMD_logmx[,j]
  
  # exposures
  
  tmp = read.table(paste0('../HMD-input-data/Exposures_1x10/',this.pop,
                          '.Exposures_1x10.txt'), 
                   stringsAsFactors = FALSE,
                   skip=2, header=TRUE) %>%
    filter(Year==this.period)
  
  tmp$Age = 0:110
  
  HMD_N[,j] = filter(tmp, Age < 100)$Female   # female exposure 0..99
  
  # deaths
  
  tmp = read.table(paste0('../HMD-input-data/Deaths_1x10/',this.pop,
                          '.Deaths_1x10.txt'), 
                   stringsAsFactors = FALSE,
                   skip=2, header=TRUE) %>%
    filter(Year==this.period)
  
  tmp$Age = 0:110
  
  HMD_D[,j] = filter(tmp, Age < 100)$Female   # female deaths 0..99
  
  
} # for j (HMD schedule)


##################################################################
# There are 3 D-spline models to test, each with a different
# set of residuals. The cubic spline model will always be 
#      log mortality at ages 0..99 = S * theta
# where S is a cubic B-spline matrix with lots of knots.
#
# Residual measures are
#   1. (1st diffs of S * theta) - (HMD avg 1st diffs) 
#   2. (2nd diffs of S * theta) - (HMD avg 2nd diffs) 
#   3. (projection residuals of [S * theta - Lee-Carter a] on
#           [Lee-Carter b]
#
# That is, the 3rd model should reward splines 
# that are Lee-Carter-ish 
#
# For all 3 measures the residuals have the form
#       eps = A*S*theta - c
# with difft A matrices and c vectors
#
# 1. A = 99x100 1st-diff matrix,  c = 99x1 HMD avg 1st diffs
# 2. A = 98x100 2nd-diff matrix,  c = 98x1 HMD avg 2nd diffs
# 3. A = 100x100 residual proj matrix M=(I-bb') using LC "b" vector,
#    c = 100x1 M * a 
#    This comes from M(S*theta-a) = 0 in an exact LC model
##################################################################

# construct a list with A,c,etc

Dspline_constants        = vector('list',3)
names(Dspline_constants) = c('diff1','diff2','LC')

#-------------------------------
# diff1 matrix and vector of
# HMD mean differences
#-------------------------------
D1 = diff(diag(100),diff=1)
c1 = rowMeans(D1 %*% HMD_logmx)

Dspline_constants[['diff1']] = list(A=D1, c=c1)

#-------------------------------
# diff2 matrix and vector of
# HMD mean differences
#-------------------------------
D2 = diff(diag(100),diff=2)
c2 = rowMeans(D2 %*% HMD_logmx)

Dspline_constants[['diff2']] = list(A=D2, c=c2)

#-------------------------------
# LC matrix and vector
#-------------------------------
a = rowMeans(HMD_logmx)

tmp = svd( sweep(HMD_logmx,1,a,'-') )
b   = tmp$u[,1]

Mb  = diag(length(b)) - tcrossprod(b)
Mba = as.vector(Mb %*% a)

Dspline_constants[['LC']] = list(A=Mb, c=Mba, LCa=a, LCb=b)


#################################################################
# calculate the HMD residuals A[empirical log mx] - c 
# for each penalty, and the (co)variances of those residuals in
# the HMD
#################################################################

for (this.model in names(Dspline_constants)) {

   this.A = Dspline_constants[[this.model]]$A
   this.c = Dspline_constants[[this.model]]$c
   
   eps = this.A %*% HMD_logmx - this.c
   
   SIGMA = var(t(eps))
   
   # GENERALIZED inverse avoids some numerical problems   
   SIGMA.INV = MASS::ginv(SIGMA)  
   
   Dspline_constants[[this.model]]$SIGMA.INV = SIGMA.INV
}

save(HMD_logmx, HMD_D, HMD_N, file='../data/HMD.RData')
save(Dspline_constants      , file='../data/Dspline_constants.RData')

