#-------------------------------------------------------------------
# creates 8-panel plot illustrating the P-spline, D1, D2, DLC fits
# to a single simulated data set using Portugal 1970-1979 female 
# log mortality rates
# 
# includes pointwise confidence intervals
#-------------------------------------------------------------------

library(tidyverse)
library(MortalitySmooth)

rm(list=ls())

# use some utility functions (code is in other files)
source('vech-xpnd-functions.R')
source('Dspline_fit_function.R')

#load Dspline_constants ( +a few other dfs) from previous analysis ----

time_stamp = '2020-07-26-1925'

fname = paste0('../data/BIG-cross-validated-HMD-test-',
               time_stamp,'-part-3.Rdata')

load(fname)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# hard-coded data ---- 
#' These are simulated samples from the Portugal 1970-1979
#' HMD population age structure and mortality rates

true_mx =
  c(`0` = 0.03591, `1` = 0.00406, `2` = 0.00189, `3` = 0.00131,
    `4` = 0.00094, `5` = 0.00074, `6` = 0.00067, `7` = 0.00053, `8` = 0.00047,
    `9` = 5e-04,  `10` = 4e-04,  `11` = 0.00045, `12` = 4e-04, `13` = 0.00043,
    `14` = 0.00048, `15` = 0.00047, `16` = 0.00053, `17` = 0.00055,
    `18` = 0.00059, `19` = 6e-04, `20` = 6e-04, `21` = 0.00069, `22` = 0.00063,
    `23` = 0.00068, `24` = 0.00073, `25` = 0.00072, `26` = 0.00078,
    `27` = 0.00081, `28` = 0.00082, `29` = 0.00083, `30` = 0.00089,
    `31` = 0.00098, `32` = 0.00105, `33` = 0.00109, `34` = 0.00115,
    `35` = 0.00129, `36` = 0.00136, `37` = 0.0014,  `38` = 0.00167,
    `39` = 0.00177, `40` = 0.00187, `41` = 0.0019,  `42` = 0.0021,
    `43` = 0.00234, `44` = 0.00242, `45` = 0.00261, `46` = 0.00297,
    `47` = 0.00316, `48` = 0.00357, `49` = 0.00368, `50` = 0.00416,
    `51` = 0.00425, `52` = 0.00474, `53` = 0.00517, `54` = 0.00538,
    `55` = 0.00593, `56` = 0.00654, `57` = 0.00684, `58` = 0.0079,
    `59` = 0.00839, `60` = 0.00962, `61` = 0.01041, `62` = 0.01173,
    `63` = 0.01309, `64` = 0.01459, `65` = 0.01586, `66` = 0.01788,
    `67` = 0.02012, `68` = 0.02289, `69` = 0.02609, `70` = 0.02989,
    `71` = 0.03326, `72` = 0.03805, `73` = 0.04374, `74` = 0.05071,
    `75` = 0.05811, `76` = 0.06464, `77` = 0.07137, `78` = 0.08402,
    `79` = 0.09513, `80` = 0.09406, `81` = 0.10159, `82` = 0.11712,
    `83` = 0.12813, `84` = 0.1441,  `85` = 0.15551, `86` = 0.17847,
    `87` = 0.19435, `88` = 0.21617, `89` = 0.23707, `90` = 0.25522,
    `91` = 0.26763, `92` = 0.28701, `93` = 0.31319, `94` = 0.34524,
    `95` = 0.37964, `96` = 0.4079,  `97` = 0.43677, `98` = 0.46609,
    `99` = 0.49563)

N10K =
  c(195L, 208L, 146L, 180L, 158L, 185L, 178L, 179L, 176L, 179L,
    154L, 154L, 167L, 182L, 176L, 197L, 195L, 166L, 152L, 180L, 169L,
    140L, 158L, 139L, 150L, 141L, 110L, 121L, 119L, 116L, 122L, 124L,
    124L, 118L, 126L, 109L, 141L, 107L, 110L, 113L, 122L, 107L, 128L,
    120L, 126L, 120L, 133L, 126L, 121L, 139L, 128L, 109L, 84L, 109L,
    103L, 100L, 94L, 97L, 97L, 84L, 124L, 97L, 114L, 94L, 90L, 103L,
    85L, 82L, 74L, 82L, 70L, 79L, 57L, 55L, 57L, 49L, 48L, 47L, 40L,
    38L, 35L, 18L, 32L, 21L, 19L, 19L, 21L, 14L, 10L, 2L, 1L, 4L,
    0L, 2L, 1L, 3L, 0L, 1L, 1L, 0L)

D10K =
  c(9L, 3L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
    1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L,
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
    0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 2L,
    0L, 1L, 0L, 2L, 0L, 3L, 0L, 1L, 5L, 1L, 1L, 2L, 1L, 1L, 5L, 8L,
    3L, 2L, 1L, 4L, 1L, 7L, 3L, 1L, 2L, 1L, 1L, 0L, 1L, 0L, 1L, 0L,
    2L, 0L, 0L, 0L, 0L)

N100K =
  c(1718L, 1652L, 1584L, 1762L, 1691L, 1692L, 1715L, 1705L, 1654L,
    1747L, 1716L, 1783L, 1665L, 1725L, 1667L, 1672L, 1655L, 1604L,
    1602L, 1599L, 1487L, 1568L, 1492L, 1546L, 1405L, 1420L, 1340L,
    1258L, 1260L, 1244L, 1201L, 1203L, 1197L, 1172L, 1137L, 1209L,
    1234L, 1258L, 1324L, 1219L, 1281L, 1312L, 1244L, 1277L, 1269L,
    1254L, 1205L, 1218L, 1179L, 1200L, 1204L, 1150L, 1148L, 1100L,
    1035L, 1030L, 1083L, 1023L, 993L, 1037L, 991L, 953L, 951L, 904L,
    917L, 887L, 894L, 864L, 827L, 835L, 791L, 717L, 738L, 634L, 612L,
    556L, 473L, 452L, 388L, 331L, 390L, 317L, 284L, 221L, 216L, 177L,
    163L, 119L, 101L, 79L, 54L, 29L, 36L, 26L, 21L, 12L, 6L, 5L,
    2L, 3L)

D100K =
  c(70L, 5L, 4L, 3L, 1L, 2L, 1L, 1L, 2L, 0L, 0L, 2L, 0L, 0L, 1L,
    1L, 0L, 2L, 2L, 0L, 0L, 0L, 1L, 0L, 0L, 2L, 0L, 0L, 0L, 1L, 0L,
    2L, 0L, 1L, 0L, 1L, 1L, 6L, 3L, 2L, 6L, 4L, 3L, 5L, 2L, 2L, 4L,
    3L, 8L, 5L, 7L, 5L, 3L, 5L, 6L, 7L, 8L, 8L, 6L, 9L, 11L, 9L,
    11L, 16L, 16L, 11L, 19L, 18L, 16L, 32L, 26L, 22L, 34L, 21L, 32L,
    32L, 27L, 30L, 37L, 26L, 34L, 31L, 37L, 33L, 37L, 22L, 25L, 26L,
    17L, 21L, 9L, 10L, 16L, 11L, 6L, 1L, 4L, 2L, 3L, 2L)

# get Dspline constants for PRT ----
tmp = filter(Dspline_constants, country=='PRT')

D1 = diff( diag(100), diff=1)
D2 = diff( diag(100), diff=2)

A1 = D1
c1 = unlist( tmp$c1)
V1 = xpnd( unlist(tmp$V1))

A2 = D2
c2 = unlist( tmp$c2)
V2 = xpnd( unlist(tmp$V2))

A3 = xpnd( unlist(tmp$A3))
c3 = unlist( tmp$c3)
V3 = xpnd( unlist(tmp$V3))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# n=10,000 example first ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N = N10K
D = D10K

zero = (D==0 & N>0)
good = (D>0)
age  = 0:99

raw_data = tibble(age,true_mx,good,zero,N,D)

w = 1*(N>0) # disregard ages with zero exposure

# P-spline curve for this data, min-BIC

MS_default = Mort1Dsmooth( x      = 0:99, 
                  y      = D,
                  w      = w,
                  deg    = 3,
                  offset = log(N),
                  ndx    = 33,
                  method = 1,
                  control = list(TOL2=.20,
                                 RANGE=c(1,1e6)))

# D-spline fits for this data

fit_D1 = Dspline_fit(N,D,
                  Amatrix=A1, 
                  cvector=c1, 
                  SIGMA.INV = V1,
                  use_MS_basis = TRUE,
                  details=TRUE,
                  max_iter = 100,
                  theta_tol = .005)

fit_D2 = Dspline_fit(N,D,
                     Amatrix=A2, 
                     cvector=c2, 
                     SIGMA.INV = V2,
                     use_MS_basis = TRUE,
                     details=TRUE,
                     max_iter = 100,
                     theta_tol = .005)

fit_LC = Dspline_fit(N,D,
                     Amatrix=A3, 
                     cvector=c3, 
                     SIGMA.INV = V3,
                     use_MS_basis = TRUE,
                     details=TRUE,
                     max_iter = 100,
                     theta_tol = .005)

#=====================================================================
# for each fit to this one sample, use the inverse of the negative 
# Hessian to calculate the pointwise median, 2.5%ile, and 97.5%ile
# of the fitted spline values at ages 0..99
#=====================================================================

pointwise_CI = function(this_basis, this_theta, this_Hessian,
                        p = c(.025, .975)) 
{
  spline_fit   = this_basis %*% this_theta
  covar_theta  = solve(-this_Hessian)
  covar_spline = this_basis %*% covar_theta %*% t(this_basis) 
  sdev       = sqrt( diag( covar_spline))
  return( cbind( as.vector(spline_fit + sdev * qnorm(p[1])),
                 as.vector(spline_fit + sdev * qnorm(p[2]))) )
}

# need a separate calculation to recover the hessian 
# for the P-spline (=MortalitySmooth default) object
MS_hessian = function(this_MS_fit) {
  K = length(this_MS_fit$coefficients)
  
  Dhat = this_MS_fit$fitted.values
  Dhat[this_MS_fit$w==0] = 0

  BWB = t(this_MS_fit$B) %*% diag(Dhat) %*% this_MS_fit$B
  LP  = this_MS_fit$lambda * crossprod( diff(diag(K), diff=2))
  
  hessian = -(BWB + LP)

  return(hessian)
}

##################################################
# calculate the confidence intervals for each fit
##################################################

QMS = pointwise_CI(MS_default$B, 
                   MS_default$coefficients, 
                   MS_hessian(MS_default))
Q1  = pointwise_CI(fit_D1$B, 
                   fit_D1$theta, 
                   fit_D1$hessian)
Q2  = pointwise_CI(fit_D2$B, 
                   fit_D2$theta, 
                   fit_D2$hessian)
QLC  = pointwise_CI(fit_LC$B, 
                   fit_LC$theta, 
                   fit_LC$hessian)

### plot data for the n=10,000 case

plot_data1 = expand.grid(age=0:99, 
                  method = c('P-spline','D-1','D-2','D-LC')) %>% 
      add_column(sample_size = 'Population = 10,000',
                 N = rep(N,4),
                 D = rep(D,4),
                 zero = rep(zero,4),
                 good = rep(good,4),
                 fit = c( exp(MS_default$logmortality),
                          exp(fit_D1$lambda.hat),
                          exp(fit_D2$lambda.hat),
                          exp(fit_LC$lambda.hat)),
                 fitL = c(QMS[,1],
                           Q1[,1],
                           Q2[,1],
                           QLC[,1]),
                 fitH = c(QMS[,2],
                           Q1[,2],
                           Q2[,2],
                           QLC[,2])
                 )
      


annotate_data1 = expand.grid(method = c('P-spline','D-1','D-2','D-LC'),
                             sample_size = 'Population = 10,000') %>% 
  add_column( df  =    c(MS_default$df,  fit_D1$df,  fit_D2$df,  fit_LC$df),
              bic =    c(MS_default$bic, fit_D1$bic, fit_D2$bic, fit_LC$bic),
              lambda = c(MS_default$lambda, NA,NA,NA))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# N = 100,000 case second ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = N100K
D = D100K

zero = (D==0 & N>0)
good = (D>0)
age  = 0:99

raw_data = tibble(age,true_mx,good,zero,N,D)

w = 1*(N>0) # disregard ages with zero exposure

# P-spline curve for this data, min-BIC

MS_default = Mort1Dsmooth( x      = 0:99, 
                           y      = D,
                           w      = w,
                           deg    = 3,
                           offset = log(N),
                           ndx    = 33,
                           method = 1,
                           control = list(TOL2=.20,
                                          RANGE=c(1,1e6)))

# D-spline fits for this data

fit_D1 = Dspline_fit(N,D,
                     Amatrix=A1, 
                     cvector=c1, 
                     SIGMA.INV = V1,
                     use_MS_basis = TRUE,
                     details=TRUE,
                     max_iter = 100,
                     theta_tol = .005)

fit_D2 = Dspline_fit(N,D,
                     Amatrix=A2, 
                     cvector=c2, 
                     SIGMA.INV = V2,
                     use_MS_basis = TRUE,
                     details=TRUE,
                     max_iter = 100,
                     theta_tol = .005)

fit_LC = Dspline_fit(N,D,
                     Amatrix=A3, 
                     cvector=c3, 
                     SIGMA.INV = V3,
                     use_MS_basis = TRUE,
                     details=TRUE,
                     max_iter = 100,
                     theta_tol = .005)



##################################################
# calculate the confidence intervals for each fit
##################################################

QMS = pointwise_CI(MS_default$B, 
                   MS_default$coefficients, 
                   MS_hessian(MS_default))
Q1  = pointwise_CI(fit_D1$B, 
                   fit_D1$theta, 
                   fit_D1$hessian)
Q2  = pointwise_CI(fit_D2$B, 
                   fit_D2$theta, 
                   fit_D2$hessian)
QLC  = pointwise_CI(fit_LC$B, 
                    fit_LC$theta, 
                    fit_LC$hessian)


  ### plot data for the n=100,000 case
  
  plot_data2 = expand.grid(age=0:99, 
                           method = c('P-spline','D-1','D-2','D-LC')) %>% 
  add_column(sample_size = 'Population = 100,000',
             N = rep(N,4),
             D = rep(D,4),
             zero = rep(zero,4),
             good = rep(good,4),
             fit = c( exp(MS_default$logmortality),
                      exp(fit_D1$lambda.hat),
                      exp(fit_D2$lambda.hat),
                      exp(fit_LC$lambda.hat)
                      ),
             fitL = c(QMS[,1],
                       Q1[,1],
                       Q2[,1],
                       QLC[,1]),
             fitH = c(QMS[,2],
                       Q1[,2],
                       Q2[,2],
                       QLC[,2])
  )


  annotate_data2 = expand.grid(method = c('P-spline','D-1','D-2','D-LC'),
                               sample_size = 'Population = 100,000') %>% 
    add_column( df  = c(MS_default$df, fit_D1$df, fit_D2$df, fit_LC$df),
                bic = c(MS_default$bic, fit_D1$bic, fit_D2$bic, fit_LC$bic),
                lambda = c(MS_default$lambda, NA,NA,NA))
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make multipanel plot ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_data = rbind(plot_data1, 
                  plot_data2) %>% 
            mutate(ssize = factor(sample_size,
                                 levels=c('Population = 100,000',
                                          'Population = 10,000'))) %>% 
            select(-sample_size)
  
annotate_data = rbind(annotate_data1, 
                      annotate_data2) %>% 
    mutate(ssize = factor(sample_size,
                          levels=c('Population = 100,000',
                                   'Population = 10,000')),
           tx = paste0(ifelse(is.na(lambda),'',
                          paste0('Lambda = ',sprintf('%.2f',lambda))),
                       '\nBIC = ',sprintf('%.1f',bic),
                       '\n df = ',sprintf('%0.1f',df))) %>% 
    select(-sample_size)     
  
#---------------------------------------------------------
# big multi-panel plot
#---------------------------------------------------------  

G = ggplot() +
  geom_line(data=plot_data,aes(x=age, y=fit,color=method),
            size=2) +
  geom_ribbon(data=plot_data,aes(x=age,ymin=exp(fitL),ymax=exp(fitH),
                                 fill=method), alpha=.25) +
  geom_line(data=raw_data,aes(x=age, y=true_mx), color='black', 
            size=0.8,lty='solid') +
  geom_point(data=filter(plot_data,good),
             aes(x=age,y=D/N), shape='+',size=2.5) +
  geom_rug(data=filter(plot_data,zero),aes(x=age)) +
  labs(x='Age',y='Mortality Rate') +
  scale_y_log10(labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(breaks=seq(0,100,10),minor_breaks = NULL) +
  scale_color_manual(values=c('red','darkgreen','blue','brown')) +
  scale_fill_manual(values=c('red','darkgreen','blue','brown')) +
  guides(color=FALSE,fill=FALSE) +
  theme_bw() +
  theme(strip.background = element_rect(fill=grey(.90)),
        strip.text       = element_text(face='bold',size=14)) +
  facet_grid(ssize~method)


  G = G + 
  geom_text(data=annotate_data, x=5, y=log(.85), 
            aes(label= annotate_data$tx), size=5, hjust=0)
G


ggsave('../plots/four-methods-Portugal-data.pdf',
       plot=G,width=11, height=8.5)

ggsave('../plots/four-methods-Portugal-data.eps',
       device=cairo_ps,
       plot=G,width=11, height=8.5)
