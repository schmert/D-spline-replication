#----------------------------------------------------------------
# creates multipanel plot with alternative P-spline fits
# to simulated small-area data using Portugal 1970-1979
# female log mortality rates
# 
# six panels: 2 sample sizes x 3 lambdas 
#----------------------------------------------------------------

library(tidyverse)
library(MortalitySmooth)

rm(list=ls())

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

MS_default1 = Mort1Dsmooth( x      = 0:99, 
                  y      = D,
                  w      = w,
                  deg    = 3,
                  offset = log(N),
                  ndx    = 33,
                  method = 1,
                  control=list(TOL2=.20,
                               RANGE=c(1,1e6)))

# find the approx deg of freedom for the minimum-BIC estimator
df_star = MS_default1$df

# rerun with 1/2 and 2x as many dof

MS_half =   Mort1Dsmooth( x      = 0:99, 
                          y      = D,
                          w      = w,
                          deg    = 3,
                          offset = log(N),
                          ndx    = 33,
                          method = 4,
                          df = 0.5*df_star,
                          control=list(TOL2=.20))

MS_double =   Mort1Dsmooth( x      = 0:99, 
                          y      = D,
                          w      = w,
                          deg    = 3,
                          offset = log(N),
                          ndx    = 33,
                          method = 4,
                          df = 2*df_star,
                          control=list(TOL2=.20))

### plot data for the n=10,000 case

plot_data1 = expand.grid(age=0:99, 
                  smoothing_level=c('Lower Penalty',
                                    'Minimum-BIC Penalty',
                                    'Higher Penalty')) %>% 
      add_column(sample_size = 'Population = 10,000',
                 N = rep(N,3),
                 D = rep(D,3),
                 zero = rep(zero,3),
                 good = rep(good,3),
                 fit = c( exp(MS_double$logmortality),
                          exp(MS_default1$logmortality),
                          exp(MS_half$logmortality)))

annotate_data1 = expand.grid(smoothing_level = c('Lower Penalty',
                                                 'Minimum-BIC Penalty',
                                                 'Higher Penalty'),
                             sample_size = 'Population = 10,000') %>% 
                add_column( df  =    c(MS_double$df,     MS_default1$df,     MS_half$df),
                            bic =    c(MS_double$bic,    MS_default1$bic,    MS_half$bic),
                            lambda = c(MS_double$lambda, MS_default1$lambda, MS_half$lambda))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# n=100,000 example second ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N = N100K
D = D100K

zero = (D==0 & N>0)
good = (D>0)
age  = 0:99

raw_data = tibble(age,true_mx,good,zero,N,D)

w = 1*(N>0) # disregard ages with zero exposure


MS_default2 = Mort1Dsmooth( x      = 0:99, 
                           y      = D,
                           w      = w,
                           deg    = 3,
                           offset = log(N),
                           ndx    = 33,
                           method = 1,
                           control=list(TOL2=.20,
                                        RANGE=c(1,1e6)))

df_star = MS_default2$df

# rerun with 1/2 and 2x as many dof

MS_half =   Mort1Dsmooth( x      = 0:99, 
                          y      = D,
                          w      = w,
                          deg    = 3,
                          offset = log(N),
                          ndx    = 33,
                          method = 4,
                          df = 0.5*df_star,
                          control=list(TOL2=.20))

MS_double =   Mort1Dsmooth( x      = 0:99, 
                            y      = D,
                            w      = w,
                            deg    = 3,
                            offset = log(N),
                            ndx    = 33,
                            method = 4,
                            df = 2*df_star,
                            control=list(TOL2=.20))

### plot data for the n=100,000 case

plot_data2 = expand.grid(age=0:99, 
                         smoothing_level=c('Lower Penalty',
                                           'Minimum-BIC Penalty',
                                           'Higher Penalty')) %>% 
  add_column(sample_size = 'Population = 100,000',
             N = rep(N,3),
             D = rep(D,3),
             zero = rep(zero,3),
             good = rep(good,3),
             fit = c( exp(MS_double$logmortality),
                      exp(MS_default2$logmortality),
                      exp(MS_half$logmortality)))


annotate_data2 = expand.grid(smoothing_level = c('Lower Penalty',
                                                'Minimum-BIC Penalty',
                                                 'Higher Penalty'),
                             sample_size = 'Population = 100,000') %>% 
  add_column( df  = c(MS_double$df, MS_default2$df, MS_half$df),
              bic = c(MS_double$bic, MS_default2$bic, MS_half$bic),
              lambda = c(MS_double$lambda, MS_default2$lambda, MS_half$lambda))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make multipanel plot ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_data = rbind(plot_data1, 
                  plot_data2) %>% 
            mutate(ssize = factor(sample_size,
                                 levels=c('Population = 100,000',
                                          'Population = 10,000')))
annotate_data = rbind(annotate_data1, 
                      annotate_data2) %>% 
  mutate(ssize = factor(sample_size,
                        levels=c('Population = 100,000',
                                 'Population = 10,000')),
         tx = paste0('Lambda = ',sprintf('%.2f',lambda),
                     '\nBIC = ',sprintf('%.1f',bic),
                     '\ndf = ',sprintf('%0.1f',df)))


G = ggplot() +
  scale_y_log10(labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(breaks=seq(0,100,10)) +
  geom_line(data=raw_data,aes(x=age, y=true_mx), color='black', 
            size=3,alpha=.40) +
  geom_line(data=plot_data,aes(x=age, y=fit,color=ssize),
            size=2) +
  scale_color_manual(values=c('red','royalblue')) +
  guides(color=FALSE) +
  theme_bw() +
  labs(x='Age',y='Mortality Rate') +
  theme(strip.background = element_rect(fill=grey(.90)),
        strip.text = element_text(face='bold',size=14)) +
  geom_point(data=filter(plot_data,good),
             aes(x=age,y=D/N),
             shape='+',size=4) +
  geom_rug(data=filter(plot_data,zero),
           aes(x=age)) +
  facet_grid(ssize~smoothing_level)

G = G + 
      geom_text(data=annotate_data, x=5, y=log(.85), 
                aes(label= annotate_data$tx), size=5, hjust=0)
G

ggsave('../plots/alternative-P-splines-Portugal-data-different-lambdas.pdf',
       plot=G,width=11, height=8.5)

ggsave('../plots/alternative-P-splines-Portugal-data-different-lambdas.eps',
       device=cairo_ps,
       plot=G,width=11, height=8.5)
