#----------------------------------------------
# Examine the Mort1DSmooth (and other) fits 
# from the big cross-validated experiment
# 
# NOTE: 
# 1. this program may take a LONG time to run
# 2. you MUST replace the value of time_stamp below with
#    the appropriate name for the file that YOU created
#    in STEP-2. It will look like YYYY-MM-DD-HHMM.
----------------------------------------------
library(tidyverse)

rm(list=ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load estimates and data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nparts = 3
time_stamp = '...your time stamp goes here...' # formatted like: 2020-07-26-1925  

fname = paste0('../data/BIG-cross-validated-HMD-test-',
               time_stamp,'-part-1.Rdata')

load(fname)
result = subresult

for (pt in 2:(nparts-1)) {
  fname = paste0('../data/BIG-cross-validated-HMD-test-',
                 time_stamp,'-part-',pt,'.Rdata')
  load(fname)
  result = bind_rows(result, subresult)
}

fname = paste0('../data/BIG-cross-validated-HMD-test-',
               time_stamp,'-part-',nparts,'.Rdata')

load(fname)

load('../data/HMD.Rdata')

# because the "ix" indices in the results now refer
# to only the limited set of HMD schedules used in the 
# analysis, we have to use the "crosswalk" to 
# grab the right columns from HMD_logmx as the true values
#
# For ex. if the first schedule actually used was the
# 6th column of HMD_logmx, then ix=1 actually refers to 
# HMD column 6. 

# re-order method and sample size to be consistent with
# what we want in the output tables (100K first, P-spline first)
df = left_join(result, crosswalk) %>% 
       left_join(true_values, by='ix') %>% 
       mutate( method = 
                 factor(method, 
                   levels=c('MS'      ,'D1','D2','LC'),
                   labels=c('P-spline','D-1','D-2','D-LC')),
               ssize = 
                 factor(sample_size,
                   levels=c(100000, 10000),
                   labels=c('Population = 100,000','Population = 10,000'))
       )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # check whether all estimators converged ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ok_table = 
  df %>% 
  group_by(ssize, method) %>% 
  summarize(pct_ok = 100*mean(!is.na(e0_hat)),
            n_ok = sum(!is.na(e0_hat)))

ok_table

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# much faster unnest() using
# data.table functions
# from https://www.r-bloggers.com/much-faster-unnesting-with-data-table/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rlang)
library(data.table)

unnest_dt <- function(tbl, col) {
   tbl   <- as.data.table(tbl)
   col   <- ensyms(col)
   clnms <- syms(setdiff(colnames(tbl), as.character(col)))
   tbl   <- as.data.table(tbl)
   tbl   <- eval(
    expr(tbl[, as.character(unlist(!!!col)), by = list(!!!clnms)])
   )
   colnames(tbl) <- c(as.character(clnms), as.character(col))
   tbl
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# expand by age ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long_df = df %>% 
  filter(is.finite(e0_hat)) %>%     # filter out the 85 non-converging cases
  dplyr::select(trial:logmx_hat,method,ssize)

# use faster unnesting
long_df = unnest_dt(long_df,logmx_hat) %>% 
          as_tibble() %>% 
          transform(age=0:99) %>% 
          mutate(logmx_hat = as.numeric(logmx_hat))
  
L = true_values %>% 
  dplyr::select(ix,logmx) %>% 
  unnest(cols='logmx') %>% 
  transform(age=0:99) 

long_df = left_join(long_df, L, by=c('ix','age')) %>% 
            mutate(err = logmx_hat - logmx)  
                            

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASMR errors ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mae_table = 
  long_df %>% 
  group_by(ssize,method) %>% 
  summarize( me   = mean(err),
             mae  = mean(abs(err))
             )

# small summary table             
mae_table %>% 
  pivot_longer(cols=c('mae','me')) %>% 
  arrange(ssize, name, method) %>% 
  pull(value) %>% 
  array(dim=c(4,2,2),
        dimnames=list(levels(mae_table$method),
                      c('MAE','ME'),
                      levels(mae_table$ssize))) %>% 
  aperm(perm=c(2,1,3)) %>% 
  round(2)

# plot mean errors by method and age
# (ugly code... :/)

#devtools::install_github("zeehio/facetscales")
library(facetscales)

scales_y <- list(
  `Population = 100,000` = scale_y_continuous(limits = c(-2, 1.6)),
  `Population = 10,000`  = scale_y_continuous(limits = c(-4.5, 1.7))
)

G = long_df %>% 
  group_by(ssize, method, age) %>% 
  summarize( me = mean(err),
             me10 = quantile(err,.10),
             me50 = quantile(err,.50),
             me90 = quantile(err,.90)
  ) %>% 
  ggplot() +
  aes(x=age, y=me50, color=method) +
  scale_x_continuous(breaks=seq(0,100,20)) +
  geom_line(size=1.5) +
  geom_ribbon(aes(x=age,ymin=me10,ymax=me90, fill=method),
              color=NA,alpha=.25) +
  facet_grid_sc(ssize~method, scales=list(y=scales_y)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values=c('red','darkgreen','blue','brown')) +
  scale_fill_manual(values=c('red','darkgreen','blue','brown')) +
  guides(color=FALSE, fill=FALSE) +
  labs(y="Median Error and 80% interval") +
  theme_bw() +
  theme(strip.background = element_rect(fill=grey(.90)),
        strip.text=element_text(face='bold',size=15))

ggsave(filename='../plots/mean-errors-by-method-and-age.pdf',
       plot=G, width=11, height=8.5)

ggsave(filename='../plots/mean-errors-by-method-and-age.eps',
       device=cairo_ps,
       plot=G, width=11, height=8.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BIC and effective df info ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bic_table =
  df %>% 
  filter(is.finite(e0_hat)) %>% 
  group_by(ssize, method) %>% 
  summarize( median_bic = median(bic),
             median_df = median(df)) %>% 
  dplyr::select(ssize,method,median_bic,median_df) %>% 
  data.frame() %>% 
  mutate_if(is.numeric,round,1)

# small summary table             
bic_table %>% 
  pivot_longer(cols=c('median_bic','median_df')) %>% 
  arrange(ssize, name, method) %>% 
  pull(value) %>% 
  array(dim=c(4,2,2),
        dimnames=list(levels(bic_table$method),
                      c('BIC','DF'),
                      levels(bic_table$ssize))) %>% 
  aperm(perm=c(2,1,3)) %>% 
  round(2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 'winners' for total error
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print('Winners - Sum of Absolute Errors (all ages)')

W = long_df %>% 
      group_by(ssize,ix,trial,method) %>% 
      summarize( e  = sum(abs(err))) %>% 
      group_by(ssize,ix,trial) %>% 
      summarize( winner = method[which.min(e)]) 

table(W$ssize, W$winner) %>% 
  prop.table(margin=1) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BIC by schedule
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print('Winners - BIC')

W = df %>% 
  group_by(ssize,ix,trial) %>% 
  summarize( winner = method[which.min(bic)]) 

table(W$ssize, W$winner) %>% 
  prop.table(margin=1) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# e0 errors ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

e0_table = 
  df %>% 
  filter(is.finite(e0_hat)) %>% 
  mutate(err = e0_hat - e0) %>% 
  group_by(ssize,method) %>%
  summarize( me   = mean(err),
             mae  = mean(abs(err)))

e0_table

# small summary table             
e0_table %>% 
  pivot_longer(cols=c('mae','me')) %>% 
  arrange(ssize, name, method) %>% 
  pull(value) %>% 
  array(dim=c(4,2,2),
        dimnames=list(levels(e0_table$method),
                      c('MAE','ME'),
                      levels(e0_table$ssize))) %>% 
  aperm(perm=c(2,1,3)) %>% 
  round(2)

# e0 winners
print('Winners - e0')

W = df %>% 
  mutate( e = abs(e0_hat - e0)) %>% 
  group_by(ssize,ix,trial) %>% 
  summarize( winner = method[which.min(e)]) 

table(W$ssize, W$winner) %>% 
  prop.table(margin=1) 



# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # e60 errors ----
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


e60_table = 
  df %>% 
  filter(is.finite(e60_hat)) %>% 
  mutate(err = e60_hat - e60) %>% 
  group_by(ssize,method) %>%
  summarize( me   = mean(err),
             mae  = mean(abs(err)))

e60_table

# small summary table             
e60_table %>% 
  pivot_longer(cols=c('mae','me')) %>% 
  arrange(ssize, name, method) %>% 
  pull(value) %>% 
  array(dim=c(4,2,2),
        dimnames=list(levels(e60_table$method),
                      c('MAE','ME'),
                      levels(e60_table$ssize))) %>% 
  aperm(perm=c(2,1,3)) %>% 
  round(2)

# e0 winners
print('Winners - e60')

W = df %>% 
  mutate( e = abs(e60_hat - e60)) %>% 
  group_by(ssize,ix,trial) %>% 
  summarize( winner = method[which.min(e)]) 

table(W$ssize, W$winner) %>% 
  prop.table(margin=1) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 45q20 experiments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# true values
nqx = function(logmx, x=20, n=45) {
  H = cumsum(c(0, exp(logmx)))
  lx = exp(-H)
  10000*(1-lx[x+1+n]/lx[x+1])  
}

this_x = 20
this_n = 45

true_values = true_values %>% 
               mutate(q = map_dbl(logmx,nqx,x=this_x,n=this_n))



tmp = df %>% 
        filter(is.finite(e0_hat)) %>% 
        mutate(qhat = map_dbl(logmx_hat, nqx, x=this_x, n=this_n)) %>% 
        dplyr::select(trial,ix,method,ssize, qhat) %>% 
        left_join(dplyr::select(true_values,ix,q), by='ix')

q_table = tmp %>% 
     mutate(err=qhat-q) %>% 
     group_by(ssize, method) %>% 
     summarize( me = mean(err),
                mae=mean(abs(err)),
                mape = mean(100*abs(err/q)))

q_table


# small summary table             
q_table %>% 
  pivot_longer(cols=c('mae','me','mape')) %>% 
  arrange(ssize, name, method) %>% 
  pull(value) %>% 
  array(dim=c(4,3,2),
        dimnames=list(levels(q_table$method),
                      c('MAE','MAPE','ME'),
                      levels(q_table$ssize))) %>% 
  aperm(perm=c(2,1,3)) %>% 
  round(1)


print('Winners - 45q20')

W = tmp %>%
  mutate( e = abs(qhat-q)) %>% 
  group_by(ssize,ix,trial) %>% 
  summarize( winner = method[which.min(e)]) 

table(W$ssize, W$winner) %>% 
  prop.table(margin=1) 
