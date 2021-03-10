#-------------------------------------------------------------
# creates a four-panel figure illustrating a cubic spline
# fit to Portugal the 1970-1990 female logmx schedule
# 
# Panels are
# (A) B-spline basis for the fits
# (B) Portugal 1970-1979 single-yr female rate schedule
# (C) A spline fit to PRT data, with cubic subsections
# (D) Zoom in on ages 0-15 and show cubics
#-------------------------------------------------------------

library(splines)
library(tidyverse)
library(gridExtra)

rm(list=ls())

theme_carl <- function () { 
  theme_bw() %+replace% 
    theme(axis.text = element_text(size=14, face='bold',color='black'),
          axis.title = element_text(size=14, face='bold'))
}



xcoarse  = seq(0,99,1)
xfine    = seq(0,99, .25)

Bfine    = splines::bs(xfine,   knots=seq(3,96,3), degree=3, intercept=TRUE)
Bcoarse  = splines::bs(xcoarse, knots=seq(3,96,3), degree=3, intercept=TRUE)


#------------------------------------------------------------------------
# 1: PRT raw data ----
#------------------------------------------------------------------------

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


logmx = log(true_mx)


#------------------------------------------------------------------------
# A: B-spline basis for the fits  ----
#------------------------------------------------------------------------

my_colors = c('#e41a1c','#377eb8','#4daf4a','#984ea3')

basis_df = expand.grid(x=xfine, fn=factor(seq(ncol(Bfine)))) %>%
  as_tibble() %>% 
  transform( y = as.vector(Bfine)) %>% 
  filter( y > 0)

GA = ggplot(data=basis_df) +
  aes(x=x, y=y, group=fn, color=fn) +
  geom_line(size=0.75) +
  labs(x='x', y='B(x)') +
  theme_carl() +
  guides(color=FALSE) +
  scale_color_manual(values=rep(my_colors,9)) +
  geom_text(aes(x=6,y=0.94), label='A.', fontface='bold',
            color='black',size=6, hjust=0) 


GA


#------------------------------------------------------------------------
# B: log mortality rates ----
#------------------------------------------------------------------------

GB = ggplot() +
      geom_point(aes(x=xcoarse, y=logmx), shape=1, size=1) +
      theme_carl() +
      scale_y_continuous(limits=c(-10,0), breaks=seq(-10,0,2)) +
      labs(x='x', y='Log Mortality Rate') +
      geom_text(aes(x=5,y=-1), label='B.', 
                fontface='bold',color='black',size=6, hjust=0)

GB

#------------------------------------------------------------------------
# C: spline fit to PRT data  ----
#------------------------------------------------------------------------

theta = solve( crossprod(Bcoarse), crossprod(Bcoarse,logmx))

piecewise_df = tibble(x    = xfine,
                      y    = as.vector(Bfine %*% theta),
                      xcat = cut(xfine, seq(0,99,3), right=FALSE))

# tweak last xcat so that 99 is included, then add hue and width
piecewise_df$xcat[length(xfine)] = piecewise_df$xcat[length(xfine)-1]

piecewise_df = piecewise_df %>% 
  mutate( hue  = 'navy',
          width= ifelse(hue==my_colors[1],0.6,0.6)
  )

GC = ggplot() +
  geom_point(aes(x=xcoarse, y=logmx), shape=1, size=1) +
  theme_carl() +
  labs(x='x', y='Log Mortality Rate') +
  scale_y_continuous(limits=c(-10,0), breaks=seq(-10,0,2)) +
  geom_vline(xintercept = seq(0,99,3), color='lightgrey') +
  geom_text(aes(x=5,y=-1.2), label='C.', fontface='bold',
            color='black',size=6,hjust=0) 

GC = GC + 
  geom_line(data=piecewise_df, 
            aes(x=x,y=y),
            color=piecewise_df$hue, 
            size=piecewise_df$width)

GC

#------------------------------------------------------------------------
# D: zoom in on spline fit at ages < 15  ----
#------------------------------------------------------------------------

# zoom in on first few spline segments

nseg = 4
xnew = seq(0,3*(nseg+1),.10)

this_linetype  = c('solid','dashed','dotted','solid')
this_linewidth = c(6, 3, 6, 3)

tmp = NULL

for (this_seg in 1:nseg) {
  valid_x = (xfine >= 3*(this_seg-1)) & (xfine <= 3*this_seg)
  xx = (piecewise_df %>% filter(valid_x) %>% pull(x))
  yy = (piecewise_df %>% filter(valid_x) %>% pull(y))
  
  fit = lm(yy ~ xx + I(xx^2) + I(xx^3))
  beta = coef(fit)
  
  ynew = beta[1] + beta[2] * xnew + beta[3] * xnew^2 + beta[4] * xnew^3
  
  tmp = rbind(tmp,
              tibble(x=xnew, y=ynew, seg=this_seg, 
                     linetype = this_linetype[this_seg],
                     linewidth = this_linewidth[this_seg],
                     linecolor = my_colors[1+ this_seg %% nseg],
                     active = (xnew >= 3*(this_seg-1)) & (xnew <= 3*this_seg)
              ))
}

tmp$seg = factor(tmp$seg)


clean_df = tmp %>% 
            filter(y > -10, y < -3) %>% 
            filter( !(seg==3 & x <1.4))

GD = ggplot(data=clean_df) + 
  geom_line(size=0.8,aes(x=x,y=y,group=seg,lty=linetype, color=linecolor)) +
  theme_carl() +
  labs(x='x',y='Log Mortality Rate') +
  scale_x_continuous(breaks=seq(0,15,3),minor_breaks = NULL) +
  guides(color=FALSE, lty=FALSE) +
  geom_text(aes(x=1,y=-1.2), label='D.', fontface='bold',
            color='black',size=6, hjust=0) +
  scale_color_identity() +
  scale_linetype_identity() +
  theme(panel.grid.major.x = element_line(color = grey(.25)))  +
  scale_y_continuous(limits=c(-10,0), breaks=seq(-10,0,2))


for (this_seg in 1:nseg) {
  tmp = filter(clean_df,seg==this_seg,active)
  this_linewidth = tmp$linewidth[1]
  GD = GD + 
    geom_line(data=tmp, 
              aes(x=x,y=y,group=seg,color=linecolor),
              size=this_linewidth, alpha=.70)
}

## this is a ridiculous way to add points, but ggplot is
## extremely finicky with the aes() othewise

GD = GD +
  geom_point(aes(x=  0, y=logmx[1])) +
  geom_point(aes(x=  1, y=logmx[2])) +
  geom_point(aes(x=  2, y=logmx[3])) +
  geom_point(aes(x=  3, y=logmx[4])) +
  geom_point(aes(x=  4, y=logmx[5])) +
  geom_point(aes(x=  5, y=logmx[6])) +
  geom_point(aes(x=  6, y=logmx[7])) +
  geom_point(aes(x=  7, y=logmx[8])) +
  geom_point(aes(x=  8, y=logmx[9])) +
  geom_point(aes(x=  9, y=logmx[10])) +
  geom_point(aes(x= 10, y=logmx[11])) +
  geom_point(aes(x= 11, y=logmx[12])) +
  geom_point(aes(x= 12, y=logmx[13])) 
  
  
GD


## arrange in a 2x2 grid and save in two formats
## (.pdf and .eps)
  
G = grid.arrange(GA,GC,
                 GB,GD, nrow=2)

ggsave('../plots/Figure-spline-example-Portugal-data.pdf',
      plot=G,height=7, width=10)

ggsave('../plots/Figure-spline-example-Portugal-data.eps',
       device=cairo_ps,
       plot=G,height=7, width=10)