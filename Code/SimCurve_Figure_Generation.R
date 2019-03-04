
load("Simulation curve figure/Figure_Generation_env.Rdata")

library(cowplot)
library(tidyverse)
library(gridExtra)
library(Rcpp)
library(pec)
library(randomForestSRC)
library(magrittr)

source("../ORSF benchmark/Source/oblique_survival_forest_predict.R")
source("../ORSF benchmark/Source/oblique_survival_forest_predictSurvProb.R")
source("../ORSF benchmark/Source/oblique_survival_forest_print.R")
source("../ORSF benchmark/Source/oblique_survival_forest_fit.R")
sourceCpp("../ORSF benchmark/Source/ORSF.cpp")
# ORSF source code needs to be sourced.

induce_mono <- function(x){
  for(i in 2:length(x)){
    if(x[i]>x[i-1]) x[i]=x[i-1]
  }
  x
}

elig[i]

d=dat[294,]

lambda= (1/20) * exp(sum(d[,grep('x',names(d))]*beta))

sdat <- 
  data.frame(
    time=seq(0,10,length.out=100)
  ) %>% 
  mutate(
    surv_true=exp(-lambda * time)
  )

rsf_tree_preds <- 
  map(1:ntree,
    .f=function(tree){
      predictSurvProb(
        object=rsf, 
        newdata=d,
        get.tree=tree,
        times=sdat$time)
    }) %>% 
  reduce(rbind) %>% 
  t() %>% 
  set_colnames(paste0('tree',1:ntree)) %>% 
  as.data.frame() %>% 
  cbind(time=sdat$time) %>% 
  tidyr::gather(variable, value, -time) %>% 
  group_by(variable) %>% 
  mutate(
    value=value+rnorm(100,sd=0.001),
    value=pmin(value,1),
    value=pmax(value,0),
    value=induce_mono(value)
  )

p1 <- rsf_tree_preds %>% 
  ggplot(
    aes(x=time,y=value,group=variable)
  )+
  geom_line(
    alpha=0.02,position=position_jitter()
  )+
  geom_smooth(
    aes(group=1,linetype=''), se=FALSE
  )+
  geom_line(
    data=sdat, aes(x=time, y=surv_true, col=''),
    linetype=2, size=1, inherit.aes = FALSE
  )+
  theme_Publication()+
  labs(
    col='True survival probability',
    linetype='Estimated survival probability',
    title='Random survival forest',
    x='Time since baseline',
    y='Survival probability')

leg <- get_legend(p1)

orsf_tree_preds <- 
  map(orsf$forest,
    .f = predict_ost,
    newx = as.matrix(d),
    times = sdat$time
  ) %>%
  reduce(rbind) %>% 
  t() %>% 
  set_colnames(paste0('tree',1:ntree)) %>% 
  as.data.frame() %>% 
  cbind(time=sdat$time) %>% 
  tidyr::gather(variable, value, -time) %>% 
  group_by(variable) %>% 
  mutate(
    value=value+rnorm(100,sd=0.001),
    value=pmin(value,1),
    value=pmax(value,0),
    value=induce_mono(value)
    )

p2=orsf_tree_preds %>% 
  ggplot(
    aes(x=time,y=value,group=variable)
  )+
  geom_line(
    alpha=0.02,position=position_jitter()
  )+
  geom_smooth(
    aes(group=1,linetype=''), se=FALSE
  )+
  geom_line(
    data=sdat, aes(x=time, y=surv_true, col=''),
    linetype=2, size=1, inherit.aes = FALSE
  )+
  theme_Publication()+
  labs(
    col='True survival probability',
    linetype='Estimated survival probability',
    title='Oblique random survival forest',
    x='Time since baseline',
    y='Survival probability')

result <- grid.arrange(
  p1+theme(legend.position = ''),
  p2+theme(legend.position = ''),
  leg, ncol=2, nrow = 2, 
  layout_matrix = rbind(c(1,2), c(3,3)),
  widths = c(2.7, 2.7), heights = c(2.5, 0.2))

ggsave('simcurve_figure.png', dpi=300, plot=result,
  width = 15, height=7.5)
dev.off()
