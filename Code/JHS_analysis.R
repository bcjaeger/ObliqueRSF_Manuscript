
library(tidyverse)
library(magrittr)
library(survival)
library(VIM)
library(party)
library(mlr)
library(obliqueRSF)
library(Rcpp)
library(pec)
library(ggpubr)
library(broman)

custom_cols <- c(
  brocolors('crayons')['Neon Carrot'],
  brocolors('crayons')['Cerulean'],
  brocolors('crayons')['Wisteria']
) %>% 
  as.character()

jhs <- readRDS(
  'Datasets/benchmark_data.RDS'
) %>% 
  .[["jhs"]]

jhs_10yrisk <- read.csv(
  '../Datasets/JHS_10_year_PCR.csv'
) %>% 
  dplyr::select(subjid,risk) 

anly <- jhs %>% 
  dplyr::left_join(jhs_10yrisk, by='subjid') %>% 
  dplyr::select(-subjid) %>% 
  mlr::removeConstantFeatures(perc=0.01) %>% 
  dplyr::filter(!is.na(chd_strk01))

dim(anly)

tmp = anly %>% 
  dplyr::filter(
    BPmeds=='No', strokeHx==1, CHDHx==1, !is.na(risk)
  ) %>% 
  dplyr::mutate(
    cvd_in10 = ifelse(chd_strk01==1 & time_chd_strk < 10, 1, 0)
  ) %>% 
  dplyr::select(risk, cvd_in10) %>% 
  na.omit()

dim(tmp)
summary(tmp)

t1=table(tmp$risk>0.10, tmp$cvd_in10)
t1
t2=prop.table(t1,margin=2)
t2

set.seed(329)

ntrn=3000
trn_indx <- sample(1:nrow(anly),ntrn)
tst_indx <- setdiff(1:nrow(anly), trn_indx)

data <- list(
  train = VIM::kNN(anly[trn_indx,] %>% dplyr::select(-risk)),
  test  = VIM::kNN(anly[tst_indx,] %>% dplyr::select(-risk))
) %>% 
  map(.f=function(df){
    df %>% 
      dplyr::select(
        -contains("hf"),
        -contains("dth")
      ) %>% 
      dplyr::rename(
        time=time_chd_strk,
        status=chd_strk01
      ) 
  })

data$train$lvm[data$train$lvm>300]%<>%divide_by(10)
data$test$lvm[data$test$lvm>300]%<>%divide_by(10)

cif <- pec::pecCforest(
  formula = Surv(time, status)~.,
  data = data$train,
  controls=cforest_unbiased(ntree=1000)
)

orsf <- ORSF(
  data=data$train, 
  ntree=1000
)

prds <- list(orsf=orsf, cif=cif) %>% 
  map(
    ~predictSurvProb(
      object = ., 
      newdata = data$test, 
      times = 10)
  ) %>% 
  map(as.numeric)

prds_survival_curve <- list(orsf=orsf, cif=cif) %>% 
  map(
    ~predictSurvProb(
      object = ., 
      newdata = data$test, 
      times = seq(0,10,length.out=500))
  ) 

survival_indx = c(4, 28, 18, 14)

merge_in <- data$test[survival_indx,] %>% 
  dplyr::select(
    status,
    time_event=time, 
    age, 
    lvm, 
    sbp, 
    eGFRckdepi
  ) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  mutate(
    id=factor(
      paste0(
        "Participant profile no. ", 1:4, ": \n",
        "ASCVD event ", 
        ifelse(
          status==1, 
          'occurred after ', 
          "did not occur in "
        ), 
        time_event, ' years.'
      )
    ),
    lab=paste0(
      "Characteristics: \n",
      "eGFR: ", eGFRckdepi, ' ml/min/1.73 m2 \n',
      "SBP: ", sbp, ' mm Hg \n',
      "Age: ", age, ' years \n',
      "LVM: ", lvm, ' g/m2 \n'
    ),
    time=2,
    value=0.725
  ) 


survival_curves <- prds_survival_curve %>% 
  map(~data.frame(t(.x[survival_indx,]))) %>% 
  bind_cols() %>% 
  as_tibble() %>% 
  set_names(c(paste0("orsf_",1:4), paste0('cif_',1:4))) %>% 
  mutate(time=seq(0,10,length.out=500)) %>% 
  gather(variable, value, -time) %>%
  separate(variable, into=c('variable','id')) %>%
  mutate(
    id=factor(
      id,
      levels = paste(1:4),
      labels = levels(merge_in$id)
    ),
    variable = fct_recode(
      variable, 
      'Oblique random survival forest'='orsf',
      'Conditional inference forest'='cif'
    )
  ) %>% 
  ggplot(aes(x=time,y=1-value))+
  facet_wrap(~id)+
  geom_line(aes(col=variable), size=1.2)+
  geom_text(data=merge_in, aes(label=lab))+
  theme_Publication()+
  scale_color_manual(
    values=custom_cols[1:2]
  ) + 
  labs(
    col='',
    y='Predicted risk of ASCVD event',
    x='Time, years'
    )




prds$risk[is.na(prds$risk)]<-mean(prds$risk,na.rm=T)

pec::cindex(
  object=prds %>% map(as.matrix), 
  data=data$test, 
  eval.times=10,
  formula = Surv(time,status)~1)

calErr <- pec::calPlot(
  object=prds %>% map(as.matrix), 
  data=data$test,
  plot=FALSE,
  eval.times=10,
  method='quantile',
  formula = Surv(time,status)~1)

calErr$plotFrames %>%
  dplyr::bind_rows(.id='model') %>% 
  ggplot(aes(x=Pred,y=Obs,col=model))+
  geom_line()+
  geom_point(size=2)+
  theme_Publication()+
  geom_abline(slope=1,intercept=0,col='black',linetype=2)

calErr$plotFrames %>% 
  dplyr::bind_rows(.id='model') %>% 
  dplyr::group_by(model) %>% 
  mutate(
    quantile=factor(1:n()),
    calErr = abs(Pred-Obs)
  ) %>% group_by(model) %>% 
  summarise(err=mean(calErr))

cols_to_bind <- prds %>% map(as.numeric)
names(cols_to_bind) <- c("orsf","cif","risk")

ggdat <- data$test %>% 
  dplyr::bind_cols(cols_to_bind)



xvars <- list(
  lvm=list(limits=c(40,100),label='Left ventricular mass, g/m2'),
  age=list(limits=c(20,80),label='Age, years'),
  sbp=list(limits=c(100,180),label='Systolic blood pressure, mm Hg'),
  eGFRckdepi = list(
    limits=c(20,100),
    label='eGFR, ml/min/1.73 m2'
  )
)

desc_plot <- data %>% 
  bind_rows(.id='set') %>% 
  .[,c("set",names(xvars))] %>% 
  gather(variable,value,-set) %>% 
  mutate(
    variable=factor(
      variable, 
      levels=names(xvars), 
      labels=map_chr(xvars,~.$label)
    ),
    set=fct_recode(set,
      "Testing"="test",
      "Training"="train")
  ) %>% 
  ggplot(aes(x=value,fill=set))+
  geom_histogram(alpha=0.80,col='black',bins=30)+
  facet_wrap(~variable, scales='free')+
  theme_Publication()+
  scale_fill_manual(values=custom_cols[1:2])+
  labs(
    fill='Dataset',
    x='Value',
    y='Count'
  )


plts <- map2(xvars,names(xvars),.f=function(.xvar,.lab){
  
  ggdat %>% 
    dplyr::select(!!.lab, orsf, cif, risk) %>% 
    .[.[[.lab]]>=.xvar$limits[1] & .[[.lab]]<=.xvar$limits[2],] %>% 
    tidyr::gather(variable, value, -!!.lab) %>% 
    dplyr::mutate(
      value = 1-value,
      variable = forcats::fct_recode(
        variable,
        'Conditional inference forest'='cif',
        'Oblique random survival forest'='orsf',
        'Pooled Cohort risk equation'='risk'
      ) 
    ) %>% 
    ggplot(aes_string(x=.lab,y='value',col='variable'))+
    labs(x=.xvar$label)+
    geom_smooth(method='gam',formula=y~s(x,bs='cs'),se=FALSE)+
    theme_Publication()+
    scale_color_manual(values=custom_cols)+
    labs(col='',y="10-year predicted risk of ASCVD event")
  
})


vdep_plot <- 
  ggarrange(
    plotlist=plts, 
    common.legend = TRUE,
    legend='bottom'
  )

# Partial variable dependence ---------------------------------------------

xgrid <- xvars %>% 
  map(~seq(
    .$limits[1],.$limits[2],length.out=10
  )
  ) %>% bind_cols() %>% 
  tidyr::gather() %>% 
  group_by(key) %>% 
  nest() %>%
  mutate(pdep = map2(data,key, gen_partial_preds))

pd=position_dodge(width=2)

xgrid2 = xgrid %>% 
  mutate(
    plts = map2(
      key, pdep, .f=function(lab,df){
        df %>% 
          mutate(
            model = forcats::fct_recode(
              model,
              'Conditional inference forest'='cif',
              'Oblique random survival forest'='orsf'
            )
          ) %>% 
          ggplot(
            aes(
              x=xval,
              y=1-mean,
              ymin=1-mean-2*stde,
              ymax=1-mean+2*stde,
              col=model
            )
          ) +
          geom_point(size=2, position=pd)+
          geom_errorbar(position=pd)+
          labs(
            col='',
            x=xvars[[lab]]$label,
            y='10-year predicted risk of ASCVD event')+
          theme_Publication()+
          scale_color_manual(values=custom_cols[1:2])
      }
    )
  )

pdep_plot <- 
  ggarrange(
    plotlist=xgrid2$plts, 
    common.legend = TRUE,
    legend='bottom'
  )


gen_partial_preds <- function(grid, lab){
  
  map(
    grid[[1]], 
    .f=function(xval){
      
      tmp=data$test
      tmp[[lab]]=xval
      
      list(orsf=orsf, cif=cif) %>% 
        map(
          ~predictSurvProb(
            object = ., 
            newdata = tmp, 
            times = 10
          )
        ) %>% 
        map(as.numeric)%>% 
        bind_cols() %>% 
        gather() %>% 
        group_by(key) %>% 
        summarise(
          mean=mean(value),
          stde=sd(value)/sqrt(n()),
          xval=xval
        ) %>% 
        rename(model=key)
    }
  ) %>% 
    bind_rows()
  
}


# Means and standard devs
library(gt)
desc_tab = data %>% 
  bind_rows(.id='set') %>% 
  .[, c("set",names(xvars))] %>% 
  as_tibble() %>% 
  gather(variable, value, -set) %>% 
  group_by(set, variable) %>% 
  summarise(
    Minimum=min(value),
    `25% Quartile`=quantile(value, probs = 1/4),
    Median=quantile(value, probs = 2/4),
    `75% Quartile`=quantile(value, probs = 3/4),
    Maximum=max(value),
    Mean=mean(value),
    `Standard Deviation` =sd(value)
  ) %>% ungroup() %>%  
  mutate(
    variable=factor(
      variable,
      levels=names(xvars),
      labels=map_chr(xvars, ~.$label)
    ),
    set=factor(
      set,
      levels=c('test','train'),
      labels=c("Testing data", "Training data")
    )
  ) %>% 
  group_by(set) %>% 
  gt() %>% 
  fmt_number(
    columns=vars(
      Minimum, 
      `25% Quartile`, 
      Median, 
      `75% Quartile`, 
      Maximum,
      Mean, 
      `Standard Deviation`), 
    decimals = 2
  )

saveRDS(survival_curves,"figure/survival_curves_plot.RDS")
saveRDS(pdep_plot,'figure/pdep_plot.RDS')
saveRDS(desc_plot,'figure/desc_plot.RDS')
saveRDS(vdep_plot, "figure/vdep_plot.RDS")
saveRDS(desc_tab,'figure/desc_tab.RDS')  

#