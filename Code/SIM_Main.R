

# Setup -------------------------------------------------------------------


require(survival)
require(randomForestSRC)
require(pec)
require(rms)
require(magrittr)
require(glmnet)
require(purrr)
require(dplyr)
require(penalized)
require(ggRandomForests)
require(haven)
require(MASS)
require(c060)
require(peperr)
require(party)
require(CoxBoost)
require(mlr)
require(xgboost)
require(rBayesianOptimization)
require(prodlim)
require(simsurv)
require(parallel)
require(obliqueRSF)

# sample sizes
n=1500; ntrn=1000; ntst=n-ntrn
# random forest hyper-parameters
minsplit=10; ntree=1000; nsplit=25; verbose=T; mpval=0.50; alpha=c(0.95)
# Numbers of predictor variables / interactions / junk variables
ncov=25; nint=ncov; njnk=ncov
# Correlation (autoregressive constant)
rho=0.33
# Administrative censoring time
maxt=10
# Scenario: 1, 2, or 3.
scenario=2

# Effect sizes ------------------------------------------------------------

# Main effects 
xnames=paste0('x',1:ncov)
beta <- seq(-1,1,length.out = ncov) %>% 
  divide_by(sqrt(ncov))%>%
  set_names(xnames)

# Interaction effects
int_ind <- vector(mode = 'list',length = nint)
for(i in 1:nint) int_ind[[i]]= if(i<ncov) c(i,i+1) else c(i,i+1-ncov)

icoefs <- seq(-1,1,length.out = nint) %>% 
  divide_by(sqrt(ncov)) %>% 
  set_names(paste0('i',1:nint))

# scenario 1 does not include interaction effects
# so only add the interactions if scenario = 2 or 3
if(scenario!=1) beta %<>% c(icoefs)

# Covariance matrix between fixed effects
Sigma=matrix(0,ncol=ncov,nrow=ncov)

# Impose an autoregressive covariance structure
# between the fixed effects. Since each x variable
# has unit variance, cov(X) = corr(X)
for(i in 1:ncov){
  for(j in 1:ncov){
    Sigma[i,j] = rho^(abs(i-j))
  }
}


# Random stream initiation ------------------------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(329)
niter <- 250 
results <- vector(mode='list',length=niter)


# Main loop ---------------------------------------------------------------

for(i in 1:niter){
  
  .Random.seed <- nextRNGStream(.Random.seed)
  # Data Generation ---------------------------------------------------------
  
  # Generate multivariate normal design matrix with
  # specified covariance and zero means.
  X=x=MASS::mvrnorm(n=n,mu=rep(0,ncov),Sigma=Sigma) %>%
    data.frame() %>% 
    set_names(xnames)
  
  if(scenario!=1){
    
    intr = map(int_ind,~x[,.[1]]*x[,.[2]]) %>% 
      reduce(cbind) %>% data.frame() %>% 
      set_names(names(icoefs))
    
    X %<>% cbind(intr)
    
  }
  
  # simulation parameters
  gammas=1
  
  s1<-simsurv(dist='exponential',lambdas=1/20,betas=beta,x=X,maxt=maxt)
  
  #hist(s1$eventtime)
  
  if(scenario==3){
    
    s2<-simsurv(lambdas=1/100,gammas=gammas,betas=beta,x=X,maxt=maxt)
    s3<-simsurv(lambdas=1/2, gammas=gammas, betas=beta,x=X,maxt=maxt)
    s4<-simsurv(lambdas=1/10, gammas=gammas, betas=beta,x=X,maxt=maxt)
    
    # hist(s2$eventtime)
    # hist(s3$eventtime)
    # hist(s4$eventtime)
    
    lc1=apply(x[,1:(ncov/2)],1,sum)
    lc2=apply(x[,(ncov/2):ncov],1,sum)
    
    grp1=which(lc1 <0 & lc2 <0) # lc1 < 0 => healthy, lc2 < 0 => unhealthy
    grp2=which(lc1 <0 & lc2>=0) # healthy healthy
    grp3=which(lc1>=0 & lc2 <0) # unhealthy unhealthy
    grp4=which(lc1>=0 & lc2>=0) # unhealthy healthy
    
  }
  
  dat=cbind(s1,x) %>% data.frame() %>%
    dplyr::select(-starts_with('i')) %>% 
    dplyr::rename(time=eventtime)
  
  if(njnk>0){
    dat[,paste0('z',1:njnk)]=rnorm(nrow(dat)*njnk,mean=0,sd=1)
  }
  
  if(scenario==3){
    
    dat[grp2,c('time','status')]=s2[grp2,c('eventtime','status')]
    dat[grp3,c('time','status')]=s3[grp3,c('eventtime','status')]
    dat[grp4,c('time','status')]=s4[grp4,c('eventtime','status')]
    
  }
  
  ftrs=setdiff(names(dat),c('time','status'))
  mat=model.matrix(~.,data=dat[,ftrs])
  
  # Data splitting ----------------------------------------------------------
  
  trn_indx <- sample(1:nrow(dat), round(nrow(dat)*1/2))
  tst_indx <- setdiff(1:nrow(dat), trn_indx)
  trn = dat[trn_indx,]
  tst = dat[tst_indx,]
  
  # Evaluation times --------------------------------------------------------
  
  mintime=max(
    c(min(tst$time[tst$status==1]),
      min(trn$time[trn$status==1]))
  )
  
  maxtime=min(
    c(median(tst$time[tst$status==1]),
      median(trn$time[trn$status==1]))
  )
  
  ptimes = seq(mintime,maxtime,length.out=100)
  
  # Random Survival Forests -------------------------------------------------
  
  start=Sys.time()
  rsf = rfsrc(Surv(time,status)~.,
    data=trn,
    ntree=ntree,
    nsplit=nsplit,
    nodesize=minsplit)
  stop=Sys.time()
  rsf_time=stop-start
  
  # Oblique random survival forests -----------------------------------------
  
  start=Sys.time()
  orsf.cv=ORSF(data=trn, alpha=alpha,
    ntree=ntree, nsplit=nsplit,
    eval_times=ptimes,
    max_pval_to_split_node=mpval,
    min_obs_to_split_node=minsplit,
    min_obs_in_leaf_node = 5,
    min_events_to_split_node = 5,
    gamma = 1/3,
    min_events_in_leaf_node = 1,
    verbose=verbose, use.cv=TRUE)
  stop=Sys.time()
  orsf_time_cv=stop-start
  
  start=Sys.time()
  orsf=ORSF(data=trn, alpha=0.50, 
    ntree=ntree, nsplit=nsplit,
    eval_times=ptimes, 
    max_pval_to_split_node=mpval,
    min_obs_to_split_node=minsplit,
    min_obs_in_leaf_node = 5,
    min_events_to_split_node = 5,
    gamma = 1/3,
    min_events_in_leaf_node = 1,
    verbose=verbose, use.cv=FALSE)
  stop=Sys.time()
  orsf_time=stop-start
  
  
  # Forward Selection -------------------------------------------------------
  
  if(ncol(trn)>1000){
    step=NULL
    step_time=0
  } else {
    scope=list(lower=Surv(time,status)~1,
      upper=paste('Surv(time,status)~',
        paste(ftrs,collapse='+')) %>%
        as.formula())
    start=Sys.time()
    step=stepAIC(coxph(Surv(time,status)~1,data=trn,x=TRUE), 
      scope = scope,direction='forward', 
      trace=verbose, steps=30)
    stop=Sys.time()
    step_time=stop-start
  }
  
  # Regularized Cox PH models -----------------------------------------------
  
  trnx = model.matrix(~., data=trn[,ftrs])[,-1L]
  tstx = model.matrix(~., data=tst[,ftrs])[,-1L]
  
  start=Sys.time()
  cox_lasso = cv.glmnet(trnx, Surv(trn$time,trn$status),
    family="cox", alpha=0.90)
  stop=Sys.time()
  lasso_time=stop-start
  
  start=Sys.time()
  cox_ridge = cv.glmnet(trnx, Surv(trn$time,trn$status),
    family="cox", alpha=0.10)
  stop=Sys.time()
  ridge_time=stop-start
  
  start=Sys.time()
  cox_bst = pec::coxboost(Hist(time,status)~.,data=trn)
  stop=Sys.time()
  bst_time=stop-start
  
  
  # Conditional Inference Forest --------------------------------------------
  
  start=Sys.time()
  cif_cntrl = cforest_unbiased(minsplit=minsplit,ntree=ntree)
  cif <- pec::pecCforest(Surv(time,status)~., data=trn, controls=cif_cntrl)
  stop=Sys.time()
  cif_time=stop-start
  
  
  # Xgboost -----------------------------------------------------------------
  
  xgb.label=trn$time
  xgb.label[trn$status==0]%<>%multiply_by(-1)
  
  mat=list(trn=trn,tst=tst)%>%
    map(.f=function(d){
      d$time=NULL;d$status=NULL
      d=model.frame(~., data=d, na.action='na.pass')
      model.matrix(~.,data=d)[,-1L]
    })
  
  dmat=list(
    train=xgb.DMatrix(data=mat$trn,label=xgb.label),
    test=xgb.DMatrix(data=mat$tst)
  )
  
  start=Sys.time()
  
  xgb<-xgb_surv_cv(data=dmat$train,nfold=3,
    eval.times=ptimes,init_points=20,n_iter=0,
    eta=0.05)
  
  stop=Sys.time()
  xgb_time=stop-start
  
  # Prediction --------------------------------------------------------------
  
  mdls = list('orsf.cv'=orsf.cv,
    'orsf'=orsf,
    'bst'=cox_bst,
    'cif'=cif, 
    'rsf'=rsf,
    'step'=step)
  
  if(any(map_lgl(mdls,is.null)))mdls=mdls[-which(map_lgl(mdls,is.null))]
  
  start=Sys.time()
  prbs = lapply(mdls, predictSurvProb, newdata=tst, times=ptimes)
  stop=Sys.time()
  prbs_time=stop-start
  
  lasso_indx=which(cox_lasso$lambda==cox_lasso$lambda.min)
  nunique=1
  ntry=0
  
  while(lasso_indx>=1 & nunique==1){
    prbs$lasso=c060:::fit.glmnet(x=trnx,response=Surv(trn$time,trn$status),
      family="cox",alpha=0.90,
      cplx=cox_lasso$lambda[lasso_indx]) %>% 
      predictProb(response=Surv(trn$time,trn$status), 
        x=tstx, times=ptimes, 
        complexity=cox_lasso$lambda.min)
    
    nunique=length(unique(prbs$lasso[,length(ptimes)]))
    ntry=ntry+1
    lasso_indx=lasso_indx-1
  }
  
  if(lasso_indx==1) prbs$lasso=NULL
  
  prbs$ridge=c060:::fit.glmnet(x=trnx,response=Surv(trn$time,trn$status),
    family="cox",alpha=0.10,
    cplx=cox_ridge$lambda.min) %>% 
    predictProb(response=Surv(trn$time,trn$status), 
      x=tstx, times=ptimes, 
      complexity=cox_ridge$lambda.min)
  
  bh=gbm::basehaz.gbm(t=trn$time,
    delta=trn$status,
    f.x=xgb$preds,
    t.eval=ptimes,
    smooth=FALSE,
    cumulative=TRUE)
  
  xgb.prds <- predict(xgb$booster,newdata=dmat$test,outputmargin=T)
  prbs$xgbst <- matrix(0, nrow=nrow(tst),ncol=length(ptimes))
  
  for(i in 1:length(bh)){
    prbs$xgbst[,i]=exp(-exp(xgb.prds) * (bh[i]))
  }
  
  # Evaluation --------------------------------------------------------------
  
  grab_last=function(x) x[length(x)]
  
  cstat = pec::cindex(
    prbs, data=tst, formula=Surv(time,status)~1,eval.times=ptimes) %>%
    use_series('AppCindex') %>% map2(names(.),function(x,y){
      data.frame(mdl=y,time=ptimes,err=1-x)
    }) %>%
    reduce(rbind) %>%
    mutate(seed=seed,measure='cstat')
  
  .cstat = cstat %>% group_by(mdl) %>% 
    dplyr::filter(time==max(time)) %>% 
    dplyr::select(-time) %>% data.frame()
  
  rmsqd<-function(data,columns=c(1,2)){
    sqr_diff=(data[,columns[1]]-data[,columns[2]])^2
    sqrt(mean(sqr_diff))
  }
  
  .calib=try(map(prbs,.f=function(pb){
    pec::calPlot(pb[,length(ptimes)],plot=T,
      time=ptimes[length(ptimes)],
      method='quantile',q=5,
      formula=Surv(time,status)~1,
      data=tst)[[1]]%>%
      magrittr::use_series('Model.1')%>%
      rmsqd(columns=c(1,2))
  }) %>% reduce(rbind) %>%
      data.frame() %>%
      set_names("err")%>%
      mutate(mdl=names(prbs),seed=seed,
        measure='calibration'))
  
  perr = pec::pec(prbs,data=tst,times=ptimes,exact=F,
    start=ptimes[1],maxtime=ptimes[length(ptimes)],
    formula=Surv(time,status)~1,
    cens.model="cox") %>% 
    use_series("AppErr") %>% map2(names(.),function(x,y){
      data.frame(mdl=y,time=ptimes,err=x)
    }) %>% reduce(rbind) %>% 
    mutate(seed=seed, measure='brier_score')
  
  .perr=crps(pec::pec(prbs,data=tst,times=ptimes,exact=F,
    start=ptimes[1],maxtime=ptimes[length(ptimes)],
    formula=Surv(time,status)~1,
    cens.model="cox"))
  .perr=data.frame(mdl=rownames(.perr),
    err=as.numeric(.perr),
    seed=seed,
    measure='brier_score')
  
  
  pdata=map(prbs,~.x[,ncol(x)])%>%bind_cols()
  
  # Store results -----------------------------------------------------------
  
  if(class(.calib)[1]=='try-error'){
    smry=rbind(.cstat,.perr)
  } else {
    smry=rbind(.cstat,.perr,.calib)
  }
  
  results[[i]] <- 
    list(
      full=rbind(cstat,perr),
      smry=smry,
      time=c(rsf_time,orsf_time_cv,orsf_time,cif_time,step_time,
        lasso_time,ridge_time,bst_time,xgb_time)
    )
  
  
}


