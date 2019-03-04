
library(tidyverse)
library(Rcpp)
library(survival)
library(protoclust)
library(magrittr)

files <- list.files(
  file.path(
    '..','ORSF benchmark','Code','ORSF Source'
  ),
  full.names = TRUE
)

for(f in seq_along(files)){
  
  if(f==5){
    sourceCpp(files[f])
  } else {
    source(files[f])
  }
  
}

data = survival::pbc %>% 
  dplyr::select(-id) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::mutate(
    status=case_when(
      status==2 ~ 1L,
      TRUE ~ status
    )
  ) 

clst_indx <- data %>% 
  model.matrix(~.-1, data=.) %>% 
  stats::dist() %>% 
  protoclust::protoclust() %>% 
  protoclust::protocut(k = 3) %>% 
  magrittr::use_series("protos")

obs_to_predict <- data[clst_indx,]

ptimes <- seq(
  quantile(data$time[-clst_indx], probs = 0.10),
  quantile(data$time[-clst_indx], probs = 0.90),
  length.out = 50
)

nboots <- seq(1,250) %>% 
  set_names(paste0("run_",.))
ntree_values <- c(10, 100, 1000) %>% 
  set_names(paste("Trees =",.))

results <- 
  ntree_values %>% 
  map(.f=function(ntrees){
    
    map(nboots, .f=function(iter){
      
      ORSF(
        data=data[-clst_indx,], 
        alpha=0.50, 
        ntree=ntrees, 
        nsplit=25,
        eval_times=ptimes, 
        max_pval_to_split_node=0.50,
        min_obs_to_split_node=20,
        min_obs_in_leaf_node=5,
        min_events_to_split_node=5,
        min_events_in_leaf_node=1,
        verbose=TRUE, 
        use.cv=FALSE
      ) %>% 
        predict.orsf(
          newdata = data[clst_indx, ], 
          times = ptimes
        ) %>% 
        base::t() %>% 
        magrittr::set_colnames(
          paste("Patient class",1:nclust)
        ) %>% 
        base::cbind(time=ptimes) %>% 
        dplyr::as_tibble()
      
    }) %>% 
      dplyr::bind_rows(
        .id = 'run'
      )
    
  }) %>% 
  dplyr::bind_rows(
    .id = 'ntree'
  )

fig <- results %>% 
  tidyr::gather(variable, value, `Patient class 1`:`Patient class 3`) %>% 
  ggplot(
    aes(
      x = time,
      y = value,
      group = run
    )
  ) + 
  geom_line(col='grey70') + 
  geom_smooth(
    col='black', aes(group=1),
    method='gam',formula=y~s(x)) +
  facet_grid(ntree~variable)+
  labs(
    x='Time, days',
    y='Predicted survival probability'
  )

saveRDS(fig, "figure/seed_err_fig.RDS")
