
library(tidyverse)
library(magrittr)
library(haven)
library(magrittr)
library(mlr)

dta = vector(mode='list',length=6)
names(dta)=c('pbc','breast','gbsg2','regards','jhs','vdv')

data("pbc",package='survival')
pbc$status[pbc$status>=1]=pbc$status[pbc$status>=1]-1
pbc$id=NULL
pbc$trt[is.na(pbc$trt)]=0
fctrs<-c('trt','ascites','spiders','edema','hepato','stage')
for(f in fctrs)pbc[[f]]=as.factor(pbc[[f]])

dta[[1]]=pbc

data("Breast",package='biospear')
data=Breast
names(data)[4:ncol(data)]=paste('var',4:ncol(data),sep='_')

dta[[2]]=data

data('GBSG2',package='pec')
data=GBSG2 %>% dplyr::rename(status=cens) %>%
  mutate(tgrade=factor(tgrade,ordered = FALSE))

dta[[3]]=data

data=read.csv('../Datasets/REGARDS_DTH_HF_ASCVD_FU.csv') %>%
  dplyr::rename(dth=death14, time_dth=years_fu_death,
                hf=hf_event, time_hf=years_fu_hf,
                ascvd=CVD, time_ascvd=years_fu) %>% 
  na.omit()

head(data)

cvars = read.csv('../Datasets/REGARDS_calcvars.csv',na.strings = '') %>%
  dplyr::select(id_num=Id_num, Age, Race, Gender,
                Income, ED_Cat, Exercise_cat, Alc_Use, Alc_Drinks_Wk,
                Relationshipstatus, PCS, MCS, Smoke, Packyears,
                starts_with('CAD_'),Diab_SRMed_glu, Dialysis_SR,
                Hyper_SRmeds_BP, Hyper_Meds_SR_now, Insurance,
                Reg_Asa,Gen_SR_Health,KidneyFailure_SR,
                starts_with('Symp_'),CESD,PSS,
                MI_SR_ECG, Afib_SR_ECG,Heartrate,Hct,Hgb,
                Hdl,Ldl,Mch,Mchc,Mcv,Pltc,Rbc,Rdwcv,Glucose,Trigly,
                Wbc,Bun,Crp,Cysc,Height,Weight,Fasted,SBP,DBP,BMI,
                Waist_cm,Albumin_urine,Creatinine_urine,
                #Albumin_serum,
                #Cancer,
                HRT_Meds_SR,
                BirthControl_Meds_SR_ever,
                EGFR_CKDEPI,lvh_main,REGION,Urbangrp)

cvars$Heartrate[cvars$Heartrate==0]=NA

data %<>% 
  left_join(cvars) %>% 
  dplyr::select(-id_num)

for(i in names(data)){
  if(length(unique(data[[i]]))==2) data[[i]] %<>% 
    as.factor() %>% as.numeric() %>% subtract(1)
  if(is.factor(data[[i]])) data[[i]] %<>% 
    factor(labels=gsub(" ", "",levels(data[[i]]),fixed=TRUE))
  if(is.character(data[[i]])) data[[i]] %<>% factor()
  data[[i]][data[[i]]=='']=NA
}

dta[[4]]=droplevels(data)


data = read.csv('../Datasets/JHS_baseline_visit.csv',
                     na.strings = c('.',''))

echo<-haven::read_sas("../Datasets/echa.sas7bdat")%>%
  set_names(tolower(names(.)))%>%mutate(
    lvm_2d=0.8*(1.04*(((echa50/10)+(echa52/10)+
                         (echa54/10))^3-(echa52/10)^3))+0.6,
    lvm_mm=0.8*(1.04*(((echa41/10)+(echa43/10)+
                         (echa45/10))^3-(echa43/10)^3))+0.6)%>%
  dplyr::select(subjid,lvm_2d)%>%
  na.omit()

data %<>% left_join(echo, by = 'subjid') %>%
  mutate(bsa=sqrt((height*weight)/3600),
         lvm=lvm_2d/bsa)

fill = with(data, is.na(AlbuminU24hr) & !is.na(AlbuminUSpot) )
data$AlbuminU24hr[fill]=data$AlbuminUSpot[fill]

fill = with(data, is.na(CreatinineU24hr) & !is.na(CreatinineUSpot) )
data$CreatinineU24hr[fill]=data$CreatinineUSpot[fill]

data %<>% mutate(acr = AlbuminU24hr / CreatinineU24hr) %>% 
  dplyr::select(subjid,age,sex,alc,currentSmoker,
                weight,height,waist,neck,BPmeds,hrtMeds,
                sbp,dbp,abi,HbA1c,FPG,FastHours,
                Diabetes,ldl,hdl,trigs,totchol,LEPTIN,HSCRP,QTcFram,
                ENDOTHELIN,ALDOSTERONE,sCort,SCrCC,
                eGFRckdepi,CHDHx,strokeHx,MedicaidIns,MedicareIns,
                perceivedStress, starts_with('nb'), activeIndex,
                eggs, fish, darkgrnVeg,ecgHR,CV,QRS,
                hyIndex,sportIndex, FEV6, lvm) 

for(i in names(data)){
  if(length(unique(data[[i]]))==2) data[[i]] %<>% as.numeric() 
  if(is.factor(data[[i]])) data[[i]] %<>% 
    factor(labels=gsub(" ", "",levels(data[[i]]),fixed=TRUE))
}


dth="../Datasets/JHS_mortality_2016.csv" %>%
  read.csv() %>% dplyr::select(-X,time_dth=time,dth)

hf="../Datasets/JHS_events_2014_V2.csv" %>% 
  read.csv() %>%
  dplyr::filter(no_fu_consent==0)%>%
  dplyr::select(subjid,hf01,time_hf)

ascvd="../Datasets/JHS_events_2014_V2.csv" %>% 
  read.csv() %>%
  dplyr::filter(no_fu_consent==0)%>%
  dplyr::select(subjid,chd_strk01,time_chd_strk)

data%<>%left_join(dth)%>%left_join(hf)%>%left_join(ascvd)

dta[[5]]=data

data(vdv,package='randomForestSRC')
data=vdv%>%dplyr::rename(status=Censoring,time=Time)
dta[[6]]=data

saveRDS(dta, 'benchmark_data.RDS')

