# Multinational evaluation of AnthropoAge as a measure of biological age
# in the USA, England, Mexico, Costa Rica, and China:
# a population-based longitudinal study
## Data Analysis: Carlos Alberto Fermin-Martinez & Omar Yaxmehen Bello-Chavolla
## Latest version of Analysis: May 2024
## For any question regarding analysis contact:
## Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#PENDIENTES:
#5. OUTCOMES ADICIONALES EN LAS VARIABLES EN LAS QUE SEA POSIBLE... ¿TELOMERE LENGHT? ¿LABS? ¿PHENOAGE?

####-------------------------------- DATABASES --------------------#### ----####
#Preparation----- ####
### Working directory ###
wd.CAFM <- paste0( #CAFM working directory
  "~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/",
  "AnthropoAge - Secundarios/AnthropoAge - G2A")
wd.OYBC <- paste0( #OYBC working directory
  "~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/",
  "AnthropoAge - Secundarios/AnthropoAge - G2A")
#FUN: Evaluates multiple working directories and sets the first existing one
set_existing_wd <- function(wds) {
  for (wd in wds) {
    if (file.exists(wd) && file.info(wd)$isdir) {
      setwd(wd); cat("Working directory set to:", wd, "\n")
      return(invisible())}}
  stop("None of the provided directories exist or are accessible.\n")}
set_existing_wd(c(wd.OYBC,wd.CAFM))

### Packages ###
#devtools::install_github("davidsjoberg/ggsankey")
#remotes::install_github("chainsawriot/hongkong")
pacman::p_load(tidyverse, survival, haven, flexsurv, survey,   flextable,
               ggpubr,blandr, BlandAltmanLeh, jtools,ggsankey, patchwork,
               sjmisc, pROC, gsubfn, lubridate, gridExtra, ggpubr, furrr,
               OptimalCutpoints, timeROC, lme4, lmerTest, coxme, geepack,
               weights, officer, riskRegression, compareC, parallel, gee,
               ggalluvial, devtools, gtools,foreach,doParallel,coxme,rms,
               iterators, measurements, survminer, gtsummary, data.table)

### Race/ethnicity coding ###
race_ethnicity1 <- c(
  "Non-Hispanic White", "Non-Hispanic Black", "Hispanic/Latino", "Other")
race_ethnicity2 <- c(
  "White", "Black", "Mexican-American", "Other")

### Functions ###
source("Models/functions_predict.R")
s_anthropoage<-function(x){
  # Data frame including Age (years), Sex (coded as "Men" and "Women"),
  # Height (cm), Weight (kg), Waist (cm) and Ethnicity (coded as
  # "White", "Black", "Mexican-American" and "Other")
  x<-as.data.frame(x)
  x$tr_imc <- log(x$BMI) #BMI is log-transformed
  x$tr_ice <- x$ICE**(1/3) #WHtR is cubic root-transformed
  x$Ethnicity <- case_when( #Correct race/ethnicity coding
    x$Ethnicity==race_ethnicity1[1]~race_ethnicity2[1],
    x$Ethnicity==race_ethnicity1[2]~race_ethnicity2[2],
    x$Ethnicity==race_ethnicity1[3]~race_ethnicity2[3],
    x$Ethnicity==race_ethnicity1[4]~race_ethnicity2[4])
  ## Load models ##
  model7<-load("Models/7_SAnthropo_M.rda") #Anthropometry model (men)
  model8<-load("Models/8_SAnthropo_F.rda") #Anthropometry model (women)
  model9<-load("Models/9_SAge_M.rda") #Age model (men)
  model10<-load("Models/10_SAge_F.rda") #Age model (women)
  ##Estimate model coefficients based on NHANES models##
  sM1<-1/((exp(coef(gomp1bM1)[1]*120)-1)/((coef(gomp1bM1)[1])))
  b0M1<-coef(gomp1bM1)[2]
  b1M1<-coef(gomp1bM1)[3]
  sW1<-1/((exp(coef(gomp1bF1)[1]*120)-1)/((coef(gomp1bF1)[1])))
  b0W1<-coef(gomp1bF1)[2]
  b1W1<-coef(gomp1bF1)[3]
  ## Loop for S-AnthropoAge estimation ##
  for(i in 1:nrow(x)){
    if(x$Sex[i]=="Women"){
      p1F<-predict(gomp1aF1, newdata = x,type="survival", ci=F, times = c(120))
      pred<-as.numeric(1-p1F$.pred)
      output<-(log(-sW1*log(1-pred))-b0W1)/b1W1}
    else if(x$Sex[i]=="Men") {
      p1M<-predict(gomp1aM1, newdata = x,type="survival", ci=F, times = c(120))
      pred<-as.numeric(1-p1M$.pred)
      output<-(log(-sM1*log(1-pred))-b0M1)/b1M1}
    ##Output
    output[is.infinite(output)]<-NA
    return(output)}}
s_anthropoage_fast<-function(x){
  # Data frame including Age (years), Sex (coded as "Men" and "Women"),
  # Height (cm), Weight (kg), Waist (cm) and Ethnicity (coded as
  # "White", "Black", "Mexican-American" and "Others")
  x<-as.data.frame(x)
  x$tr_imc <- log(x$BMI) #BMI is log-transformed
  x$tr_ice <- x$ICE**(1/3) #WHtR is cubic root-transformed
  x$Ethnicity <- case_when( #Correct race/ethnicity coding
    x$Ethnicity==race_ethnicity1[1]~race_ethnicity2[1],
    x$Ethnicity==race_ethnicity1[2]~race_ethnicity2[2],
    x$Ethnicity==race_ethnicity1[3]~race_ethnicity2[3],
    x$Ethnicity==race_ethnicity1[4]~race_ethnicity2[4])
  ## Load models ##
  model7<-load("Models/7_SAnthropo_M.rda") #Anthropometry model (men)
  model8<-load("Models/8_SAnthropo_F.rda") #Anthropometry model (women)
  model9<-load("Models/9_SAge_M.rda") #Age model (men)
  model10<-load("Models/10_SAge_F.rda") #Age model (women)
  ##Estimate model coefficients based on NHANES models##
  sM1<-1/((exp(coef(gomp1bM1)[1]*120)-1)/((coef(gomp1bM1)[1])))
  b0M1<-coef(gomp1bM1)[2]
  b1M1<-coef(gomp1bM1)[3]
  sW1<-1/((exp(coef(gomp1bF1)[1]*120)-1)/((coef(gomp1bF1)[1])))
  b0W1<-coef(gomp1bF1)[2]
  b1W1<-coef(gomp1bF1)[3]
  ## Loop for parallel S-AnthropoAge estimation ##
  cl <- makeCluster(detectCores()); registerDoParallel(cl); chunk_size <- 1000
  it <- iter(seq(1, nrow(x), by = chunk_size), chunksize = chunk_size)
  for(i in 1:nrow(x)){
    if(x$Sex[i]=="Women"){ #AnthropoAge for women
      predict_chunk <- function(idx) {chunk <- x[idx:min(
        idx + chunk_size - 1, nrow(x)), ]
      predict(gomp1aF1, newdata = chunk, type = "survival",
              ci = FALSE, times = c(120))}
      p1F<-foreach(idx = (it[["state"]][["obj"]]),
                   .combine = rbind) %dopar% {predict_chunk(idx)}
      pred<-as.numeric(1-p1F$.pred_survival)
      output<-(log(-sW1*log(1-pred))-b0W1)/b1W1}
    else if(x$Sex[i]=="Men") { #AnthropoAge for men
      predict_chunk <- function(idx) {chunk <- x[idx:min(
        idx + chunk_size - 1, nrow(x)), ]
      predict(gomp1aM1, newdata = chunk, type = "survival",
              ci = FALSE, times = c(120))}
      p1M<-foreach(idx = (it[["state"]][["obj"]]),
                   .combine = rbind) %dopar% {predict_chunk(idx)}
      pred<-as.numeric(1-p1M$.pred_survival)
      output<-(log(-sM1*log(1-pred))-b0M1)/b1M1}
    stopCluster(cl)
    ##Output
    output[is.infinite(output)]<-NA
    return(output)}}
weighted.se.mean <- function(x, w, na.rm = T){
  ## Remove NAs 
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]}
  
  ## Calculate effective N and correction factor
  n_eff <- (sum(w))^2/(sum(w^2))
  correction = n_eff/(n_eff-1)
  
  ## Get weighted variance 
  numerator = sum(w*(x-weighted.mean(x,w))^2)
  denominator = sum(w)
  
  ## get weighted standard error of the mean 
  se_x = sqrt((correction * (numerator/denominator))/n_eff)
  return(se_x)}
fllwp.yrs <- function(x){
  y <- x/12; z <- quantile(y, na.rm=T) %>% round(1)
  a <- paste0(z[3], " (", z[2], "–", z[4], ")"); a}
fllwp.yrs2 <- function(x){y <- x/12; z1 <- mean(y, na.rm=T) %>% round(1)
z2 <- sd(y, na.rm=T) %>% round(1); a <- paste0(z1, " ± ", z2, " years"); a}
confint.geeglm <- function(object, parm, level = 0.95, ...) {
  cc <- coef(summary(object))
  mult <- qnorm((1+level)/2)
  citab <- with(as.data.frame(cc),
                cbind(lwr=Estimate-mult*Std.err,
                      upr=Estimate+mult*Std.err))
  rownames(citab) <- rownames(cc)
  citab[parm,]
}
nrowf2 <- function(x){nrow(x) %>% format(big.mark = ",")}

### Table of AnthropoAge calculation ###
#model7<-load("Models/7_SAnthropo_M.rda") #Anthropometry model (men)
#model8<-load("Models/8_SAnthropo_F.rda") #Anthropometry model (women)
#model9<-load("Models/9_SAge_M.rda") #Age model (men)
#model10<-load("Models/10_SAge_F.rda") #Age model (women)
#options(scipen=100)
#model9 %>% get %>% coefficients %>% round(5) #Men
#model10 %>% get %>% coefficients %>% round(5) #Women
#model7 %>% get %>% coefficients %>% round(5) #Men
#model8 %>% get %>% coefficients %>% round(5) #Women


#HRS------------- ####
## HRS ## -- 15 Waves: 1992-2020 (biennially)
## S-AnthropoAge can only be calculated in W8-W14 (2006,08,10,12,14,16,18)
HRS_alt <- readRDS("Databases/randhrs1992_2020v1.rds.gz")

#Way too heavy:
#HRS_h <- readRDS("Databases/HRS/H_HRS_d.rds.gz")
#HRS_h_var <- HRS_h %>% select(hhidpn, raeducl, all_of(paste0(
#  "r",8:14,"smokef"))) %>% rename("HHIDPN"=hhidpn)
#remove(HRS_h)
#saveRDS(HRS_h_var, gzfile("Databases/HRS/HRS_h_var.rds.gz"))

#Harmonized HRS to obtain demographics and lifestyle data
HRS_h_var <- readRDS("Databases/HRS/HRS_h_var.rds.gz") %>% 
  `names<-`(c("HHIDPN","raeducl", paste0("R",8:14,"SMOKEF")))
HRS_alt <- dplyr::bind_cols(HRS_alt, HRS_h_var[-1])

# Select variables and change to long format
HRS.A1.0 <- HRS_alt %>% select(
  HHIDPN, RAGENDER, RARACEM, RAHISPAN, raeducl, RADMONTH, RADYEAR,
  R15IWSTAT, HACOHORT, all_of(paste0("INW",8:14)),
  (8:14 %>% sapply(gsub, pattern="@", x=c(
    "R@WTRESP", "R@IWENDY", "R@IWENDM","R@AGEY_E",
    "R@PMHGHT", "R@PMWGHT", "R@PMWAIST", "R@SHLT", "R@BATHA",
    "R@EATA", "R@BEDA", "R@WALKRA", "R@TOILTA", "R@MONEYA",
    "R@MEDSA", "R@SHOPA", "R@MEALSA", "R@HIBPE", "R@DIABE",
    "R@CANCRE", "R@LUNGE", "R@HEARTE", "R@STROKE", "R@ARTHRE",
    "R@SMOKEV","R@SMOKEN","R@SMOKEF","R@DRINK","R@DRINKD","R@VGACTX")) %>%
     as.character), R15IWENDY, R15IWENDM); HRS.A1.0_names <- HRS.A1.0 %>%
  select(-c(1:16)) %>% names; HRS.A1.0 <- HRS.A1.0 %>% pivot_longer(cols = -c(
    1:16, (HRS.A1.0_names[!grepl(pattern="^R(\\d+)(.*)", HRS.A1.0_names)])),
    names_to = c("wave",".value"), names_pattern = "^R(\\d+)(.*)")

# Create new variables
HRS.A1.0$Ethnicity <- HRS.A1.0$RARACEM; HRS.A1.0$Ethnicity[
  HRS.A1.0$RAHISPAN==1]<-4; HRS.A1_W <- HRS.A1.0 %>%
  mutate( #Ethnicity, anthropometry, health decline
    "ID"=HHIDPN, "Sex"=factor(RAGENDER,1:2,c("Men","Women")),
    "G2ASTUDY"="HRS", "mortstat"=R15IWSTAT, "Ethnicity"=factor(
      Ethnicity, c(1,2,4,3), race_ethnicity1), "Age"=AGEY_E, "SRH"=SHLT,
    "ADL"=(BATHA + EATA + BEDA + WALKRA + TOILTA), #ADL (0-5)
    "IADL"=(MONEYA + MEDSA + SHOPA + MEALSA), #IADL (0-4)
    "Height"=PMHGHT*100, "Weight"=PMWGHT, "Waist"=PMWAIST*2.54,
    "B.IY"=IWENDY, "B.IM"=IWENDM, "HBP"=HIBPE, "T2D"=DIABE,
    "CAN"=CANCRE, "LUD"=LUNGE, "AMI"=HEARTE, "EVC"=STROKE,
    "ART"=ARTHRE, "BMI"=Weight/((Height/100)^2), "ICE"=Waist/Height) %>%
  mutate( #Date of baseline interview and date of death
    "B.DATE"=paste(B.IY,B.IM,"1",sep="-") %>% as.Date("%Y-%m-%d"),
    "D.DATE"=paste(RADYEAR,RADMONTH,"1",sep="-")%>%as.Date("%Y-%m-%d"),
    "year"=factor(wave,8:14,seq(2006,2018,2)) %>% as.character %>%
      as.numeric); HRS.A1_W<-HRS.A1_W %>%
  mutate( #AnthropoAge and additional variables
    "AnthropoAge" = s_anthropoage_fast(HRS.A1_W),
    "bi"=case_when(
      year%in%2005:2006~"2005/06", year%in%2007:2008~"2007/08",
      year%in%2009:2010~"2009/10", year%in%2011:2012~"2011/12",
      year%in%2013:2014~"2013/14", year%in%2015:2016~"2015/16",
      year%in%2017:2018~"2017/18"),
    "participation"=(INW8==1)+(INW9==1)+(INW10==1)+(INW11==1)+
      (INW12==1)+(INW13==1)+(INW14==1),
    "consecutives"=(INW8==1&INW9==1)|(INW9==1&INW10==1)|(INW10==1&INW11==1)|
      (INW11==1&INW12==1)|(INW12==1&INW13==1)|(INW13==1&INW14==1),
    "coh"=case_when(HACOHORT%in%1:5~"A", HACOHORT%in%6~"B", HACOHORT%in%7~"C"))
# Date of last interview
setDT(HRS.A1_W); HRS.A1_W[, LAST.IW := max( B.DATE, na.rm = TRUE), by = ID]
HRS.A1_W$LAST.IW[is.infinite(HRS.A1_W$LAST.IW)] <- NA
# Mortality status (0=Alive | 1=Dead)
HRS.A1_W$mortstat2<-ifelse(HRS.A1_W$mortstat %in% 5:6, 1, 0)
HRS.A1_W$mortstat2[HRS.A1_W$mortstat==0] <- NA
# Censoring date: either date of death or date of last interview
HRS.A1_W <- HRS.A1_W %>% mutate("E.DATE" = ifelse(
  mortstat2==0|is.na(D.DATE), as.character(LAST.IW),
  as.character(D.DATE)) %>% as.Date)

# Filters
HRS.A1 <- HRS.A1_W %>%
  mutate( #Time to death in months + country specific race/ethnicity
    "TTDM"=interval(start=B.DATE,end=E.DATE) %/% months(1),
    "Ethnicity2"=factor(Ethnicity, levels=race_ethnicity1,
                        labels=paste0("US-",race_ethnicity1))) %>%
  filter( #Complete follow-up data (≥0 months)
    TTDM>=0|is.na(TTDM)) %>% filter(!is.na(coh)) %>%
  filter( #Complete AnthropoAge (non-infinite) and mortality data
    !is.infinite(AnthropoAge), !is.na(AnthropoAge),
    !is.na(mortstat2), !is.na(WTRESP), WTRESP!=0) %>%
  filter( #Age 50-90, anthropometry within original NHANES-III ranges
    Age>=50, Age<90, AnthropoAge>=10, Height>=125, Height<=200,
    Weight>=30, Weight<=150, Waist>=50, Waist<=160, BMI>=15, BMI<=60)

# Education and lifestyle
#Education
HRS.A1$Education <- HRS.A1$raeducl %>%
  factor(1:3, c("Primary or less", "Secondary", "Tertiary"))
#Smoking
HRS.A1$Smoking <- HRS.A1 %>% with(case_when(
  SMOKEV==0~0, SMOKEN==0~1, SMOKEF<10~2, SMOKEF>=10~3)) %>%
  factor(0:3, c("Never smoker","Former smoker","<10/day",">=10/day"))
#Drinking
HRS.A1$Drinking <- HRS.A1 %>% with(case_when(
  DRINK==0~0, DRINKD==0~1, DRINKD%in%1:6~2, DRINKD>=7~3)) %>%
  factor(0:3, c("Never","Less than weekly","Less than daily","Daily"))
#Frequent physical activity (>1/week)
HRS.A1$VigPA1 <- with(HRS.A1, case_when(VGACTX%in%3:5~0, VGACTX%in%1:2~1))
HRS.A1$VigPA2 <- HRS.A1$VigPA1 %>% factor(0:1, c("No","Yes"))

# Set different baselines
HRS.A1_B1 <- HRS.A1 %>% filter(wave==8&INW8==1) %>%
  filter(!duplicated(ID)) #Only 2006 + remove follow-ups
HRS.A1_BX <- HRS.A1 %>% filter(!duplicated(ID)) #Remove follow-ups


#ELSA------------ ####
## ELSA --- 9 Waves: 2002-2018 (biennially)
## S-AnthropoAge can be calculated in W2,4,6 (2004, 2008, 2012)
ELSA <- readRDS("Databases/h_elsa_g2.rds.gz")

# Select variables and change to long format
ELSA.A1.0 <- ELSA %>% select(
  idauniqc, ragender, raracem, raeducl, radyear, r6iwstat, inw2, inw4, inw6,
  (c(2:6) %>% sapply(gsub, pattern="@", x=c("r@iwindm","r@iwindy")) %>%
     as.character), (c(2,4,6) %>% sapply(gsub, pattern="@", x=c(
    "r@cohort_e","r@agey", "r@mheight", "r@mweight", "r@mwaist",
    "r@shlt", "r@batha", "r@eata", "r@beda", "r@walkra", "r@toilta",
    "r@iadlfour", "r@hibpe", "r@diabe", "r@cancre", "r@lunge",
    "r@hearte", "r@stroke", "r@arthre","r@cwtresp",
    "r@smokev","r@smoken","r@smokef","r@drink","r@drinkd_e","r@vgactx_e")) %>%
      as.character)); ELSA.A1.0_names <- ELSA.A1.0 %>% select(-c(1:9)) %>%
  names; ELSA.A1.0 <- ELSA.A1.0 %>% pivot_longer(cols = -c(
  1:9, (ELSA.A1.0_names[!grepl(pattern="^r(\\d)(.*)", ELSA.A1.0_names)])),
  names_to = c("wave",".value"), names_pattern = "^r(\\d)(.*)")

ELSA.A1.0$eata <- with(ELSA.A1.0, case_when(eata==0~0, eata==1~1))
ELSA.A1.0$beda <- with(ELSA.A1.0, case_when(beda==0~0, beda==1~1))
ELSA.A1.0$batha <- with(ELSA.A1.0, case_when(batha==0~0, batha==1~1))
ELSA.A1.0$walkra <- with(ELSA.A1.0, case_when(walkra==0~0, walkra==1~1))
ELSA.A1.0$toilta <- with(ELSA.A1.0, case_when(toilta==0~0, toilta==1~1))

# Create new variables
ELSA.A1_W <- ELSA.A1.0 %>% filter(raracem==1) %>%
  mutate(#Ethnicity, anthropometry, health decline
  "ID"=idauniqc, "Sex"=factor(ragender,1:2, c("Men","Women")),
  "Ethnicity"=race_ethnicity1[1], "G2ASTUDY" = "ELSA",
  "mortstat"=r6iwstat, "Age"=agey, "Height"=mheight*100,
  "Weight"=mweight, "Waist"=mwaist, "B.IY"=iwindy, "B.IM"=iwindm,
  "SRH"=shlt, "ADL"=(batha + eata + beda + walkra + toilta),
  "IADL"=iadlfour, "HBP"=hibpe, "T2D"=diabe, "CAN"=cancre,
  "LUD"=lunge, "AMI"=hearte, "EVC"=stroke, "ART"=arthre,
  "BMI"=Weight/((Height/100)^2), "ICE"=Waist/Height) %>%
  mutate(#Date of baseline interview
    "B.DATE"=paste(B.IY,B.IM,"1",sep="-") %>% as.Date("%Y-%m-%d")
    #Date of last interview
    ); setDT(ELSA.A1_W)[, LAST.IW := max(
        B.DATE, na.rm = TRUE), by = ID]; ELSA.A1_W$LAST.IW[is.infinite(
        ELSA.A1_W$LAST.IW)] <- NA; ELSA.A1_W <- ELSA.A1_W %>%
  filter(wave%in%c(2,4,6)) %>%
  mutate(#Date of death
    "D.DATE"=paste(
      radyear,substr(LAST.IW,6,7),"1",sep="-") %>% as.Date("%Y-%m-%d"),
    "year"=factor(wave,c(2,4,6), c(2004,2008,2012)) %>% as.character %>%
      as.numeric); ELSA.A1_W <- ELSA.A1_W %>%
  mutate(#AnthropoAge and additional variables
    "AnthropoAge" = s_anthropoage_fast(ELSA.A1_W),
    "bi"=case_when(year%in%2003:2004~"2003/04", year%in%2007:2008~"2007/08",
                   year%in%2011:2012~"2011/12", year%in%2017:2018~"2017/18"),
    "participation"=(inw2==1)+(inw4==1)+(inw6==1),
    "consecutives"=(inw2==1&inw4==1)|(inw4==1&inw6==1),
    "coh"=case_when(cohort_e%in%c(1)~"A", cohort_e%in%c(4,8)~"B",
                    cohort_e%in%c(12)~"C"))

# Mortality status (0=Alive | 1=Dead)
ELSA.A1_W$mortstat2<-ifelse(ELSA.A1_W$mortstat %in% 5:6, 1, 0)
ELSA.A1_W$mortstat2 [ELSA.A1_W$mortstat==0] <- NA
# Censoring date: either date of death or date of last interview
ELSA.A1_W <- ELSA.A1_W %>% mutate("E.DATE" = ifelse(
  mortstat2==0|is.na(D.DATE), as.character(LAST.IW),
  as.character(D.DATE)) %>% as.Date)

# Filters
ELSA.A1 <- ELSA.A1_W %>%
  mutate( #Time to death in months + country specific race/ethnicity
    "TTDM"=interval(start=B.DATE,end=E.DATE) %/% months(1),
    "Ethnicity2"=paste0("UK-",race_ethnicity1[1])) %>%
  filter( #Complete follow-up data (≥0 months)
  TTDM>=0|is.na(TTDM)) %>%
  filter( #Complete AnthropoAge (non-infinite) and mortality data
    !is.infinite(AnthropoAge), !is.na(AnthropoAge),
    !is.na(cwtresp), cwtresp!=0) %>%
  filter( #Age 50-90, anthropometry within original NHANES-III ranges
    Age>=50, Age<90, AnthropoAge>=10, Height>=125, Height<=200,
    Weight>=30, Weight<=150, Waist>=50, Waist<=160, BMI>=15, BMI<=60)

# Education and lifestyle
#Education
ELSA.A1$Education <- ELSA.A1$raeducl %>%
  factor(1:3, c("Primary or less", "Secondary", "Tertiary"))
#Smoking
ELSA.A1$Smoking <- ELSA.A1 %>% with(case_when(
  smokev==0~0, smoken==0~1, smokef<10~2, smokef>=10~3)) %>%
  factor(0:3, c("Never smoker","Former smoker","<10/day",">=10/day"))
#Drinking
ELSA.A1$Drinking <- ELSA.A1 %>% with(case_when(
  drink==0~0, drinkd_e==0~1, drinkd_e%in%1:6~2, drinkd_e>=7~3)) %>%
  factor(0:3, c("Never","Less than weekly","Less than daily","Daily"))
#Frequent physical activity (>1/week)
ELSA.A1$VigPA1 <- with(ELSA.A1,case_when(vgactx_e%in%3:5~0,vgactx_e==2~1))
ELSA.A1$VigPA2 <- ELSA.A1$VigPA1 %>% factor(0:1, c("No","Yes"))

# Set different baselines
ELSA.A1_B1 <- ELSA.A1 %>% filter(wave==2&inw2==1) %>%
  filter(!duplicated(ID)) #Only 2004 + remove follow-ups
ELSA.A1_BX <- ELSA.A1 %>% filter(!duplicated(ID)) #Remove follow-ups


#MHAS------------ ####
#Harmonized data
MHAS <- readRDS("Databases/H_MHAS_c.rds.gz")

#2021 follow-up data
MHAS2021.h <- read_sav("Databases/MHAS/2021/MHAS_TRH_2021.sav") #Households
MHAS2021.c <- read_sav("Databases/MHAS/2021/MHAS_Core_2021.sav") #Individuals
MHAS2021.d <- read_sav("Databases/MHAS/2021/MHAS_EOL_2021.sav") #Next-of-kin
#Add interview dates
MHAS2021.hc <- MHAS2021.c %>% left_join(
  MHAS2021.h %>% filter(!duplicated(cunicah)) %>% select(
    cunicah,int_date_21), by = "cunicah"); MHAS2021.hc <- MHAS2021.hc %>%
  mutate(
    "int_date_21"=na_if(int_date_21, ""),
    "r6iwy"=MHAS2021.hc$int_date_21 %>% substr(7,8),
    "r6iwm"=MHAS2021.hc$int_date_21 %>% substr(4,5))
#Create new variables
MHAS2021.hc <- MHAS2021.hc %>% mutate( #ID = cunicah + NP
  "unhhidnp"=paste0(cunicah,np) %>% as.numeric, 
  "r6iwstat"=1, "rady21"=NA, "radm21"=NA) #Alive
MHAS2021.d <- MHAS2021.d %>% mutate( #ID = cunicah + NP
  "r6iwm"=NA,"r6iwy"=NA,"unhhidnp"=paste0(cunicah,np) %>% as.numeric, 
  "r6iwstat"=5, "rady21"=sa8_2_21, "radm21"=sa8_1_21) #Died this wave
#Select variables and rbind
MHAS2021 <- rbind(
  MHAS2021.hc %>% select(unhhidnp, r6iwy, r6iwm, rady21, radm21, r6iwstat),
  MHAS2021.d %>% select(unhhidnp, r6iwy, r6iwm, rady21, radm21, r6iwstat))
#Merge with harmonized data
MHAS_fin <- MHAS %>% left_join(MHAS2021, by = "unhhidnp")
#Interview status at wave 6
MHAS_fin$r6iwstat <- with( #1=Respondent, alive
  MHAS_fin, case_when(!is.na(r6iwstat)~r6iwstat, #5=Died this wave
                      is.na(r6iwstat)&r5iwstat%in%5:6~6, #6=Died prev. wave
                      is.na(r6iwstat)&!r5iwstat%in%5:6~9)) #9=Unknown
#Updated wave 6 interview year and month
MHAS_fin$r6iwy <- (MHAS_fin$r6iwy %>% as.numeric)+2000
MHAS_fin$r6iwm <- (MHAS_fin$r6iwm %>% as.numeric)
#Updated year and month of death
MHAS_fin$radyear <- with(MHAS_fin, ifelse((!is.na(rady21)),rady21,radyear))
MHAS_fin$radmonth <- with(MHAS_fin, ifelse((!is.na(radm21)),radm21,radmonth))

# Select variables and change to long format
MHAS.A1.0 <- MHAS_fin %>% select(
  unhhidnp, ragender, raeducl, radmonth, radyear,
  hacohort, r6iwstat, inw1, inw2, inw3,
  (1:6 %>% sapply(gsub, pattern="@", x=c("r@iwy","r@iwm")) %>% as.character),
  (1:3 %>% sapply(gsub, pattern="@", x=c(
    "r@agey", "r@mheight", "r@mweight", "r@mwaist", "r@shlt",
    "r@adla_m", "r@toilta", "r@iadlfour", "r@hibpe", "r@diabe", "r@cancre",
    "r@lunge_m", "r@hrtatte", "r@stroke", "r@arthre", "r@wtresp",
    "r@smokev","r@smoken","r@smokef","r@drink","r@drinkd","r@vigact")) %>%
     as.character)); MHAS.A1.0_names <- MHAS.A1.0 %>% select(-c(1:10)) %>%
  names; MHAS.A1.0 <- MHAS.A1.0 %>% pivot_longer(cols = -c(
    1:10, (MHAS.A1.0_names[!grepl(pattern="^r(\\d)(.*)", MHAS.A1.0_names)])),
    names_to = c("wave",".value"), names_pattern = "^r(\\d)(.*)")

# Create new variables
MHAS.A1_W <- MHAS.A1.0 %>%
  mutate( #Ethnicity, anthropometry, health decline
    "ID"=unhhidnp, "Sex"=factor(ragender,1:2, c("Men","Women")),
    "Ethnicity"=race_ethnicity1[3], "G2ASTUDY" = "MHAS",
    "mortstat"=r6iwstat, "Age"=agey, "Height"=mheight*100,
    "Weight"=mweight, "Waist"=mwaist, "B.IY"=iwy, "B.IM"=iwm, "SRH"=shlt,
    "ADL"=adla_m+toilta, #ADL: bath, eat, transfer, walk, toilet (0-5)
    "IADL"=iadlfour, "HBP"=hibpe, "T2D"=diabe, "CAN"=cancre, #IADL (0-4)
    "LUD"=lunge_m, "AMI"=hrtatte, "EVC"=stroke, "ART"=arthre,
    "BMI"=Weight/((Height/100)^2), "ICE"=Waist/Height) %>%
  mutate( #Date of baseline interview and date of death
    "B.DATE"=paste(B.IY,B.IM,"1",sep="-") %>% as.Date("%Y-%m-%d"),
    "D.DATE"=paste(radyear,radmonth,"1",sep="-")%>%as.Date("%Y-%m-%d"),
    "year"=factor(wave,1:5, c(2001,2003,2012,2015,2018)) %>% as.character %>%
      as.numeric); MHAS.A1_W<-MHAS.A1_W %>%
  mutate( #AnthropoAge and additional variables
    "AnthropoAge" = s_anthropoage_fast(MHAS.A1_W),
    "bi"=case_when(year%in%2001:2002~"2001/02", year%in%2003:2004~"2003/04",
                   year%in%2011:2012~"2011/12", year%in%2017:2018~"2017/18"),
    "participation"=(inw1==1)+(inw2==1)+(inw3==1),
    "consecutives"=(inw1==1&inw2==1)|(inw2==1&inw3==1),
    "coh"=case_when(hacohort%in%1~"A", hacohort%in%2~"B"))

# Date of last interview
setDT(MHAS.A1_W); MHAS.A1_W[, LAST.IW := max(B.DATE, na.rm = TRUE), by = ID]
MHAS.A1_W$LAST.IW[is.infinite(MHAS.A1_W$LAST.IW)] <- NA
# Mortality status (0=Alive | 1=Dead)
MHAS.A1_W$mortstat2<-ifelse(MHAS.A1_W$mortstat %in% 5:6, 1, 0)
MHAS.A1_W$mortstat2[MHAS.A1_W$mortstat==0] <- NA
# Censoring date: either date of death or date of last interview
MHAS.A1_W <- MHAS.A1_W %>% mutate("E.DATE" = ifelse(
  mortstat2==0|is.na(D.DATE), as.character(LAST.IW),
  as.character(D.DATE)) %>% as.Date)

# Filters
MHAS.A1 <- MHAS.A1_W %>%
  mutate( #Time to death in months + country specific race/ethnicity
    "TTDM"=interval(start=B.DATE,end=E.DATE) %/% months(1),
    "Ethnicity2"=paste0("MX-",race_ethnicity1[3])) %>%
  filter( #Complete follow-up data (≥0 months)
    wave%in%1:3) %>% filter(TTDM>=0|is.na(TTDM)) %>%
  filter( #Complete AnthropoAge (non-infinite) and mortality data
    !is.infinite(AnthropoAge), !is.na(AnthropoAge),
    !is.na(mortstat2), !is.na(B.DATE), !is.na(wtresp), wtresp!=0) %>%
  filter( #Age 50-90, anthropometry within original NHANES-III ranges
    Age>=50, Age<90, AnthropoAge>=10, Height>=125, Height<=200,
    Weight>=30, Weight<=150, Waist>=50, Waist<=160, BMI>=15, BMI<=60)

# Education and lifestyle
#Education
MHAS.A1$Education <- MHAS.A1$raeducl %>%
  factor(1:3, c("Primary or less", "Secondary", "Tertiary"))
#Smoking
MHAS.A1$Smoking <- MHAS.A1 %>% with(case_when(
  smokev==0~0, smoken==0~1, smokef<10~2, smokef>=10~3)) %>%
  factor(0:3, c("Never smoker","Former smoker","<10/day",">=10/day"))
#Drinking
MHAS.A1$Drinking <- MHAS.A1 %>% with(case_when(
  drink==0~0, drinkd==0~1, drinkd%in%1:6~2, drinkd>=7~3)) %>%
  factor(0:3, c("Never","Less than weekly","Less than daily","Daily"))
#Frequent physical activity (>1/week)
MHAS.A1$VigPA1 <- with(MHAS.A1,case_when(vigact==0~0,vigact==1~1))
MHAS.A1$VigPA2 <- MHAS.A1$VigPA1 %>% factor(0:1, c("No","Yes"))

# Set different baselines
MHAS.A1_B1 <- MHAS.A1 %>% filter(wave==1&inw1==1) %>%
  filter(!duplicated(ID)) #Only 2001 + remove follow-ups
MHAS.A1_B2 <- MHAS.A1 %>% filter(wave==3&inw3==1) %>%
  filter(!duplicated(ID)) #Only 2012 + remove follow-ups
MHAS.A1_BX <- MHAS.A1 %>% filter(!duplicated(ID)) #Remove follow-ups


#CRELES---------- ####
## CRELES ## -- 3 Waves: 2005 (1), 2007 (2), 2009 (3)
## S-AnthropoAge can be calculated in W1-W3
CRELES_W1 <- readRDS("Databases/CRELES/W1_CRELES.rds.gz")%>% distinct()
CRELES_W2 <- readRDS("Databases/CRELES/W2_CRELES.rds.gz")%>% distinct()
CRELES_W3 <- readRDS("Databases/CRELES/W3_CRELES.rds.gz")%>% distinct()
CRELES_deaths2 <- readRDS("Databases/CRELES/W2_CRELES_MORT.rds.gz")
CRELES_deaths3 <- readRDS("Databases/CRELES/W3_CRELES_MORT.rds.gz")
CRELES<- read_dta("Databases/H_CRELES.dta")%>% distinct()

# Wave 3: codes 997 and 998 are values missing not at random: waist/weight
# were seemingly not measured in "very obese" participants
CRELES_W3$k3[CRELES_W3$k3>900]<-NA; CRELES_W3$k3[CRELES_W3$k3==0]<-NA
CRELES_W3$k4[CRELES_W3$k4>900]<-NA; CRELES_W3$k4[CRELES_W3$k4==0]<-NA
CRELES_W3$k6[CRELES_W3$k6>900]<-NA; CRELES_W3$k6[CRELES_W3$k6==0]<-NA

# Extract anthropometry from CRELES studies
CRELES_W1a<-CRELES_W1 %>% select(idsujeto, k3b, k4b, k6b) %>% distinct() %>%
  mutate(k3b=measurements::conv_unit(k3b, "lbs", "kg"))%>%
  rename(r1pmwght=k3b, r1pmhght=k4b, r1pmwaist=k6b)
CRELES_W2a<-CRELES_W2 %>% select(idsujeto, k3b, k4b, k6b) %>% distinct() %>%
  mutate(k3b=measurements::conv_unit(k3b, "lbs", "kg"))%>%
  rename(r2pmwght=k3b, r2pmhght=k4b, r2pmwaist=k6b)
CRELES_W3a<-CRELES_W3 %>% select(idsujeto, k3, k4, k6) %>% distinct() %>%
  mutate(k3=measurements::conv_unit(k3, "lbs", "kg"))%>%
  rename(r3pmwght=k3, r3pmhght=k4, r3pmwaist=k6)

# Pooled mortality dataset
CRELES_D <- rbind(
  data.frame(
    "idsujeto"=CRELES_deaths2$idsujeto, "radyear"=CRELES_deaths2$sa4c,
    "radmonth"=CRELES_deaths2$sa4b),
  data.frame(
    "idsujeto"=CRELES_deaths3$idsujeto, "radyear"=CRELES_deaths3$sa4c,
    "radmonth"=CRELES_deaths3$sa4b)) %>% distinct()

# Add anthropometry and mortality to harmonized CRELES
CRELES <- CRELES %>%
  left_join(CRELES_W1a %>% select(
    idsujeto,r1pmwght,r1pmhght,r1pmwaist), by="idsujeto") %>% 
  left_join(CRELES_W2a %>% select(
    idsujeto,r2pmwght,r2pmhght,r2pmwaist), by="idsujeto") %>%
  left_join(CRELES_W3a %>% select(
    idsujeto,r3pmwght,r3pmhght,r3pmwaist), by="idsujeto") %>%
  left_join(CRELES_D, by="idsujeto")
# Mortality status (0=Alive | 1=Dead)
CRELES$mortstat<-ifelse(
  CRELES$r1iwstat%in%5:6 | CRELES$r2iwstat%in%5:6 | CRELES$r3iwstat%in%5:6,
  1,0); CRELES$mortstat[CRELES$mortstat==0]<-NA

# Select variables and change to long format
CRELES.A1.0 <- CRELES %>% select(
  idsujeto, ragender, raeducl, mortstat, radmonth, radyear, hacohort,
  inw1, inw2, inw3, (1:3 %>% sapply(gsub, pattern="@", x=c(
    "r@wtresp", "r@iwy", "r@iwm","r@agey", "r@iwstat", "r@shlt",
    "r@adla_cr", "r@iadla_cr","r@hibpe", "r@diabe", "r@cancre", "r@lunge",
    "r@hearte", "r@stroke", "r@arthre", "r@pmhght","r@pmwght","r@pmwaist",
    "r@smokev","r@smoken","r@smokef","r@drink","r@drinkd_cr","r@vigact")) %>%
      as.character)); CRELES.A1.0_names <- CRELES.A1.0 %>% select(-c(1:10)) %>%
  names; CRELES.A1.0 <- CRELES.A1.0 %>% pivot_longer(cols = -c(1:10,(
    CRELES.A1.0_names[!grepl(pattern="^r(\\d+)(.*)",CRELES.A1.0_names)])),
    names_to = c("wave",".value"), names_pattern = "^r(\\d+)(.*)")

# Create new variables
CRELES.A1_W <- CRELES.A1.0 %>%
  mutate(
    "ID"=idsujeto, "Sex"=factor(ragender,1:2, c("Men","Women")),
    "G2ASTUDY"="CRELES", "mortstat"=mortstat, "Ethnicity"=race_ethnicity1[3],
    "Age"=agey, "SRH"=shlt, "ADL"= adla_cr, "IADL"=iadla_cr, "Height"=pmhght,
    "Weight"=pmwght, "Waist"=pmwaist, "B.IY"=iwy, "B.IM"=iwm, "HBP"=hibpe,
    "T2D"=diabe, "CAN"=cancre, "LUD"=lunge, "AMI"=hearte, "EVC"=stroke,
    "ART"=arthre, "BMI"=Weight/((Height/100)^2), "ICE"=Waist/Height) %>%
  mutate(
    "B.DATE"=paste(B.IY,B.IM,"1",sep="-") %>% as.Date("%Y-%m-%d"),
    "D.DATE"=paste(radyear,radmonth,"1",sep="-")%>%as.Date("%Y-%m-%d"),
    "year"=factor(wave,1:3,c(2005,2007,2009)) %>% as.character %>%
      as.numeric); CRELES.A1_W<-CRELES.A1_W %>%
  mutate(
    "AnthropoAge" = s_anthropoage_fast(CRELES.A1_W),
    "bi"=case_when(wave==1~"2005/06", wave==2~"2007/08", wave==3~"2009/10"),
    "participation"=(inw1==1)+(inw2==1)+(inw3==1),
    "consecutives"=(inw1==1&inw2==1)|(inw2==1&inw3==1),
    "coh"=case_when(hacohort%in%1~"A", hacohort%in%2~"B")) %>%
  filter(hacohort==1)

# Date of last interview
setDT(CRELES.A1_W); CRELES.A1_W[, LAST.IW := max(B.DATE, na.rm = T), by = ID]
CRELES.A1_W$LAST.IW[is.infinite(CRELES.A1_W$LAST.IW)] <- NA
# Mortality status (0=Alive | 1=Dead)
CRELES.A1_W$mortstat2<-ifelse(CRELES.A1_W$mortstat %in% 1, 1, 0)
# Censoring date: either date of death or date of last interview
CRELES.A1_W <- CRELES.A1_W %>% mutate (E.DATE = ifelse(
  mortstat2==0|is.na(D.DATE), as.character(LAST.IW),
  as.character(D.DATE)) %>% as.Date)

# Filters
CRELES.A1 <- CRELES.A1_W %>%
  mutate( #Time to death in months + country specific race/ethnicity
    "TTDM"=interval(start=B.DATE,end=E.DATE) %/% months(1),
    "Ethnicity2"=paste0("CR-", race_ethnicity1[3])) %>%
  filter( #Complete follow-up data (≥0 months)
    wave%in%1:3) %>% filter(TTDM>=0|is.na(TTDM)) %>%
  filter( #Complete AnthropoAge (non-infinite) and mortality data
    !is.infinite(AnthropoAge), !is.na(AnthropoAge),
    !is.na(mortstat2), !is.na(B.DATE), !is.na(wtresp), wtresp!=0) %>%
  filter( #Age 50-90, anthropometry within original NHANES-III ranges
    Age>=50, Age<90, AnthropoAge>=10, Height>=125, Height<=200,
    Weight>=30, Weight<=150, Waist>=50, Waist<=160, BMI>=15, BMI<=60)

# Education and lifestyle
#Education
CRELES.A1$Education <- CRELES.A1$raeducl %>%
  factor(1:3, c("Primary or less", "Secondary", "Tertiary"))
#Smoking
CRELES.A1$Smoking <- CRELES.A1 %>% with(case_when(
  smokev==0~0, smoken==0~1, smokef<10~2, smokef>=10~3)) %>%
  factor(0:3, c("Never smoker","Former smoker","<10/day",">=10/day"))
#Drinking
CRELES.A1$Drinking <- CRELES.A1 %>% with(case_when(
  drinkd_cr==0~0, drinkd_cr==1~1, drinkd_cr==2~2, drinkd_cr==3~3)) %>%
  factor(0:3, c("Never","Less than weekly","Less than daily","Daily"))
#Frequent physical activity (>1/week)
CRELES.A1$VigPA1 <- with(CRELES.A1,case_when(vigact==0~0,vigact==1~1))
CRELES.A1$VigPA2 <- CRELES.A1$VigPA1 %>% factor(0:1, c("No","Yes"))

# Set different baselines
CRELES.A1_B1 <- CRELES.A1 %>% filter(wave==2&inw2==1) %>%
  filter(!duplicated(ID)) #Only 2007 + remove follow-ups
CRELES.A1_B2 <- CRELES.A1 %>% filter(wave==3&inw3==1) %>%
  filter(!duplicated(ID)) #Only 2009 + remove follow-ups
CRELES.A1_BX <- CRELES.A1 %>% filter(!duplicated(ID)) #Remove follow-ups


#CHARLS---------- ####
## CHARLS ## -- 5 Waves: 2011/12 (1), 2013/14 (2),
#                        2015/16 (3), 2018 (4), 2020 (5)
## S-AnthropoAge can be calculated only in W1-W3
#Harmonized data
CHARLS <- readRDS("Databases/H_CHARLS_D.rds.gz")

#2020 follow-up data
CHARLS20.c <- read_dta("Databases/CHARLS/2020/Demographic_Background.dta")#Core
CHARLS20.h <- read_dta("Databases/CHARLS/2020/Sample_Infor.dta")#Sample Info
CHARLS20.d <- read_dta("Databases/CHARLS/2020/Exit_Module.dta")#Next-of-kin
#Join status + interview date (died==0) + death date (died==1)
CHARLS20 <- CHARLS20.h %>% select(ID, died, iyear, imonth) %>% left_join(
  CHARLS20.d %>% select(ID, exb001_1, exb001_2, exb001_3, exb002),
  by="ID") %>% rename("r5iwy"=iyear, "r5iwm"=imonth) %>%
  mutate("r5iwy"=ifelse(died==1,NA,r5iwy) %>% as.numeric,
         "r5iwm"=ifelse(died==1,NA,r5iwm) %>% as.numeric,
         "r5iwstat"=ifelse(died==1,5,1))
#Convert dates from lunar to solar calendar
lunar <- CHARLS20 %>% filter(exb002==2); lunar <- lunar[c(1,5:7)] %>% na.omit
lunar$exb001_1 <- gsub("\\b(\\d)\\b", "0\\1", lunar$exb001_1)
lunar$exb001_2 <- gsub("\\b(\\d)\\b", "0\\1", lunar$exb001_2)
lunar$exb001_3 <- gsub("\\b(\\d)\\b", "0\\1", lunar$exb001_3)
lunar <- lunar %>% mutate("ld"=paste(exb001_1,exb001_2,exb001_3,sep="-"))
lunar2solar <- function(x){ #Vector of YYYY-MM-DD strings
  z<-c()
  for(i in 1:length(x)){
    a<-tryCatch(hongkong::as.lunar(x[i]), error=function(e)return(NA))
    b<-tryCatch(hongkong::lunarCal(a), error=function(e)return(NA))
    z<-c(z,b)}
  return(as.Date(z))
}
lunar$sd <- lunar2solar(lunar$ld)
lunar2 <- lunar %>% transmute(
  ID, "rady20"=lunar$sd %>% substr(1,4) %>% as.numeric(),
  "radm20"=lunar$sd %>% substr(6,7) %>% as.numeric())
#Add updated death year and death month to the full 2020 dataset
solar <- CHARLS20 %>% filter(exb002==1); solar <- solar[c(1,5:7)] %>% na.omit
solar2 <- solar %>% transmute(ID, "rady20"=exb001_1, "radm20"=exb001_2)
CHARLS20_fin <- CHARLS20 %>% left_join(rbind(solar2,lunar2), by="ID") %>%
  select(ID, r5iwy, r5iwm, rady20, radm20, r5iwstat)

#Merge 2020 follow-up with harmonized data
CHARLS_fin <- CHARLS %>% left_join(CHARLS20_fin, by = "ID")
#Interview status at wave 5
CHARLS_fin$r5iwstat <- with( #1=Respondent, alive
  CHARLS_fin, case_when(!is.na(r5iwstat)~r5iwstat, #5=Died this wave
                        is.na(r5iwstat)&r4iwstat%in%5:6~6, #6=Died prev. wave
                        is.na(r5iwstat)&!r4iwstat%in%5:6~9)) #9=Unknown
#Updated year and month of death
CHARLS_fin$radyear<-with(CHARLS_fin,ifelse((!is.na(rady20)),rady20,radyear))
CHARLS_fin$radmonth<-with(CHARLS_fin,ifelse((!is.na(radm20)),radm20,radmonth))

# Select variables and change to long format
CHARLS.A1.0 <- CHARLS_fin %>% select(
  ID, ragender, raeducl, hacohort_c, r5iwstat, radmonth, radyear, inw1, inw2,
  inw3, r4iwy, r4iwm, r5iwy, r5iwm, (1:3 %>% sapply(gsub, pattern="@", x=c(
    "r@iwy", "r@iwm","r@agey", "r@mheight", "r@mweight", "r@mwaist", "r@shlt",
    "r@adlfive", "r@moneya", "r@medsa", "r@shopa", "r@mealsa", "r@hibpe",
    "r@diabe", "r@cancre", "r@lunge", "r@hearte", "r@stroke", "r@arthre",
    "r@psyche", "r@dyslipe", "r@livere", "r@kidneye", "r@digeste", "r@wtresp",
    "r@smokev","r@smoken","r@smokef","r@drinkl","r@drinkn_c","r@vgactx_c")) %>%
      as.character)); CHARLS.A1.0_names <- CHARLS.A1.0 %>% select(-c(1:10)) %>%
  names; CHARLS.A1.0 <- CHARLS.A1.0 %>% pivot_longer(cols = -c(1:10,(
    CHARLS.A1.0_names[!grepl(pattern="^r(\\d)(.*)", CHARLS.A1.0_names)])),
    names_to = c("wave",".value"), names_pattern = "^r(\\d)(.*)")

# Create new variables
CHARLS.A1_W <- CHARLS.A1.0 %>%
  mutate( #Ethnicity, anthropometry, health decline
    "ID"=ID, "Sex"=factor(ragender,1:2, c("Men","Women")),
    "Ethnicity"=race_ethnicity1[4], "G2ASTUDY"="CHARLS", "mortstat"=r5iwstat,
    "Age"=agey, "Height"=mheight*100, "Weight"=mweight, "Waist"=mwaist,
    "B.IY"=iwy, "B.IM"=iwm, "HBP"=hibpe, "T2D"=diabe, "CAN"=cancre,
    "LUD"=lunge, "AMI"=hearte, "EVC"=stroke, "ART"=arthre,
    "SRH"=shlt, "ADL"=adlfive, #ADL (0-5)
    "IADL"=(moneya + medsa + shopa + mealsa), #IADL (0-4)
    "BMI"=Weight/((Height/100)^2), "ICE"=Waist/Height) %>%
  mutate( #Date of baseline interview and date of death
    "B.DATE"=paste(B.IY,B.IM,"1",sep="-") %>% as.Date("%Y-%m-%d"),
    "D.DATE"=paste(radyear,radmonth,"1",sep="-")%>%as.Date("%Y-%m-%d"),
    "year"=factor(wave,1:3, c(2011,2013,2015))%>%as.character%>%
      as.numeric); CHARLS.A1_W<-CHARLS.A1_W %>%
  mutate( #AnthropoAge and additional variables
    "AnthropoAge" = s_anthropoage_fast(CHARLS.A1_W),
    "bi"=case_when(year%in%2011:2012~"2011/12", year%in%2013:2014~"2013/14",
                   year%in%2015:2016~"2015/16", year%in%2017:2018~"2017/18"),
    "participation"=(inw1==1)+(inw2==1)+(inw3==1),
    "consecutives"=(inw1==1&inw2==1)|(inw2==1&inw3==1),
    "coh"=case_when(hacohort_c%in%1~"A", hacohort_c%in%2~"B",
                    hacohort_c%in%3~"C"))

# Mortality status and censoring
setDT(CHARLS.A1_W); CHARLS.A1_W[, LAST.IW := max(B.DATE, na.rm = T), by = ID]
CHARLS.A1_W$LAST.IW[is.infinite(CHARLS.A1_W$LAST.IW)] <- NA
# Mortality status (0=Alive | 1=Dead)
CHARLS.A1_W$mortstat2<-ifelse(CHARLS.A1_W$mortstat %in% 5:6, 1, 0)
CHARLS.A1_W$mortstat2[CHARLS.A1_W$mortstat==0] <- NA
# Censoring date: either date of death or date of last interview
CHARLS.A1_W <- CHARLS.A1_W %>% mutate("E.DATE" = ifelse(
  mortstat2==0|is.na(D.DATE), as.character(LAST.IW),
  as.character(D.DATE)) %>% as.Date)

# Filters
CHARLS.A1 <- CHARLS.A1_W %>%
  mutate( #Time to death in months + country specific race/ethnicity
    "TTDM"=interval(start=B.DATE,end=E.DATE) %/% months(1),
    "Ethnicity2"="CH-Other") %>%
  filter( #Complete follow-up data (≥0 months)
    TTDM>=0|is.na(TTDM)) %>% filter(!is.na(coh)) %>% 
  filter( #Complete AnthropoAge (non-infinite) and mortality data
    !is.infinite(AnthropoAge), !is.na(AnthropoAge),
    !is.na(mortstat2), !is.na(wtresp), wtresp!=0) %>%
  filter( #Age 50-90, anthropometry within original NHANES-III ranges
    Age>=50, Age<90, AnthropoAge>=10, Height>=125, Height<=200,
    Weight>=30, Weight<=150, Waist>=50, Waist<=160, BMI>=15, BMI<=60)

# Education and lifestyle
#Education
CHARLS.A1$Education <- CHARLS.A1$raeducl %>%
  factor(1:3, c("Primary or less", "Secondary", "Tertiary"))
#Smoking
CHARLS.A1$Smoking <- CHARLS.A1 %>% with(case_when(
  smokev==0~0, smoken==0~1, smokef<10~2, smokef>=10~3)) %>%
  factor(0:3, c("Never smoker","Former smoker","<10/day",">=10/day"))
#Drinking
CHARLS.A1$Drinking <- CHARLS.A1 %>% with(case_when(
  drinkn_c==0~0, drinkn_c%in%1:3~1, drinkn_c%in%4:6~2, drinkn_c%in%7:9~3)) %>%
  factor(0:3, c("Never","Less than weekly","Less than daily","Daily"))
#Frequent physical activity (>1/week)
CHARLS.A1$VigPA1 <- with(CHARLS.A1,case_when(vgactx_c<=3~0,vgactx_c>3~1))
CHARLS.A1$VigPA2 <- CHARLS.A1$VigPA1 %>% factor(0:1, c("No","Yes"))

# Set different baselines
CHARLS.A1_B1 <- CHARLS.A1 %>% filter(wave==1&inw1==1) %>%
  filter(!duplicated(ID)) #Only 2011/12 + remove follow-ups
CHARLS.A1_BX <- CHARLS.A1 %>% filter(!duplicated(ID)) #Remove follow-ups


#Pooled---------- ####
### BASELINE + FOLLOW-UPS ###
#Race/ethnicity recoding
race_ethnicity3 <- c(
  paste0("US-",race_ethnicity1)[1],
  paste0("UK-",race_ethnicity1[1]),
  paste0("US-",race_ethnicity1)[2:3],
  paste0(c("MX-","CR-"),race_ethnicity1[3]),
  paste0("US-",race_ethnicity1)[4], "CH-Other")
#Variable selection
mega_var <- c(
  "ID", "G2ASTUDY", "coh", "year", "wave", "bi", "participation",
  "consecutives", "Ethnicity", "Ethnicity2", "Sex", "Age", "mortstat2",
  "TTDM", "B.DATE", "D.DATE", "E.DATE", "AnthropoAge", "BMI","ICE","SRH",
  "ADL", "IADL", "HBP", "T2D", "CAN", "LUD", "AMI", "EVC", "ART",
  "Education","Smoking","Drinking","VigPA1","VigPA2")
HRS.A1 <- HRS.A1 %>% mutate(
  "baseline"=haven::labelled(ifelse(INW8+INW9>0, 1, 0))
  ); G2A_mega <- rbind(
  MHAS.A1 %>% select(all_of(mega_var), "Weights"=wtresp,"baseline"=inw1),
  HRS.A1 %>% select(all_of(mega_var), "Weights"=WTRESP, baseline),
  CHARLS.A1 %>% select(all_of(mega_var), "Weights"=wtresp, "baseline"=inw1),
  ELSA.A1 %>% select(all_of(mega_var), "Weights"=cwtresp, "baseline"=inw2),
  CRELES.A1 %>% select(all_of(mega_var), "Weights"=wtresp, "baseline"=inw1)
  ) %>% mutate("JID" = paste0(G2ASTUDY, ID),
         "Ethnicity2" = ordered(Ethnicity2, levels=race_ethnicity3))
#Comorbidities
health_var <- c(
  "SRH", "ADL", "IADL", "HBP", "T2D", "CAN", "LUD", "AMI", "EVC", "ART",
  "Education","Smoking","Drinking","VigPA1","VigPA2")
G2A_mega$n_comorb <- G2A_mega %>% select(
  all_of(health_var[4:10])) %>% apply(1,sum, na.rm=T)
#AnthropoAgeAccel
G2A_mega$AnthropoAgeAccel[
  with(G2A_mega,!is.na(AnthropoAge)&Sex=="Men")] <- lmer(
    AnthropoAge ~ Age + (1|G2ASTUDY) + (1|ID),
    G2A_mega %>% filter(Sex=="Men")) %>% residuals
G2A_mega$AnthropoAgeAccel[
  with(G2A_mega,!is.na(AnthropoAge)&Sex=="Women")] <- lmer(
    AnthropoAge ~ Age + (1|G2ASTUDY) + (1|ID),
    G2A_mega %>% filter(Sex=="Women")) %>% residuals
#New variables
G2A_mega <- G2A_mega %>% group_by(year) %>%
  mutate(
    "accelQ4"=quantcut(AnthropoAgeAccel,4) %>%
      factor(labels=c("Q1","Q2","Q3","Q4"))) %>% ungroup %>% 
  mutate(
    "accel"=ifelse(AnthropoAgeAccel>0,1,0),
    "comorb_cat"=cut(n_comorb, breaks=c(-Inf,0,1,Inf)) %>%
      ordered(labels=c("0","1",">=2")),
    "bi"=factor(bi, labels=c("01/02","03/04","05/06","07/08","09/10",
                             "11/12","13/14","15/16","17/18")),
    "Country"=ordered(G2ASTUDY, levels=c(
      "HRS","ELSA","MHAS","CRELES","CHARLS"), labels=c(
        "USA","England","Mexico","Costa Rica","China")),
    "Country3"=ordered(G2ASTUDY, levels=c(
      "HRS","ELSA","MHAS","CRELES","CHARLS"), labels=c(
        "EEUU","Inglaterra","México","Costa Rica","China")),
    "wave2" = as.numeric(wave), "TTDY" = TTDM/12, "Ethnicity3"=factor(
      Ethnicity2, ordered=F, levels=levels(Ethnicity2)))
#Standardized weights
G2A_mega<-G2A_mega %>% filter(Weights!=0) %>% group_by(Country, year) %>%
  mutate(std_weights=Weights*(n()/sum(Weights))) %>% ungroup()
#Time between interviews
setDT(G2A_mega); G2A_mega[,O.DATE := min(B.DATE, na.rm=T), by=ID] #First IW
G2A_mega[,Days := difftime(B.DATE, O.DATE, units="days"), by=ID] #Time in days
G2A_mega[, YearsDiff := (as.numeric(Days)/(365*2))] #Time in biennial cycles
G2A_mega[, YearsDiff2 := (as.numeric(Days)/(365))] #Time in years
#By Study
G2A_mega0 <- G2A_mega
G2A_mega1 <- filter(G2A_mega, G2ASTUDY=="HRS")
G2A_mega2 <- filter(G2A_mega, G2ASTUDY=="ELSA")
G2A_mega3 <- filter(G2A_mega, G2ASTUDY=="MHAS")
G2A_mega4 <- filter(G2A_mega, G2ASTUDY=="CRELES")
G2A_mega5 <- filter(G2A_mega, G2ASTUDY=="CHARLS")


### BASELINE ###
MHAS.A1_BX$n_comorb <- MHAS.A1_BX %>% select(
  all_of(health_var[4:10])) %>% apply(1,sum, na.rm=T)
HRS.A1_BX$n_comorb <- HRS.A1_BX %>% select(
  all_of(health_var[4:10])) %>% apply(1, sum, na.rm=T)
ELSA.A1_BX$n_comorb <- ELSA.A1_BX %>% select(
  all_of(health_var[4:10])) %>%
  apply(2, function(x){ifelse(x==1,1,0)}) %>% apply(1,sum, na.rm=T)
CRELES.A1_BX$n_comorb <- CRELES.A1_BX %>% select(
  all_of(health_var[4:10])) %>% apply(1,sum, na.rm=T)
CHARLS.A1_BX$n_comorb <- CHARLS.A1_BX %>% select(
  all_of(health_var[4:10])) %>% apply(1,sum, na.rm=T)
#Variable selection
pooled <- rbind(
  HRS.A1_BX %>% rename("Weights"=WTRESP) %>% select( #HRS
    ID, G2ASTUDY, mortstat2, AnthropoAge, Age, Sex, TTDM, BMI, ICE,
    Ethnicity, Ethnicity2, n_comorb, Weights, all_of(mega_var)),
  ELSA.A1_BX %>% rename("Weights"=cwtresp) %>% select( #ELSA
    ID, G2ASTUDY, mortstat2, AnthropoAge, Age, Sex, TTDM, BMI, ICE,
    Ethnicity, Ethnicity2, n_comorb, Weights, all_of(mega_var)),
  MHAS.A1_BX %>% rename("Weights"=wtresp) %>% select( #MHAS
    ID, G2ASTUDY, mortstat2, AnthropoAge, Age, Sex, TTDM, BMI, ICE,
    Ethnicity, Ethnicity2, n_comorb, Weights, all_of(mega_var)),
  CRELES.A1_BX %>% rename("Weights"=wtresp) %>% select( #CRELES
    ID, G2ASTUDY, mortstat2, AnthropoAge, Age, Sex, TTDM, BMI, ICE,
    Ethnicity, Ethnicity2, n_comorb, Weights, all_of(mega_var)),
  CHARLS.A1_BX %>% rename("Weights"=wtresp) %>% select( #CHARLS
    ID, G2ASTUDY, mortstat2, AnthropoAge, Age, Sex, TTDM, BMI, ICE,
    Ethnicity, Ethnicity2, n_comorb, Weights, all_of(mega_var))) %>%
  mutate("Ethnicity2" = ordered(Ethnicity2, levels=race_ethnicity3))
#AnthropoAgeAccel
#We do not need to add a random intercept for ID because this
#pooled dataset only contains cross-sectional surveys
pooled$AnthropoAgeAccel[
  with(pooled,!is.na(AnthropoAge)&Sex=="Men")] <- lmer(
    AnthropoAge ~ Age + (1|G2ASTUDY),# + (1|ID),
    pooled %>% filter(Sex=="Men")) %>% residuals
pooled$AnthropoAgeAccel[
  with(pooled,!is.na(AnthropoAge)&Sex=="Women")] <- lmer(
    AnthropoAge ~ Age + (1|G2ASTUDY),# + (1|ID),
    pooled %>% filter(Sex=="Women")) %>% residuals
#New variables + standardized weights
pooled<- pooled %>% filter(Weights!=0) %>%
  mutate(
    "accelQ4"=quantcut(AnthropoAgeAccel,4) %>%
      factor(labels=c("Q1","Q2","Q3","Q4")),
    "accel"=ifelse(AnthropoAgeAccel>=0,1,0),
    "comorb_cat"=cut(n_comorb, breaks=c(-Inf,0,1,Inf)) %>%
      ordered(labels=c("0","1",">=2")),
    "Ethnicity"= Ethnicity %>% factor(labels=c(
      "White","Black","Hispanic/Latino","Chinese")),
    "G2ASTUDY"=factor(G2ASTUDY, levels=c(
      "HRS","ELSA","MHAS","CRELES","CHARLS"))) %>% group_by(G2ASTUDY) %>%
  mutate(
    "std_weights"=Weights*(n()/sum(Weights))) %>% ungroup()
#By Study
pooled0 <- pooled
pooled1 <- filter(pooled, G2ASTUDY=="HRS")
pooled2 <- filter(pooled, G2ASTUDY=="ELSA")
pooled3 <- filter(pooled, G2ASTUDY=="MHAS")
pooled4 <- filter(pooled, G2ASTUDY=="CRELES")
pooled5 <- filter(pooled, G2ASTUDY=="CHARLS")
nrow(pooled); nrow(G2A_mega)


#Fig 1: Flowchart ####
##-------------------------------- HRS --------------------------------##
## Number of participants
HRS_alt %>% nrowf2 #1
HRS.A1_W %>% filter( 
  !is.na(AnthropoAge), !is.na(mortstat2), !is.na(WTRESP), WTRESP!=0) %>%
  filter(!duplicated(ID)) %>%  nrowf2 #2
HRS.A1 %>% filter(!duplicated(ID)) %>% nrowf2 #3
## Number of visits
HRS_alt %>% select( 
  (1:14 %>% sapply(gsub, pattern="@", x=c("R@AGEY_E")) %>% as.character)) %>%
  pivot_longer(cols=1:14, names_to = c("wave",".value"),
               names_pattern = "^R(\\d)(.*)") %>% nrowf2 #1
HRS.A1_W %>% filter(!is.na(AnthropoAge), !is.na(mortstat2),
                    !is.na(WTRESP), WTRESP!=0) %>%  nrowf2 #2
HRS.A1 %>% nrowf2 #3
## Deaths
(HRS.A1_BX$mortstat2 %>% table)[2] %>% format(big.mark=",")
(HRS.A1_BX$mortstat2 %>% table %>% prop.table)[2]
## Follow-up
with(HRS.A1_BX, fllwp.yrs(TTDM))
(with(HRS.A1_BX, sum(TTDM))/12) %>% floor


##-------------------------------- ELSA --------------------------------##
## Number of participants
ELSA %>% nrowf2 #1
ELSA.A1_W %>% filter(
  !is.na(AnthropoAge), !is.na(mortstat2), !is.na(cwtresp),
  cwtresp!=0) %>% filter(!duplicated(ID)) %>% nrowf2 #2
ELSA.A1 %>% filter(!duplicated(ID)) %>% nrowf2 #3
## Number of visits
ELSA %>% select(
  (1:9 %>% sapply(gsub, pattern="@", x=c("r@agey")) %>% as.character)) %>%
  pivot_longer(cols=1:9, names_to = c("wave",".value"),
               names_pattern = "^r(\\d)(.*)") %>% nrowf2 #1
ELSA.A1_W %>% filter(!is.na(AnthropoAge), !is.na(mortstat2),
                     !is.na(cwtresp), cwtresp!=0) %>%  nrowf2 #2
ELSA.A1 %>% nrowf2 #3
## Deaths
(ELSA.A1_BX$mortstat2 %>% table)[2] %>% format(big.mark=",")
## Follow-up
with(ELSA.A1_BX, fllwp.yrs(TTDM))
(with(ELSA.A1_BX, sum(TTDM))/12) %>% floor


##-------------------------------- MHAS --------------------------------##
## Number of participants
MHAS_fin %>% nrowf2 #1
MHAS.A1_W %>% filter( 
  !is.na(AnthropoAge), !is.na(mortstat2), !is.na(wtresp)) %>%
  filter(!duplicated(ID)) %>%  nrowf2 #2
MHAS.A1 %>% filter(!duplicated(ID)) %>% nrowf2 #3
## Number of visits
MHAS_fin %>% select(
  (1:5 %>% sapply(gsub, pattern="@", x=c("r@agey")) %>% as.character)) %>%
  pivot_longer(cols=1:5, names_to = c("wave",".value"),
               names_pattern = "^r(\\d)(.*)") %>% nrowf2 #1
MHAS.A1_W %>% filter(!is.na(AnthropoAge), !is.na(mortstat2),
                     !is.na(wtresp), wtresp!=0) %>%  nrowf2 #2
MHAS.A1 %>% nrowf2 #3
## Deaths
(MHAS.A1_BX$mortstat2 %>% table)[2] %>% format(big.mark=",")
## Follow-up
with(MHAS.A1_BX, fllwp.yrs(TTDM))
(with(MHAS.A1_BX, sum(TTDM))/12) %>% floor
MHAS.A1_BX$radyear %>% table
##--- 2012 subset ---##
MHAS.A1_B2 %>% nrowf2
(MHAS.A1_B2$mortstat2 %>% table)[2] %>% format(big.mark=",")
with(MHAS.A1_B2, fllwp.yrs(TTDM))


##-------------------------------- CRELES --------------------------------##
## Number of participants
CRELES %>% nrowf2 #1
CRELES.A1_W %>% filter(
  !is.na(AnthropoAge), !is.na(mortstat2), !is.na(wtresp),
  wtresp!=0) %>% filter(!duplicated(ID)) %>% nrowf2 #2
CRELES.A1 %>% filter(!duplicated(ID)) %>% nrowf2 #3
## Number of visits
CRELES %>% select(
  (1:4 %>% sapply(gsub, pattern="@", x=c("r@agey")) %>% as.character)) %>%
  pivot_longer(cols=1:4, names_to = c("wave",".value"),
               names_pattern = "^r(\\d)(.*)") %>% nrowf2 #1
CRELES.A1_W %>% filter(!is.na(AnthropoAge), !is.na(mortstat2),
                       !is.na(wtresp), wtresp!=0) %>%  nrowf2 #2
CRELES.A1 %>% nrowf2 #3
## Deaths
(CRELES.A1_BX$mortstat2 %>% table)[2] %>% format(big.mark=",")
## Follow-up
with(CRELES.A1_BX, fllwp.yrs(TTDM))
(with(CRELES.A1_BX, sum(TTDM))/12) %>% floor


##-------------------------------- CHARLS --------------------------------##
## Number of participants
CHARLS_fin %>% nrowf2 #1
CHARLS.A1_W %>% filter(
  !is.na(AnthropoAge), !is.na(mortstat2), !is.na(wtresp),
  wtresp!=0) %>% filter(!duplicated(ID)) %>% nrowf2 #2
CHARLS.A1 %>% filter(!duplicated(ID)) %>% nrowf2 #3
## Number of visits
CHARLS_fin %>% select(
  (1:4 %>% sapply(gsub, pattern="@", x=c("r@agey")) %>% as.character)) %>%
  pivot_longer(cols=1:4, names_to = c("wave",".value"),
               names_pattern = "^r(\\d)(.*)") %>% nrowf2 #1
CHARLS.A1_W %>% filter(!is.na(AnthropoAge), !is.na(mortstat2),
                       !is.na(wtresp), wtresp!=0) %>%  nrowf2 #2
CHARLS.A1 %>% nrowf2 #3
## Deaths
(CHARLS.A1_BX$mortstat2 %>% table)[2] %>% format(big.mark=",")
## Follow-up
with(CHARLS.A1_BX, fllwp.yrs(TTDM))
(with(CHARLS.A1_BX, sum(TTDM))/12) %>% floor

##-------------------------------- OVERALL --------------------------------##
#Overall sample
G2A_mega %>% filter(!duplicated(ID)) %>% nrowf2 #3
#Overal follow up
with(pooled, fllwp.yrs(TTDM))

#Participants followed-up from first wave with at least 2 data points
G2A_mega0 %>% filter(baseline==1, participation>1, !duplicated(ID)) %>% nrowf2
G2A_mega1 %>% filter(baseline==1, participation>1, !duplicated(ID)) %>% nrowf2
G2A_mega2 %>% filter(baseline==1, participation>1, !duplicated(ID)) %>% nrowf2
G2A_mega3 %>% filter(baseline==1, participation>1, !duplicated(ID)) %>% nrowf2
G2A_mega4 %>% filter(baseline==1, participation>1, !duplicated(ID)) %>% nrowf2
G2A_mega5 %>% filter(baseline==1, participation>1, !duplicated(ID)) %>% nrowf2

#STab1: Pop.table ####
#--- BASELINE POPULATION CHARACTERISTICS ---#
#Select variables
pooled_tab <- pooled %>% transmute(
  G2ASTUDY, Age, Sex, Ethnicity, Education, BMI, ICE, AnthropoAge,
  AnthropoAgeAccel, accel, Smoking, Drinking, VigPA1, ADL,IADL,SRH, T2D, HBP,
  AMI, CAN, LUD, EVC, ART, comorb_cat, mortstat2, "TTDY"=TTDM/12) %>% mutate(
    "T2D"=case_when(T2D==1~1, T2D==0~0), "HBP"=case_when(HBP==1~1, HBP==0~0),
    "AMI"=case_when(AMI==1~1, AMI==0~0), "CAN"=case_when(CAN==1~1, CAN==0~0),
    "LUD"=case_when(LUD==1~1, LUD==0~0), "EVC"=case_when(EVC==1~1, EVC==0~0),
    "ART"=case_when(ART==1~1, ART==0~0), "SRH"=case_when(
      SRH%in%1:5~as.numeric(SRH)), "G2ASTUDY"=ordered(
        G2ASTUDY, c("HRS","ELSA","MHAS","CRELES","CHARLS")))

#Set labels
labs1 <- paste0(names(pooled_tab[-1])); nlab <- labs1 %>% length
labs2 <- c(
  "Age (years)", "Sex (%)", "Race/ethnicity (%)", "Education level (%)",
  "BMI (kg/m^2)", "WHtR", "AnthropoAge (years)", "AnthropoAgeAccel (years)",
  "Accelerated aging (%)", "Smoking (%)", "Alcohol intake (%)",
  "Frequent vigorous activiy (%)", "ADL deficits (%)", "IADL deficits (%)",
  "Self-reported health (%)","Diabetes mellitus (%)",
  "Arterial hypertension (%)", "Myocardial infarction (%)",
  "Cancer (%)", "Chronic lung disease (%)", "Stroke (%)",
  "Arthritis (%)", "Number of comorbidities (%)",
  "Number of deahts (%)", "Follow-up time (years)")
for(i in 1:nlab){with(pooled_tab, setattr(get(labs1[i]), "label", labs2[i]))}

#Create table
tab1 <- tbl_summary(pooled_tab, by=G2ASTUDY, missing_text = "Missing") %>%
  add_p() %>% bold_labels() %>% add_overall() %>%
  as_flex_table() %>% align(align = "center",part = "all") %>% autofit()

set_flextable_defaults(
  split=F, table_align="center", table.layout="autofit"); read_docx() %>%
  body_add_flextable(value=tab1) %>% print(target = "Tables/TableS1.docx")

#Remove to save disk memory
remove(HRS_alt, MHAS, ELSA, CHARLS, CRELES)
remove(HRS.A1_W, MHAS.A1_W, ELSA.A1_W, CHARLS.A1_W, CRELES.A1_W)


####-------------------------------- METRIC PERFORMANCE -----------#### ----####
#Uno's c-statistic----- ####
strat.unoc <- function(x){
  #AnthropoAge
  a <- x %>% group_by(Q5) %>% summarise(
    "c"=concordance(
      Surv(TTDM, mortstat2) ~ AnthropoAge+strata(Sex)+strata(Ethnicity),
      reverse = TRUE,timewt = "n/G2", weights=std_weights)$concordance,
    "c_var"=concordance(
      Surv(TTDM, mortstat2) ~ AnthropoAge+strata(Sex)+strata(Ethnicity),
      reverse = TRUE,timewt = "n/G2", weights=std_weights)$var,
    "c_l"=c - sqrt(c_var)*qnorm(0.975), "c_u"=c + sqrt(c_var)*qnorm(0.975),
    "Par"="AnthropoAge") %>% select(Q5, c, c_l, c_u, Par)
  #Age
  b <- x %>% group_by(Q5) %>% summarise(
    "c"=concordance(
      Surv(TTDM, mortstat2) ~ Age+strata(Sex)+strata(Ethnicity),
      reverse = TRUE,timewt = "n/G2", weights=std_weights)$concordance,
    "c_var"=concordance(
      Surv(TTDM, mortstat2) ~ Age+strata(Sex)+strata(Ethnicity),
      reverse = TRUE,timewt = "n/G2", weights=std_weights)$var,
    "c_l"=c - sqrt(c_var)*qnorm(0.975), "c_u"=c + sqrt(c_var)*qnorm(0.975),
    "Par"="Age") %>% select(Q5, c, c_l, c_u, Par)
  rbind(a,b)}
evpergroup <- function(x,y){
  ev1 <- (x %>% select(Q5, mortstat2) %>% table)[,2]
  ev2 <- (x %>% select(Q5, mortstat2) %>% table %>% prop.table(1))[,2]
  ev <- paste0("Number of deaths per group: ", paste0(
    y,"=",ev1," (",round(ev2*100,1),"%)") %>% paste(collapse = ", "),".")
  ev}

## Overall
c1 <- with(pooled, concordance(
  Surv(TTDM, mortstat2) ~ AnthropoAge+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=std_weights, timewt = "n/G2"))
c2 <- with(pooled, concordance(
  Surv(TTDM, mortstat2) ~ Age+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=std_weights, timewt = "n/G2"))
c(c1$concordance, c1$concordance + (c(-1,1)*sqrt(c1$var)*qnorm(0.975)))
c(c2$concordance, c2$concordance + (c(-1,1)*sqrt(c2$var)*qnorm(0.975)))

## By G2A study
pooled %>% mutate("Q5"=G2ASTUDY) %>% strat.unoc %>%
  mutate("Q5"=reorder(Q5, -c)) %>% 
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=2) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1))+ labs(
    x="G2A study", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Fig2a; Fig2a
#Events per group
pooled %>% mutate("Q5"=G2ASTUDY) %>%
  evpergroup(y=levels(pooled$G2ASTUDY)) -> ev_G2A

## By Ethnicity
pooled %>% filter(Ethnicity2!="US-Other") %>%
  mutate("Q5"=Ethnicity) %>% strat.unoc %>%
  mutate("Q5"=reorder(Q5, -c)) %>% 
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="Race/ethnicity", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Fig2b; Fig2b
#Events per group
pooled %>% filter(Ethnicity2!="US-Other") %>% mutate("Q5"=Ethnicity) %>%
  evpergroup(y=levels(pooled$Ethnicity)) -> ev_Eth


#Time-dependent AUC---- ####
### Overall sample
pooled$years<-pooled$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=pooled, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=pooled, y=TRUE, x=TRUE)

follow <- 13; ROC.df.1 <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=pooled, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))

Fig2c <- autoplot(ROC.df.1, conf.int=T, ylim=c(0.7,0.8), lwd = 1.2) +
  theme_pubclean() + theme(legend.position = "top") +
  scale_color_manual(values = c("#9683ec","#35034f")) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)",
       x="Years of follow-up"); Fig2c

Fig2d <- as.data.frame(ROC.df.1$AUC$contrasts) %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) +
  labs(y="Delta tAUC of AnthropoAge - CA", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5); Fig2d


#Stratified c---------- ####
## By age quintiles
pooled %>% mutate("Q5"=Age %>% quantcut(seq(0,1,0.2))) %>% strat.unoc %>%
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="Age quintile (years)", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Figs1a; Figs1a
#Events per group
pooled %>% mutate("Q5"=Age %>% quantcut(seq(0,1,0.2))) %>%
  evpergroup(y=paste0("Q",1:5)) -> ev_Age

## By BMI quintiles
pooled %>% mutate("Q5"=BMI %>% quantcut(seq(0,1,0.2))) %>% strat.unoc %>%
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="BMI quintile (kg/m^2)", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Figs1b; Figs1b
#Events per group
pooled %>% mutate("Q5"=BMI %>% quantcut(seq(0,1,0.2))) %>%
  evpergroup(y=paste0("Q",1:5)) -> ev_BMI

## By WHtR quintiles
pooled %>% mutate("Q5"=ICE %>% round(2) %>% quantcut(seq(0,1,0.2))) %>%
  strat.unoc %>%
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="WHtR quintile", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Figs1c; Figs1c
#Events per group
pooled %>% mutate("Q5"=ICE %>% quantcut(seq(0,1,0.2))) %>%
  evpergroup(y=paste0("Q",1:5)) -> ev_ICE

## By Sex
pooled %>% mutate("Q5"=Sex) %>% strat.unoc %>%
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="Sex", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Figs1d; Figs1d
#Events per group
pooled %>% mutate("Q5"=Sex) %>% evpergroup(y=c("Men","Women")) -> ev_Sex

## By number of comorbidities
pooled %>% mutate("Q5"=n_comorb %>% cut(c(-Inf,0,1,Inf)) %>% factor(
  labels=c("0","1","≥2"))) %>% strat.unoc %>%
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="Number of comorbidities",
    y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> Figs1e; Figs1e
#Events per group
pooled %>% mutate("Q5"=n_comorb %>% cut(c(-Inf,0,1,Inf)) %>% factor(
  labels=c("0","1","≥2"))) %>% evpergroup(y=c("No","Yes")) -> ev_MMB

ev_Age
ev_BMI
ev_ICE
ev_Sex
ev_MMB


#Stratified tAUC------- ####
### HRS
HRS.A1_BX$years<-HRS.A1_BX$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=HRS.A1_BX, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=HRS.A1_BX, y=TRUE, x=TRUE)
follow <- 13; ROC.HRS <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=HRS.A1_BX, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.HRS<-as.data.frame(ROC.HRS$AUC$contrasts)

### ELSA
ELSA.A1_BX$years<-ELSA.A1_BX$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=ELSA.A1_BX, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=ELSA.A1_BX, y=TRUE, x=TRUE)
follow <- 13; ROC.ELSA <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=ELSA.A1_BX, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.ELSA<-as.data.frame(ROC.ELSA$AUC$contrasts)

### MHAS
MHAS.A1_BX$years<-MHAS.A1_BX$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_BX, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_BX, y=TRUE, x=TRUE)
follow <- 13; ROC.MHAS <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=MHAS.A1_BX, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.MHAS<-as.data.frame(ROC.MHAS$AUC$contrasts)

### CRELES
CRELES.A1_BX$years<-CRELES.A1_BX$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=CRELES.A1_BX, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=CRELES.A1_BX, y=TRUE, x=TRUE)
follow <- 13; ROC.CRELES <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=CRELES.A1_BX, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.CRELES<-as.data.frame(ROC.CRELES$AUC$contrasts)

### CHARLS
CHARLS.A1_BX$years<-CHARLS.A1_BX$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=CHARLS.A1_BX, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=CHARLS.A1_BX, y=TRUE, x=TRUE)
follow <- 13; ROC.CHARLS <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=CHARLS.A1_BX, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.CHARLS<-as.data.frame(ROC.CHARLS$AUC$contrasts)


#Fig 2: overall c+tAUC- ####
Fig2<-ggarrange(
  Fig2a, Fig2b, Fig2c, Fig2d, labels = LETTERS[1:4],
  nrow=2, ncol=2, common.legend = T)

ggsave(
  Fig2, file="Figures/Figure2.jpg", bg="transparent",
  width=25, height=20, units=c("cm"), dpi=600, limitsize=F)

#SFig1: c by quintiles- ####
Figs1 <- ggarrange(
  ggarrange(Figs1a, Figs1d, Figs1e, labels = LETTERS[1:3],
            ncol=3, nrow=1, common.legend = T),
  ggarrange(Figs1b, Figs1c, labels = LETTERS[4:5],
            ncol=2, nrow=1, legend = "none"), nrow=2)
ggsave(Figs1, file="Figures/FigureS1.jpg", bg="transparent",
       width=9*2.5, height=8*2, units=c("cm"), dpi=600, limitsize = FALSE)

#SFig2: tAUC by survey- ####
paste0("ROC.",c("HRS","ELSA","MHAS","CRELES","CHARLS")) %>% as.list %>%
  lapply(get) %>% lapply(function(x){x$AUC[["score"]][["lower"]] %>% min}) %>%
  as.numeric %>% min-0.01 -> min1
paste0("ROC.",c("HRS","ELSA","MHAS","CRELES","CHARLS")) %>% as.list %>%
  lapply(get) %>% lapply(function(x){x$AUC[["score"]][["upper"]] %>% max}) %>%
  as.numeric %>% max+0.01 -> max1
paste0("data.",c("HRS","ELSA","MHAS","CRELES","CHARLS")) %>% as.list %>%
  lapply(get) %>% lapply(function(x){x$lower %>% min}) %>%
  as.numeric %>% min-0.005 -> min2
paste0("data.",c("HRS","ELSA","MHAS","CRELES","CHARLS")) %>% as.list %>%
  lapply(get) %>% lapply(function(x){x$upper %>% max}) %>%
  as.numeric %>% max+0.005 -> max2
colors <- c("#9683ec","#35034f")

#HRS
FS2a_1 <- autoplot(ROC.HRS, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)",
       x="Years of follow-up", title="HRS") + theme_pubclean() + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FS2a_1
FS2a_2 <- data.HRS %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC (AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FS2a_2

#ELSA
FS2b_1 <- autoplot(ROC.ELSA, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)",
       x="Years of follow-up", title="ELSA") + theme_pubclean()  + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FS2b_1
FS2b_2 <- data.ELSA %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC (AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FS2b_2

#MHAS
FS2c_1 <- autoplot(ROC.MHAS, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)",
       x="Years of follow-up", title="MHAS") + theme_pubclean() + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FS2c_1
FS2c_2 <- data.MHAS %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC (AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FS2c_2

#CRELES
FS2d_1 <- autoplot(ROC.CRELES, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)",
       x="Years of follow-up", title="CRELES") + theme_pubclean() + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FS2d_1
FS2d_2 <- data.CRELES %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC (AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FS2d_2

#CHARLS
FS2e_1 <- autoplot(ROC.CHARLS, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)",
       x="Years of follow-up", title="CHARLS") + theme_pubclean() + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FS2e_1
FS2e_2 <- data.CHARLS %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC (AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FS2e_2

ROC.MHAS$times %>% max
ROC.HRS$times %>% max
ROC.CHARLS$times %>% max
ROC.ELSA$times %>% max
ROC.CRELES$times %>% max

ggarrange(
  ggarrange(FS2c_1,FS2a_1,FS2e_1,FS2b_1,FS2d_1, nrow=1, ncol=5, 
            common.legend = T, legend="bottom", labels=LETTERS[1:5]),
  ggarrange(FS2c_2,FS2a_2,FS2e_2,FS2b_2,FS2d_2, nrow=1, ncol=5, 
            labels=LETTERS[6:10]), nrow=2, ncol=1) -> FigS2
ggsave(
  FigS2, file="Figures/FigureS2.jpg", bg="transparent",
  width=40, height=23.5, units=c("cm"), dpi=600, limitsize = FALSE)


#SFig3: c+tAUC in MHAS- ####
## Overall MHAS
#Uno's c-statistic
c1.1 <- with(MHAS.A1_BX, concordance(
  Surv(TTDM, mortstat2) ~ Age+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=wtresp, timewt = "n/G2"))
c1.2 <- with(MHAS.A1_BX, concordance(
  Surv(TTDM, mortstat2) ~ AnthropoAge+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=wtresp, timewt = "n/G2"))
#t-AUC
MHAS.A1_BX$years<-MHAS.A1_BX$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_BX, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_BX, y=TRUE, x=TRUE)
follow <- 13; ROC.MHAS1 <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=MHAS.A1_BX, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.MHAS1<-as.data.frame(ROC.MHAS1$AUC$contrasts)

## MHAS - 2001 cohort
#Uno's c-statistic
c2.1 <- with(MHAS.A1_B1, concordance(
  Surv(TTDM, mortstat2) ~ Age+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=wtresp, timewt = "n/G2"))
c2.2 <- with(MHAS.A1_B1, concordance(
  Surv(TTDM, mortstat2) ~ AnthropoAge+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=wtresp, timewt = "n/G2"))
#tAUC
MHAS.A1_B1$years<-MHAS.A1_B1$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_B1, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_B1, y=TRUE, x=TRUE)
follow <- 13; ROC.MHAS2 <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=MHAS.A1_B1, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.MHAS2<-as.data.frame(ROC.MHAS2$AUC$contrasts)

## MHAS - 2012 cohort
c3.1 <- with(MHAS.A1_B2, concordance(
  Surv(TTDM, mortstat2) ~ Age+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=wtresp, timewt = "n/G2"))
c3.2 <- with(MHAS.A1_B2, concordance(
  Surv(TTDM, mortstat2) ~ AnthropoAge+strata(Sex)+strata(Ethnicity),
  reverse = T, weights=wtresp, timewt = "n/G2"))
#tAUC
MHAS.A1_B2$years<-MHAS.A1_B2$TTDM/12
m1 <- coxph(Surv(years,mortstat2)~ Age+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_B2, y=TRUE, x=TRUE)
m2 <- coxph(Surv(years,mortstat2)~ AnthropoAge+strata(Sex)+
              strata(Ethnicity), data=MHAS.A1_B2, y=TRUE, x=TRUE)
follow <- 9; ROC.MHAS3 <- Score(
  list("Age"=m1, "AnthropoAge"=m2), formula=Hist(years,mortstat2)~1,
  data=MHAS.A1_B2, conf.int=T, cens.method = "ipcw", summary="risks",
  metrics="auc", plots="roc", times = c(seq(1, follow, 1)))
data.MHAS3<-as.data.frame(ROC.MHAS3$AUC$contrasts)


### Plots
#Uno's c statistic
paste0("c",c(1,1,2,2,3,3),".",rep(1:2,3)) %>%
  as.list %>% lapply(get) %>% sapply(function(x){c(
    x$concordance, x$concordance + (c(-1,1)*sqrt(x$var)*qnorm(0.975)))}) %>%
  t %>% as.data.frame() %>% `names<-`(c("c","c_l","c_u")) %>% mutate(
    "Q5"=c(1,1,2,2,3,3) %>% factor(labels=c(
      "Overall\nMHAS","2001\ncohort", "2012\ncohort")),
    "Par"=rep(1:2,3) %>% factor(labels=c("Age","AnthropoAge"))) %>%
  ggplot(aes(x=Q5, y=c, group=Par, fill=Par, ymin=c_l, ymax=c_u)) +
  geom_crossbar(position = position_dodge2(), width=0.3, alpha=0.85) +
  scale_fill_manual(values = c("#9683ec","#35034f")) +
  theme_pubclean() + geom_hline(yintercept = 0.5, linetype=5) +
  scale_y_continuous(breaks=seq(.2,1,0.15), limits = c(0.3,1)) + labs(
    x="MHAS cohort", y="c-statistic (95% CI)", fill="Predictor") + 
  theme(legend.position = "top", plot.title = element_text(
    hjust=0.5, face="bold")) -> FigS3a; FigS3a

#tAUC
(c(ROC.MHAS2$AUC$score$lower, ROC.MHAS3$AUC$score$lower) %>%
    min-0.01) %>% round(2) -> min1
(c(ROC.MHAS2$AUC$score$upper, ROC.MHAS3$AUC$score$upper) %>%
    max+0.01) %>% round(2) -> max1
FigS3b <- autoplot(ROC.MHAS2, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)\n",
       x="Years of follow-up", title="MHAS 2001") + theme_pubclean() + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FigS3b
FigS3c <- autoplot(ROC.MHAS3, conf.int=T, lwd = 1.2) + ylim(min1,max1) +
  scale_color_manual(values=colors) + scale_fill_manual(values=colors) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(fill="Marker", col="Marker", y="tAUC (95%CI)\n",
       x="Years of follow-up", title="MHAS 2012") + theme_pubclean() + 
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position = "top"); FigS3c

#Delta tAUC
(c(data.MHAS2$lower,data.MHAS3$lower) %>% min-0.005) %>% round(2) -> min2 
(c(data.MHAS2$upper,data.MHAS3$upper) %>% max+0.005) %>% round(2) -> max2
FigS3d <- data.MHAS2 %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC\n(AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FigS3d
FigS3e <- data.MHAS3 %>%
  ggplot(aes(y=delta.AUC, x=times)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.2) + ylim(min2,max2) +
  scale_x_continuous(breaks=seq(0,14,2)) +
  labs(y="Delta tAUC\n(AnthropoAge - CA)", x="Years of follow-up")+
  theme_pubclean() + geom_hline(yintercept = 0, linetype=5) + theme(
    plot.title=element_text(hjust=0.5, face="bold", size=14)); FigS3e

FigS3 <- ggarrange(
  ggarrange(FigS3a, labels="A"), ggarrange(
    ggarrange(FigS3b,FigS3c,nrow=1,ncol=2,labels=LETTERS[2:3],
              common.legend = T, legend="bottom"),
    ggarrange(FigS3d,FigS3e,nrow=1,ncol=2,labels=LETTERS[4:5]),
    nrow=2, ncol=1, heights = c(1,0.6)), ncol=2, widths = c(0.6, 1))
ggsave(
  FigS3, file="Figures/FigureS3.jpg", bg="transparent",
  width=36, height=18.5, units=c("cm"), dpi=600, limitsize = FALSE)

#STab2: Harrel's c----- ####
##Harrell's c-statistic  comparison
compareC_alt <- function(x,y){
  #Load packages
  pacman::p_load(BiocManager)
  if (!requireNamespace("survcomp", quietly = TRUE)) {
    BiocManager::install("survcomp", ask = FALSE)}
  #C statistics
  index1 <- with(x, survcomp::concordance.index(
    x=AnthropoAge, surv.time=TTDM, surv.event=mortstat2, outx = F,
    weights=std_weights, strat=paste(Sex, Ethnicity), method="noether"))
  index2 <- with(x, survcomp::concordance.index(
    x=Age, surv.time=TTDM, surv.event=mortstat2, outx = F,
    weights=std_weights, strat=paste(Sex, Ethnicity), method="noether"))
  #Compare
  comp <- survcomp::cindex.comp(index1, index2)
  #Data.frame
  a <- paste0(sprintf(index1[[1]], fmt="%#.3f"), " (",
              sprintf(index1[[3]], fmt="%#.3f"), "-",
              sprintf(index1[[4]], fmt="%#.3f"),")")
  b <- paste0(sprintf(index2[[1]], fmt="%#.3f"), " (",
              sprintf(index2[[3]], fmt="%#.3f"), "-",
              sprintf(index2[[4]], fmt="%#.3f"),")")
  c <- (index1$c.index - index2$c.index )%>% sprintf(fmt="%#.4f")
  d <- comp$p.value %>% format.pval(digits=3, eps=c(0.001))
  options(scipen=10)
  data.frame(y,a,b,c,d) %>% `names<-`(c(
    "Strata", "AnthropoAge\nc-statistic (95% CI)",
    "Age\nc-statistic (95% CI)", "Difference", "P-value"))}

var1 <- pooled$Ethnicity %>% table %>% names; rbind(
  compareC_alt(pooled, "Pooled"),
  c("By G2A study","","","",""),
  compareC_alt(pooled %>% filter(G2ASTUDY=="HRS"), "HRS"),
  compareC_alt(pooled %>% filter(G2ASTUDY=="ELSA"), "ELSA"),
  compareC_alt(pooled %>% filter(G2ASTUDY=="MHAS"), "MHAS"),
  compareC_alt(pooled %>% filter(G2ASTUDY=="CRELES"), "CRELES"),
  compareC_alt(pooled %>% filter(G2ASTUDY=="CHARLS"), "CHARLS"),
  c("By race/ethnicity","","","",""),
  compareC_alt(pooled %>% filter(Ethnicity==var1[1]), var1[1]),
  compareC_alt(pooled %>% filter(Ethnicity==var1[2]), var1[2]),
  compareC_alt(pooled %>% filter(Ethnicity==var1[3]), var1[3]),
  compareC_alt(pooled %>% filter(Ethnicity==var1[4]), var1[4]),
  c("By sex","","","",""),
  compareC_alt(pooled %>% filter(Sex=="Men"), "Men"),
  compareC_alt(pooled %>% filter(Sex=="Women"), "Women"),
  c("By number of comorbidities","","","",""),
  compareC_alt(pooled %>% filter(n_comorb==0), "0"),
  compareC_alt(pooled %>% filter(n_comorb==1), "1"),
  compareC_alt(pooled %>% filter(n_comorb>=2), "≥2")) -> tab.conc1

tab.conc1 %>% flextable %>%
  merge_at(i=2, j=1:5) %>% merge_at(i=8, j=1:5) %>%
  merge_at(i=13, j=1:5) %>% merge_at(i=16, j=1:5) %>% 
  bold(part = "header") %>% bold(i=c(2,8,13,16), j=1:5) %>% 
  italic(i=c(2,8,13,16), j=1:5) %>%
  align(align="center", part="all") %>% autofit() -> tab2

set_flextable_defaults(
  split=F, table_align="center", table.layout="autofit"); read_docx() %>%
  body_add_flextable(value=tab2) %>% print(target = "Tables/TableS2.docx")


####-------------------------------- CHANGES OVER TIME ------------#### ----####
#Correlation structure- ####
#HRS 8, 10, 12, 14
G2A_mega1 %>% select(ID, AnthropoAge, year) %>% arrange(year) %>%
  pivot_wider(names_from=year, values_from=AnthropoAge, values_fill=NA) %>%
  select(2,4,6,8) %>% cor(use="pairwise.complete.obs")
#HRS 9, 11, 13
G2A_mega1 %>% select(ID, AnthropoAge, year) %>% arrange(year) %>%
  pivot_wider(names_from=year, values_from=AnthropoAge, values_fill=NA) %>%
  select(3,5,7) %>% cor(use="pairwise.complete.obs")
#ELSA
G2A_mega2 %>% select(ID, AnthropoAge, year) %>% arrange(year) %>%
  pivot_wider(names_from=year, values_from=AnthropoAge, values_fill=NA) %>%
  select(-1) %>% cor(use="pairwise.complete.obs")
#MHAS
G2A_mega3 %>% select(ID, AnthropoAge, year) %>% arrange(year) %>%
  pivot_wider(names_from=year, values_from=AnthropoAge, values_fill=NA) %>%
  select(-1) %>% cor(use="pairwise.complete.obs")
#CRELES
G2A_mega4 %>% select(ID, AnthropoAge, year) %>% arrange(year) %>%
  pivot_wider(names_from=year, values_from=AnthropoAge, values_fill=NA) %>%
  select(-1) %>% cor(use="pairwise.complete.obs")
#CHARLS
G2A_mega5 %>% select(ID, AnthropoAge, year) %>% arrange(year) %>%
  pivot_wider(names_from=year, values_from=AnthropoAge, values_fill=NA) %>%
  select(-1) %>% cor(use="pairwise.complete.obs")


#Gaussian GEE's ------- ####
mBS1 <- geeglm(
  AnthropoAge~YearsDiff2, id=ID, data=filter(G2A_mega, Country=="USA"),
  weights = std_weights, wave = wave2, corstr = "ar1", family = gaussian())
mBS2 <- geeglm(
  AnthropoAge~YearsDiff2, id=ID, data=filter(G2A_mega, Country=="England"),
  weights = std_weights, wave = wave2, corstr = "ar1", family = gaussian())
mBS3 <- geeglm(
  AnthropoAge~YearsDiff2, id=ID, data=filter(G2A_mega, Country=="Mexico"),
  weights = std_weights, wave = wave2, corstr = "ar1", family = gaussian())
mBS4 <- geeglm(
  AnthropoAge~YearsDiff2, id=ID, data=filter(G2A_mega, Country=="Costa Rica"),
  weights = std_weights, wave = wave2, corstr = "ar1", family = gaussian())
mBS5 <- geeglm(
  AnthropoAge~YearsDiff2, id=ID, data=filter(G2A_mega, Country=="China"),
  weights = std_weights, wave = wave2, corstr = "ar1", family = gaussian())

b1 <- summary(mBS1)$coefficients[2,1]; ci1 <- confint(mBS1) %>% round(3)
b2 <- summary(mBS2)$coefficients[2,1]; ci2 <- confint(mBS2) %>% round(3)
b3 <- summary(mBS3)$coefficients[2,1]; ci3 <- confint(mBS3) %>% round(3)
b4 <- summary(mBS4)$coefficients[2,1]; ci4 <- confint(mBS4) %>% round(3)
b5 <- summary(mBS5)$coefficients[2,1]; ci5 <- confint(mBS5) %>% round(3)

ann_f1 <- c(
  paste0("β = ", b1 %>% round(3), "\n(", ci1[2,1], "—", ci1[2,2], ")"),
  paste0("β = ", b2 %>% round(3), "\n(", ci2[2,1], "—", ci2[2,2], ")"),
  paste0("β = ", b3 %>% round(3), "\n(", ci3[2,1], "—", ci3[2,2], ")"),
  paste0("β = ", b4 %>% round(3), "\n(", ci4[2,1], "—", ci4[2,2], ")"),
  paste0("β = ", b5 %>% round(3), "\n(", ci5[2,1], "—", ci5[2,2], ")"))

#Fig 3: Aging trends--- ####
#CA vs BA: Longitudinal trends by year and country
set.data.plot <- function(x){
  a <- svyby(
    ~AnthropoAge, by=~bi+Country+G2ASTUDY, design=svydesign(
      data=x, ids=~JID, weights=~std_weights), FUN=svymean) %>%
    `rownames<-`(NULL) %>% rename("mean"=AnthropoAge) %>% mutate(
      "Par"="AnthropoAge", "l.ci"=mean-qnorm(0.975)*se,
      "u.ci"=mean+qnorm(0.975)*se)
  b <- svyby(
    ~Age, by=~bi+Country+G2ASTUDY, design=svydesign(
      data=x, ids=~JID, weights=~std_weights), FUN=svymean) %>%
    `rownames<-`(NULL) %>% rename("mean"=Age) %>% mutate(
      "Par"="Age", "l.ci"=mean-qnorm(0.975)*se,
      "u.ci"=mean+qnorm(0.975)*se)
  rbind(a,b)}
country1 <- c("USA","England","Mexico","Costa Rica","China")
country2 <- c("EEUU","Inglaterra","México","Costa Rica","China")
data.plot_fin <- rbind(
  set.data.plot(G2A_mega1 %>% filter(baseline==1, participation>1)),
  set.data.plot(G2A_mega2 %>% filter(baseline==1, participation>1)),
  set.data.plot(G2A_mega3 %>% filter(baseline==1, participation>1)),
  set.data.plot(G2A_mega4 %>% filter(baseline==1, participation>1)),
  set.data.plot(G2A_mega5 %>% filter(baseline==1, participation>1))) %>% mutate(
    "Country"=ordered(Country, levels=country1),
    "Country2"=ordered(Country, labels=country2),
    "Study"=ordered(G2ASTUDY, c("HRS","ELSA","MHAS","CRELES","CHARLS")),
    "ann"=case_when(
      Country=="USA"~ann_f1[1], Country=="England"~ann_f1[2],
      Country=="Mexico"~ann_f1[3], Country=="Costa Rica"~ann_f1[4],
      Country=="China"~ann_f1[5]),
    "nobs"=case_when(
      Country=="USA"~table(Country)[1],
      Country=="England"~table(Country)[2],
      Country=="Mexico"~table(Country)[3],
      Country=="Costa Rica"~table(Country)[4],
      Country=="China"~table(Country)[5]), "ycoord"=65
  ); data.plot_fin %>% 
  ggplot(aes(x=bi, y=mean, ymin=l.ci, ymax=u.ci, color=Par, group=Par)) +
  geom_line() + geom_pointrange() + geom_point() + theme_pubclean() +
  facet_wrap(~Study, ncol=5, nrow=1) +
  scale_color_manual(values=c("#9683ec", "#35034f"))  + scale_x_discrete(
    breaks=c("01/02","05/06","09/10","13/14","17/18"),
    labels=c("2001/02","2005/06","2009/10","2013/14","2017/18")) +
  labs(x=NULL, y="Weighted mean (years)\nwith 95%CI", color=NULL, group=NULL)+
  theme(
    legend.position = "bottom", legend.key.size = unit(1,"cm"),
    legend.title = element_text(hjust=0.5, face = "bold", size=12),
    legend.text = element_text(hjust=0.5, size=11),
    strip.text = element_text(face="bold.italic", vjust=-1, size=12),
    panel.spacing.x = unit(2.5,"lines"), strip.background = element_blank(),
    axis.text.x = element_text(vjust=0.7, hjust=0.8, size=12, angle=45),
    axis.text.y = element_text(size=12), axis.title=element_text(size=13)) +
  geom_text(aes(label = ifelse(nobs[nobs] > 2, ann, ""), y=ycoord),
            color="black", x="13/14", size=4.25) -> F2A; F2A

#Accelerated aging by year and country
data.p2 <- G2A_mega %>% filter(baseline==1, participation>1) %>%
  transmute(JID, Ct=Country, accel, yr=year, bi, wt=std_weights) %>%
  filter((Ct=="USA"&yr%in%(2e3+6*1:3))|(Ct!="USA"&yr%in%(2e3+1:15))) %>%
  mutate(accel=accel*100); svyby(~accel, by=~bi+Ct, design=svydesign(
           data=data.p2, ids=~JID, weights=~wt), FUN=svymean) %>%
  `rownames<-`(NULL) %>% rename("mean"=accel) %>% mutate(
    "l.ci"=mean-qnorm(0.975)*se, "u.ci"=mean+qnorm(0.975)*se) %>%
  ggplot(aes(x=bi, y=mean)) + geom_errorbar(
    aes(ymin = l.ci, ymax = u.ci), width=0.35, linewidth=0.85) +
  geom_col(fill="#35034f") + facet_wrap(~Ct, ncol=5, nrow=1, "free_x") +
  theme_pubclean() + xlab(NULL) + ylab("Accelerated\naging (%)") +
  ggbreak::scale_y_break(breaks = c(1,15), space=0.05) +
  scale_y_continuous(breaks = c(0,15,30,45,60)) + theme(
    axis.line.y.right=element_blank(), axis.ticks.y.right=element_blank(),
    axis.text.y.right=element_blank(), axis.title=element_text(size=13),
    axis.text.x=element_text(vjust=0.85, size=12), legend.position = "none",
    axis.text.y=element_text(size=12), strip.text = element_blank(),
    strip.background = element_blank()) + geom_text(
      aes(label=round(mean,1), x=bi, y=mean-4.5),
      color="white", size=5, fontface="bold.italic") -> F2B; F2B

fig3<-ggarrange(F2A, print(F2B), ncol=1, labels=LETTERS[1:2],
                nrow=2, heights = c(1.25,0.5)) #1500x750pix
ggsave(fig3, file="Figures/Figure3.jpg", bg="transparent",
       width=35, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

#CA vs BA: Scatter plots by year and country
G2A_mega %>% mutate(
  Country=ordered(G2ASTUDY, levels=c(
    "MHAS","HRS","ELSA","CHARLS", "CRELES"))) %>% 
  ggplot(aes(x=Age, y=AnthropoAge, color=Country, linetype=Country)) +
  geom_jitter(color="gray70", alpha=0.15) +
  geom_smooth(method = "lm") + facet_wrap(~bi, ncol=3, nrow=3) + 
  theme_pubclean() + scale_color_manual(values=c(
    "#280659", "#9d0610", "#f18701", "#ae60d3","#ff60d3")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold.italic", vjust=-1, size=12),
        strip.background=element_blank(), plot.background=element_blank(),
        panel.background=element_blank(), legend.background=element_blank())

#SFig4: Trends by sex-- ####
gee_gauss_sex <- function(x){
  geeglm1 = geeglm(
    AnthropoAge~YearsDiff2, id = ID, data = x %>% filter(Sex=="Women"),
    weights = std_weights, wave = wave2, corstr = "ar1",
    std.err = "san.se", family = gaussian())
  geeglm2 = geeglm(
    AnthropoAge~YearsDiff2, id = ID, data = x %>% filter(Sex=="Men"),
    weights = std_weights, wave = wave2, corstr = "ar1",
    std.err = "san.se", family = gaussian())
  return(list(geeglm1, geeglm2))}
mBS1_ <- gee_gauss_sex(G2A_mega %>% filter(
  G2ASTUDY=="HRS")); mBS1_w <- mBS1_[[1]]; mBS1_m <- mBS1_[[2]]
mBS2_ <- gee_gauss_sex(G2A_mega %>% filter(
  G2ASTUDY=="ELSA")); mBS2_w <- mBS2_[[1]]; mBS2_m <- mBS2_[[2]]
mBS3_ <- gee_gauss_sex(G2A_mega %>% filter(
  G2ASTUDY=="MHAS")); mBS3_w <- mBS3_[[1]]; mBS3_m <- mBS3_[[2]]
mBS4_ <- gee_gauss_sex(G2A_mega %>% filter(
  G2ASTUDY=="CRELES")); mBS4_w <- mBS4_[[1]]; mBS4_m <- mBS4_[[2]]
mBS5_ <- gee_gauss_sex(G2A_mega %>% filter(
  G2ASTUDY=="CHARLS")); mBS5_w <- mBS5_[[1]]; mBS5_m <- mBS5_[[2]]

#Women labels
b1w <- summary(mBS1_w)$coefficients[2,1]; ci1w <- confint(mBS1_w) %>% round(3)
b2w <- summary(mBS2_w)$coefficients[2,1]; ci2w <- confint(mBS2_w) %>% round(3)
b3w <- summary(mBS3_w)$coefficients[2,1]; ci3w <- confint(mBS3_w) %>% round(3)
b4w <- summary(mBS4_w)$coefficients[2,1]; ci4w <- confint(mBS4_w) %>% round(3)
b5w <- summary(mBS5_w)$coefficients[2,1]; ci5w <- confint(mBS5_w) %>% round(3)
ann_f1w <- c(
  paste0("β = ", b1w %>% round(3), "\n(", ci1w[2,1], "—", ci1w[2,2], ")"),
  paste0("β = ", b2w %>% round(3), "\n(", ci2w[2,1], "—", ci2w[2,2], ")"),
  paste0("β = ", b3w %>% round(3), "\n(", ci3w[2,1], "—", ci3w[2,2], ")"),
  paste0("β = ", b4w %>% round(3), "\n(", ci4w[2,1], "—", ci4w[2,2], ")"),
  paste0("β = ", b5w %>% round(3), "\n(", ci5w[2,1], "—", ci5w[2,2], ")"))

#Men labels
b1m <- summary(mBS1_m)$coefficients[2,1]; ci1m <- confint(mBS1_m) %>% round(3)
b2m <- summary(mBS2_m)$coefficients[2,1]; ci2m <- confint(mBS2_m) %>% round(3)
b3m <- summary(mBS3_m)$coefficients[2,1]; ci3m <- confint(mBS3_m) %>% round(3)
b4m <- summary(mBS4_m)$coefficients[2,1]; ci4m <- confint(mBS4_m) %>% round(3)
b5m <- summary(mBS5_m)$coefficients[2,1]; ci5m <- confint(mBS5_m) %>% round(3)
ann_f1m <- c(
  paste0("β = ", b1m %>% round(3), "\n(", ci1m[2,1], "—", ci1m[2,2], ")"),
  paste0("β = ", b2m %>% round(3), "\n(", ci2m[2,1], "—", ci2m[2,2], ")"),
  paste0("β = ", b3m %>% round(3), "\n(", ci3m[2,1], "—", ci3m[2,2], ")"),
  paste0("β = ", b4m %>% round(3), "\n(", ci4m[2,1], "—", ci4m[2,2], ")"),
  paste0("β = ", b5m %>% round(3), "\n(", ci5m[2,1], "—", ci5m[2,2], ")"))

#Function to arrange data
sdp2 <- function(x){
  fn_temp <- set.data.plot; y <- get(paste0("G2A_mega", x))
  rbind(
    fn_temp(filter(y, baseline==1, Sex=="Women")) %>%
      mutate(Sex="Women", ann=ann_f1w[x]),
    fn_temp(filter(y, baseline==1, Sex=="Men")) %>%
      mutate(Sex="Men", ann=ann_f1m[x]))}

#Trends by sex plot
st.order <- c("HRS","ELSA","MHAS","CRELES","CHARLS")
rbind(sdp2(1), sdp2(2), sdp2(3), sdp2(4), sdp2(5)) %>%
  mutate("St"=ordered(G2ASTUDY, levels=st.order), "nobs"=case_when(
    St=="HRS"~table(St)[1], St=="ELSA"~table(St)[2], St=="MHAS"~table(St)[3],
    St=="CRELES"~table(St)[4], St=="CHARLS"~table(St)[5]), "ycoord"=65,
    "Sex"=ordered(Sex, levels=c("Women","Men"))) %>%
  ggplot(aes(x=bi, y=mean, ymin=l.ci, ymax=u.ci, color=Par, group=Par)) +
  geom_line() + geom_pointrange() + geom_point() + theme_pubclean() +
  facet_grid(Sex~St) + scale_color_manual(values=c("#9683ec", "#35034f")) +
  scale_x_discrete(
    breaks=c("01/02","05/06","09/10","13/14","17/18"),
    labels=c("2001/02","2005/06","2009/10","2013/14","2017/18")) +
  labs(x=NULL, y="Weighted mean (years) with 95%CI", color=NULL, group=NULL)+
  theme(
    legend.position = "bottom", legend.key.size = unit(1,"cm"),
    legend.title = element_text(hjust=0.5, face = "bold", size=12),
    legend.text = element_text(hjust=0.5, size=11),
    strip.text = element_text(face="bold.italic", vjust=0.5, size=12),
    panel.spacing.x = unit(2.5,"lines"), panel.spacing.y = unit(1,"lines"),
    axis.text.x = element_text(vjust=0.7, hjust=0.8, size=12, angle=45),
    axis.text.y = element_text(size=12), axis.title=element_text(size=13)) +
  geom_text(aes(label = ifelse(nobs[nobs] > 2, ann, ""), y=ycoord),
            color="black", x="13/14", size=4.25) -> FigS4; FigS4

ggsave(FigS4, file="Figures/FigureS4.jpg", bg="transparent", #1500x750pix
       width=35, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

#Accelerated aging plot
data.p3 <- G2A_mega %>% filter(baseline==1, !is.na(Weights)) %>%
  transmute(JID, Ct=G2ASTUDY, accel, Sex, yr=year, bi, wt=std_weights) %>%
  filter((Ct=="HRS"&yr%in%(2e3+6*1:3))|(Ct!="HRS"&yr%in%(2e3+1:15))) %>%
  mutate(accel=accel*100); svyby(~accel, by=~bi+Ct+Sex, design=svydesign(
    data=data.p3, ids=~JID, weights=~wt), FUN=svymean) %>%
  `rownames<-`(NULL) %>% rename("mean"=accel) %>% mutate(
    "l.ci"=mean-qnorm(0.975)*se, "u.ci"=mean+qnorm(0.975)*se) %>%
  mutate("St"=ordered(Ct, levels=st.order),
         "Sex"=ordered(Sex, levels=c("Women", "Men"))) %>% 
  ggplot(aes(x=bi, y=mean)) + geom_errorbar(
    aes(ymin = l.ci, ymax = u.ci), width=0.35, linewidth=0.85) +
  geom_col(fill="#35034f") + facet_grid(Sex~St, scales="free_x") +
  theme_pubclean() + xlab(NULL) + ylab("Accelerated aging (%)") +
  scale_y_continuous(breaks = c(0,15,30,45,60)) +
  theme(
    legend.position = "bottom", legend.key.size = unit(1,"cm"),
    legend.title = element_text(hjust=0.5, face = "bold", size=12),
    legend.text = element_text(hjust=0.5, size=11),
    strip.text = element_text(face="bold.italic", vjust=0.5, size=12),
    panel.spacing.x = unit(2.5,"lines"), panel.spacing.y = unit(1,"lines"),
    axis.text.x = element_text(vjust=0.7, hjust=0.8, size=12, angle=45),
    axis.text.y = element_text(size=12), axis.title=element_text(size=13)) +
  geom_text(aes(label=round(mean,1), x=bi, y=mean-4.5),
            color="white", size=5, fontface="bold.italic") -> FigS5

ggsave(FigS5, file="Figures/FigureS5.jpg", bg="transparent",
       width=40, height=25, units=c("cm"), dpi=600, limitsize = FALSE)


#SFig6: By accelQ4----- ####
#Function to arrange data
sdp3 <- function(x){
  fn_temp <- set.data.plot; y <- get(paste0("G2A_mega", x))
  rbind(
    fn_temp(filter(y, baseline==1, accelQ4=="Q1")) %>%
      filter(Par!="Age") %>% mutate(accelQ4="Q1"),
    fn_temp(filter(y, baseline==1, accelQ4=="Q2")) %>%
      filter(Par!="Age") %>% mutate(accelQ4="Q2"),
    fn_temp(filter(y, baseline==1, accelQ4=="Q3")) %>%
      filter(Par!="Age") %>% mutate(accelQ4="Q3"),
    fn_temp(filter(y, baseline==1, accelQ4=="Q4")) %>%
      filter(Par!="Age") %>% mutate(accelQ4="Q4"),
    fn_temp(filter(y, baseline==1)) %>%
      filter(Par=="Age") %>% mutate(accelQ4="Age"))}

#Plot of trends with Age, Q1 and Q4
rbind(sdp3(1), sdp3(2), sdp3(3), sdp3(4), sdp3(5)) %>%
  mutate("St"=ordered(G2ASTUDY, levels=st.order)) %>%
  filter(accelQ4 %in% c("Age", "Q1", "Q4")) %>% mutate(accelQ4=factor(
    accelQ4, labels = c("Age", "AnthropoAge (acceleration Q1)",
                        "AnthropoAge (acceleration Q4)"))) %>%
  ggplot(aes(x=bi, y=mean,ymin=l.ci,ymax=u.ci, color=accelQ4,group=accelQ4)) +
  geom_line() + geom_pointrange() + geom_point() + theme_pubclean() +
  facet_wrap(~St, ncol=5, nrow=1) +
  scale_color_manual(values=c(scales::viridis_pal(option = "G")(4)[1:3])) +
  scale_x_discrete(breaks=c("01/02","05/06","09/10","13/14","17/18")) +
  labs(x=NULL, y="Weighted mean (years) with 95%CI", color=NULL, group=NULL)+
  theme(
    legend.position = "bottom", legend.key.size = unit(1,"cm"),
    legend.title = element_text(hjust=0.5, face = "bold", size=12),
    legend.text = element_text(hjust=0.5, size=11),
    strip.background.x = element_blank(),
    strip.text = element_text(face="bold.italic", vjust=0.5, size=12),
    panel.spacing.x = unit(2.5,"lines"), panel.spacing.y = unit(1,"lines"),
    axis.text.x = element_text(vjust=0.7, hjust=0.8, size=12, angle=45),
    axis.text.y = element_text(size=12), axis.title=element_text(size=13)
    ) -> FigS6

ggsave(FigS6, file="Figures/FigureS6.jpg", bg="transparent",
       width=35, height=15, units=c("cm"), dpi=600, limitsize = FALSE)

####-------------------------------- ACCELERATED AGING ------------#### ----####
#Kaplan Meier plot---------- ####
pooled0 <- pooled
km_by_study <- function(x,y=0.8,z=144,a="top"){
  i <- as.numeric(deparse(substitute(x)) %>% substr(7,7))+1
  studies <- c("Overall sample", "HRS", "ELSA", "MHAS", "CRELES", "CHARLS")
  km_mod <- survfit(Surv(TTDM, mortstat2) ~ (n_comorb>=2)+accel, data=x)
  kmf_st_<-ggsurvplot(
    km_mod, data=pooled, size=1, conf.int=T, risk.table=T,
    xlab="Follow-up (months)", ylab="Mortality (%)",
    legend.labs = c("Non-Accel+\nNo-Multimorb", "Accel+\nNo-Multimorb",
                    "Non-Accel+\nMultimorb","Accel+\nMultimorb"),
    surv.scale="percent", fun="cumhaz", pval = T,
    palette = c("#B8B8F7","#957ADF","#683798","#35034F"),
    risk.table.y.text.col=TRUE,  risk.table.y.text=FALSE, risk.table.title="",
    break.y.by= c(0.1), ylim=c(0,y),xlim=c(0,z),break.x.by= c(z/6),
    ggtheme = theme_light() +(theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size=15),
      text = element_text(hjust=0.5, family = "sans", size=15),
      legend.text = element_text(hjust=0.5, size=12),
      plot.subtitle = element_text(hjust=0.5, size=12)))
  ); kmf_st_$table <- kmf_st_$table + xlab("") +
    theme(panel.grid = element_blank(), title = element_blank()) +
    scale_x_continuous(labels = NULL, breaks = NULL, sec.axis = dup_axis(
      name = element_blank(), breaks = c(seq(0,z, by=z/6)),
      labels = c(seq(0,z, by=z/6)))
    ); kmf_st_$plot <- kmf_st_$plot + ggtitle(studies[i]) + theme(
      panel.grid.major.x = element_line(size = 0.2,linetype = 1),
      panel.grid.major.y = element_line(size = 0.2,linetype = 1),
      panel.grid.minor = element_blank(), legend.position = a,
    ); kmf_st<-kmf_st_$plot/kmf_st_$table +
    plot_layout(heights = c(5, 1.5))
  if (i==1){return (list(kmf_st, kmf_st_$plot))}
  else if (i!=1){return(kmf_st)}}

km_by_study(pooled0, 0.8, 144) -> fig4a #Overall


#SFig7: KM plot by study---- ####
km_by_study(pooled3, 0.8, 144, "none") -> fs7a #MHAS
km_by_study(pooled1, 0.8, 144, "none") -> fs7b #HRS
km_by_study(pooled5, 0.5, 96, "none") -> fs7c #CHARLS
km_by_study(pooled2, 0.5, 96, "none") -> fs7d #ELSA
km_by_study(pooled4, 0.5, 48, "none") -> fs7e #CRELES

FigS7 <- ggarrange(
  nrow=2, heights = c(1.2,1),
  ggarrange(
    fs7a, fs7b, nrow=1, labels = LETTERS[1:2],  common.legend = T,
    legend = "bottom", legend.grob = get_legend(fig4a[[2]])),
  ggarrange(
    fs7c, fs7d, fs7e, nrow=1, labels = LETTERS[3:5]))

ggsave(FigS7, file="Figures/FigureS7.jpg", bg="transparent", #1500x750pix
       width=35*1.1, height=20*1.1, units=c("cm"), dpi=600, limitsize = FALSE)


#Fig 4: Overall KM + Cox---- ####
cox.accel1 <- function(x){
  i <- as.numeric(deparse(substitute(x)) %>% substr(7,7))
  if(i>0){
    coxph(Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
            strata(Sex) + strata(coh) + strata(Ethnicity2),
          weights = std_weights, data = x)}
  else {
    coxph(Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
            strata(Sex) + strata(coh) + strata(Ethnicity2) +
            frailty(G2ASTUDY), weights = std_weights, data = x)}}
cox.accel2 <- function(x){
  i <- as.numeric(deparse(substitute(x)) %>% substr(7,7))
  if(i>0){
    coxph(Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
             Education + Smoking + Drinking +
            strata(Sex) + strata(coh) + strata(Ethnicity2),
          weights = std_weights, data = x)}
  else {
    coxph(Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
            Education + Smoking + Drinking +
            strata(Sex) + strata(coh) + strata(Ethnicity2) +
            frailty(G2ASTUDY), weights = std_weights, data = x)}}

#Simple adjustments
studies <- c("Overall", "HRS", "ELSA", "MHAS", "CRELES", "CHARLS")
mBS0 <- cox.accel1(pooled0); mBS1 <- cox.accel1(pooled1)
mBS2 <- cox.accel1(pooled2); mBS3 <- cox.accel1(pooled3)
mBS4 <- cox.accel1(pooled4); mBS5 <- cox.accel1(pooled5)
models <- list(mBS0, mBS1, mBS2, mBS3, mBS4, mBS5)
ann_f1 <- data.frame(
  "estimate"=c(1:6), "lower"=c(1:6), "upper"=c(1:6), "G2ASTUDY"=c(1:6),
  "Adj"=1, "events"=c(table(pooled0$mortstat2)[2], with(
    pooled0, tapply(mortstat2, G2ASTUDY, function(x){table(x)[2]})))
  ); for(i in 1:6){
    ann_f1$estimate[i] <- exp(coef(models[[i]])[1]) %>% round(3)
    ann_f1$lower[i] <- exp(confint(models[[i]])[1,1]) %>% round(2)
    ann_f1$upper[i] <- exp(confint(models[[i]])[1,2]) %>% round(2)
    ann_f1$G2ASTUDY[i] <- studies[i]}

#Full adjustments
mBS0 <- cox.accel2(pooled0); mBS1 <- cox.accel2(pooled1)
mBS2 <- cox.accel2(pooled2); mBS3 <- cox.accel2(pooled3)
mBS4 <- cox.accel2(pooled4); mBS5 <- cox.accel2(pooled5)
models <- list(mBS0, mBS1, mBS2, mBS3, mBS4, mBS5)
ann_f2 <- data.frame(
  "estimate"=c(1:6), "lower"=c(1:6), "upper"=c(1:6), "G2ASTUDY"=c(1:6),
  "Adj"=2, "events"=c(table(pooled0$mortstat2)[2], with(
    pooled0, tapply(mortstat2, G2ASTUDY, function(x){table(x)[2]})))
); for(i in 1:6){
  ann_f2$estimate[i] <- exp(coef(models[[i]])[1]) %>% round(3)
  ann_f2$lower[i] <- exp(confint(models[[i]])[1,1]) %>% round(2)
  ann_f2$upper[i] <- exp(confint(models[[i]])[1,2]) %>% round(2)
  ann_f2$G2ASTUDY[i] <- studies[i]}

rbind(ann_f1, ann_f2) %>%
  mutate("Adjustments"=ordered(Adj,1:2,c(
    "Age + multimorbidity", "+ education + lifestyle")),
    "G2ASTUDY"=ordered(G2ASTUDY, levels=studies),
    "lab1"=trunc(estimate*100)/100,
    "lab2"=ifelse(Adj==1,prettyNum(events,big.mark=","), NA)) %>%
  ggplot(aes(x=Adjustments, y=estimate, color=Adjustments)) +
  geom_pointrange(aes(size=events, ymin=lower, ymax=upper), shape=15) +
  scale_size_continuous(guide = "none", range = c(0.5, 1)) +
  facet_wrap(~G2ASTUDY, nrow=2) + scale_y_log10()+ 
  labs(x = "", y = "Adjusted HR (95%CI) for accelerated aging",
       title = "Accelerated aging and all-cause mortality") +
  geom_hline(yintercept = 1, lty = "dashed") + theme_pubclean() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#6C4CA9","#35034F")) +
  geom_text(aes(label=lab1), hjust=0, nudge_x=0.125, nudge_y=0.05,
            size=4, color=c(rep("#6C4CA9",6),rep("#35034F",6))) +
  geom_text(aes(label=lab2, y=1.15), nudge_x=0.5,
            size=3, color="black")+
  theme(plot.title = element_text(face = "bold", hjust=0.5),
        strip.text = element_text(face = "bold.italic", hjust=0.5),
        strip.background = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) -> fig4b; fig4b

fig4<-ggarrange(fig4a[[1]], fig4b, labels = LETTERS[1:2])
ggsave(fig4,file="Figures/Figure4.jpg", bg="transparent",
       width=35, height=15, units=c("cm"), dpi=600, limitsize = FALSE)

#Proportional hazards assumption
models <- paste0("mBS",1:5) %>% as.list %>% lapply(get)
#Accel and global p-values
models %>% sapply(function(x){(cox.zph(x))$table[c(1,7),3]}) %>%
  t %>% data.frame %>% `rownames<-`(studies[-1])
#Visualize
plotlist <- list(); for(i in 1:5){
  surveys = c("HRS", "ELSA", "MHAS", "CRELES", "CHARLS")
  pnew = survminer::ggcoxzph(
    (lapply(models, cox.zph))[[i]],
    ggtheme = ggpubr::theme_pubclean(), point.col = "midnightblue", 
    point.alpha = 0.5, var = "accel", se=T) + labs(subtitle=surveys[i])
  plotlist = append(plotlist, pnew)}; ggarrange(
    plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]],
    plotlist[[5]], nrow=2, ncol=3)
#Roughly constant hazards, variance increase in HRS


#SFig8: Cox by multimorb---- ####
cox.accel.comorb0 <- function(x){
  i <- as.numeric(deparse(substitute(x)) %>% substr(7,7))
  if(i>0){
    coxph(Surv(TTDM, mortstat2) ~ accel + Age +
            Education + Smoking + Drinking +
            strata(Sex) + strata(coh) + strata(Ethnicity2),
          weights = std_weights, data = x %>% filter(n_comorb<2))}
  else {
    coxph(Surv(TTDM, mortstat2) ~ accel + Age +
            Education + Smoking + Drinking +
            strata(Sex) + strata(coh) + strata(Ethnicity2) +
            frailty(G2ASTUDY), weights = std_weights,
          data = x %>% filter(n_comorb<2))}}
cox.accel.comorb1 <- function(x){
  i <- as.numeric(deparse(substitute(x)) %>% substr(7,7))
  if(i>0){
    coxph(Surv(TTDM, mortstat2) ~ accel + Age +
            Education + Smoking + Drinking +
            strata(Sex) + strata(coh) + strata(Ethnicity2),
          weights = std_weights, data = x %>% filter(n_comorb>=2))}
  else {
    coxph(Surv(TTDM, mortstat2) ~ accel + Age +
            Education + Smoking + Drinking +
            strata(Sex) + strata(coh) + strata(Ethnicity2) +
            frailty(G2ASTUDY), weights = std_weights,
          data = x %>% filter(n_comorb>=2))}}

studies <- c("Overall", "HRS", "ELSA", "MHAS", "CRELES", "CHARLS")
mAC0.0 <- cox.accel.comorb0(pooled0); mAC1.0 <- cox.accel.comorb1(pooled0)
mAC0.1 <- cox.accel.comorb0(pooled1); mAC1.1 <- cox.accel.comorb1(pooled1)
mAC0.2 <- cox.accel.comorb0(pooled2); mAC1.2 <- cox.accel.comorb1(pooled2)
mAC0.3 <- cox.accel.comorb0(pooled3); mAC1.3 <- cox.accel.comorb1(pooled3)
mAC0.4 <- cox.accel.comorb0(pooled4); mAC1.4 <- cox.accel.comorb1(pooled4)
mAC0.5 <- cox.accel.comorb0(pooled5); mAC1.5 <- cox.accel.comorb1(pooled5)

#No multimorbidity
models0 <- list(mAC0.0, mAC0.1, mAC0.2, mAC0.3, mAC0.4, mAC0.5)
pooled0.c0 <- filter(pooled0, n_comorb<2)
ann_c0 <- data.frame(
  "estimate"=c(1:6), "lower"=c(1:6), "upper"=c(1:6), "G2ASTUDY"=c(1:6),
  "Comorb"=0, "events"=c(table(pooled0.c0$mortstat2)[2], with(
    pooled0.c0, tapply(mortstat2, G2ASTUDY, function(x){table(x)[2]})))
); for(i in 1:6){
  ann_c0$estimate[i] <- exp(coef(models0[[i]])[1]) %>% round(3)
  ann_c0$lower[i] <- exp(confint(models0[[i]])[1,1]) %>% round(2)
  ann_c0$upper[i] <- exp(confint(models0[[i]])[1,2]) %>% round(2)
  ann_c0$G2ASTUDY[i] <- studies[i]}

#Multimorbidity
models1 <- list(mAC1.0, mAC1.1, mAC1.2, mAC1.3, mAC1.4, mAC1.5)
pooled0.c1 <- filter(pooled0, n_comorb>=2)
ann_c1 <- data.frame(
  "estimate"=c(1:6), "lower"=c(1:6), "upper"=c(1:6), "G2ASTUDY"=c(1:6),
  "Comorb"=1, "events"=c(table(pooled0.c1$mortstat2)[2], with(
    pooled0.c1, tapply(mortstat2, G2ASTUDY, function(x){table(x)[2]})))
); for(i in 1:6){
  ann_c1$estimate[i] <- exp(coef(models1[[i]])[1]) %>% round(3)
  ann_c1$lower[i] <- exp(confint(models1[[i]])[1,1]) %>% round(2)
  ann_c1$upper[i] <- exp(confint(models1[[i]])[1,2]) %>% round(2)
  ann_c1$G2ASTUDY[i] <- studies[i]}

rbind(ann_c0, ann_c1) %>%
  mutate("Multimorbidity"=ordered(Comorb, 0:1, c(
    "No multimorbidity", "Multimorbidity")),
    "G2ASTUDY"=ordered(G2ASTUDY, levels=studies),
    "lab1"=trunc(estimate*100)/100,
    "lab2"=prettyNum(events,big.mark=",")) %>%
  ggplot(aes(x=Multimorbidity, y=estimate, color=Multimorbidity)) +
  geom_pointrange(aes(size=events, ymin=lower, ymax=upper), shape=15) +
  scale_size_continuous(guide = "none", range = c(0.5, 1)) +
  facet_wrap(~G2ASTUDY, nrow=2) + labs(
    x = "", y = "Adjusted HR (95%CI) for accelerated aging",
    title = "Accelerated aging and all-cause mortality by multimorbidity") +
  geom_hline(yintercept = 1, lty = "dashed") + theme_pubclean() +
  theme(legend.position = "bottom") + scale_y_log10() +
  scale_color_manual(values = c("#6C4CA9", "#35034F")) +
  geom_text(aes(label=lab1), hjust=0, nudge_x=0.125, nudge_y=0.05,
            size=4, color=rep(c("#6C4CA9", "#35034F"),6)) +
  geom_text(aes(label=lab2), nudge_x=0.175, nudge_y=-0.05,
            size=3, color="black")+
  theme(plot.title = element_text(face = "bold", hjust=0.5),
        strip.text = element_text(face = "bold.italic", hjust=0.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank()) -> fs8; fs8

ggsave(fs8,file="Figures/FigureS8.jpg", bg="transparent",
       width=25, height=15, units=c("cm"), dpi=600, limitsize = FALSE)



#STab3: Cox cluster(ID)----- ####
#AnthropoAgeAccel Q1 vs Q2, Q3 and Q4
cox.Q <- function(x,y=0, adj=1){
  G2A_mega0 = G2A_mega
  x = deparse(substitute(G2A_mega)); z = get(paste0(x,y))
  if(adj==1){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accelQ4 + Age + comorb_cat +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data=z)
    return(output)
  }
  else if (adj==2){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accelQ4 + Age + comorb_cat +
        Education + Smoking + Drinking +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data=z)
    return(output)
  }
  else {
    return("Error: adjustment should be 1 or 2")}}
cQ0.1 <- cox.Q(G2A_mega, 0, 1); cQ0.2 <- cox.Q(G2A_mega, 0, 2)
cQ1.1 <- cox.Q(G2A_mega, 1, 1); cQ1.2 <- cox.Q(G2A_mega, 1, 2)
cQ2.1 <- cox.Q(G2A_mega, 2, 1); cQ2.2 <- cox.Q(G2A_mega, 2, 2)
cQ3.1 <- cox.Q(G2A_mega, 3, 1); cQ3.2 <- cox.Q(G2A_mega, 3, 2)
cQ4.1 <- cox.Q(G2A_mega, 4, 1); cQ4.2 <- cox.Q(G2A_mega, 4, 2)
cQ5.1 <- cox.Q(G2A_mega, 5, 1); cQ5.2 <- cox.Q(G2A_mega, 5, 2)

#AnthropoAgeAccel >0
cox.A <- function(x,y=0, adj=1){
  G2A_mega0 = G2A_mega
  x = deparse(substitute(G2A_mega)); z = get(paste0(x,y))
  if(adj==1){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data=z)
    return(output)
  }
  else if (adj==2){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
        Education + Smoking + Drinking +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data=z)
    return(output)
  }
  else {
    return("Error: adjustment should be 1 or 2")}}
cA0.1 <- cox.A(G2A_mega, 0, 1); cA0.2 <- cox.A(G2A_mega, 0, 2)
cA1.1 <- cox.A(G2A_mega, 1, 1); cA1.2 <- cox.A(G2A_mega, 1, 2)
cA2.1 <- cox.A(G2A_mega, 2, 1); cA2.2 <- cox.A(G2A_mega, 2, 2)
cA3.1 <- cox.A(G2A_mega, 3, 1); cA3.2 <- cox.A(G2A_mega, 3, 2)
cA4.1 <- cox.A(G2A_mega, 4, 1); cA4.2 <- cox.A(G2A_mega, 4, 2)
cA5.1 <- cox.A(G2A_mega, 5, 1); cA5.2 <- cox.A(G2A_mega, 5, 2)

#MHAS 2001
cox.MQ <- function(x, y, adj=1){
  if(adj==1){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accelQ4 + Age + comorb_cat +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data= x %>%
        filter(wave==y, !duplicated(ID)))
    return(output)
  }
  else if (adj==2){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accelQ4 + Age + comorb_cat +
        Education + Smoking + Drinking +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data= x %>%
        filter(wave==y, !duplicated(ID)))
    return(output)
  }
  else {
    return("Error: adjustment should be 1 or 2")}}
cox.MA <- function(x, y, adj=1){
  if(adj==1){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data= x %>%
        filter(wave==y, !duplicated(ID)))
    return(output)
  }
  else if (adj==2){
    output = coxph(
      Surv(TTDM, mortstat2) ~ accel + Age + comorb_cat +
        Education + Smoking + Drinking +
        strata(Sex) + strata(Ethnicity2) + strata(coh) +
        cluster(ID), weights=std_weights, data= x %>%
        filter(wave==y, !duplicated(ID)))
    return(output)
  }
  else {
    return("Error: adjustment should be 1 or 2")}}
cQ3a.1 <- cox.MQ(G2A_mega3, 1, 1); cQ3a.2 <- cox.MQ(G2A_mega3, 1, 2)
cA3a.1 <- cox.MA(G2A_mega3, 1, 1); cA3a.2 <- cox.MA(G2A_mega3, 1, 2)
#MHAS 2012
cQ3b.1 <- cox.MQ(G2A_mega3, 3, 1); cQ3b.2 <- cox.MQ(G2A_mega3, 3, 2)
cA3b.1 <- cox.MA(G2A_mega3, 3, 1); cA3b.2 <- cox.MA(G2A_mega3, 3, 2)

#Create table
getHR_ <- function(x){
  y = deparse(substitute(x)) %>% substr(2,2)
  if(y=="Q"){
    out = cbind(x$coefficients[1:3], confint(x)[1:3,]) %>%
      exp %>% `colnames<-`(LETTERS[1:3]) %>%
      as.data.frame() %>% mutate(A=sprintf(fmt="%#.2f", A),
                                 B=sprintf(fmt="%#.2f", B),
                                 C=sprintf(fmt="%#.2f", C)) %>%
      transmute("HR.CI"=paste0(A," (",B,"-",C,")"))
    return(out %>% `row.names<-`(NULL))}
  
  else if(y=="A"){
    out2 = c(x$coefficients[1], confint(x)[1,]) %>%
      exp %>% `names<-`(LETTERS[1:3]) %>% rbind() %>% 
      as.data.frame() %>% mutate(A=sprintf(fmt="%#.2f", A),
                                 B=sprintf(fmt="%#.2f", B),
                                 C=sprintf(fmt="%#.2f", C)) %>%
      transmute("HR.CI"=paste0(A," (",B,"-",C,")"))
    return(out2 %>% `row.names<-`(NULL))}}
studies1 <- c("Overall","HRS","ELSA"," Overall",
              "2001", "2012", "CRELES","CHARLS")

rbind(
  c("","","",""),
  t(rbind(getHR_(cQ0.1), getHR_(cA0.1))),
  t(rbind(getHR_(cQ1.1), getHR_(cA1.1))),
  t(rbind(getHR_(cQ2.1), getHR_(cA2.1))),
  t(rbind(getHR_(cQ3.1), getHR_(cA3.1))),
  t(rbind(getHR_(cQ3a.1), getHR_(cA3a.1))),
  t(rbind(getHR_(cQ3b.1), getHR_(cA3b.1))),
  t(rbind(getHR_(cQ4.1), getHR_(cA4.1))),
  t(rbind(getHR_(cQ5.1), getHR_(cA5.1))),
  c("","","",""),
  t(rbind(getHR_(cQ0.2), getHR_(cA0.2))),
  t(rbind(getHR_(cQ1.2), getHR_(cA1.2))),
  t(rbind(getHR_(cQ2.2), getHR_(cA2.2))),
  t(rbind(getHR_(cQ3.2), getHR_(cA3.2))),
  t(rbind(getHR_(cQ3a.2), getHR_(cA3a.2))),
  t(rbind(getHR_(cQ3b.2), getHR_(cA3b.2))),
  t(rbind(getHR_(cQ4.2), getHR_(cA4.2))),
  t(rbind(getHR_(cQ5.2), getHR_(cA5.2)))) %>% as.data.frame() %>%
  `colnames<-`(c("Q2 vs Q1","Q3 vs Q1","Q4 vs Q1","Accelerated\naging")) %>%
  `rownames<-`(c("", studies1, " ", paste0(studies1," "))
  ) -> tab_S3
rep("",3)

cbind(" "=c(
  "Adjustments 1: CA, comorbidities",studies[1:3],
  c(rep("MHAS",3),studies[1:2]),
  "Adjustments 2: + education, lifestyle",studies[1:3],
  c(rep("MHAS",3),studies[1:2])),
  tab_S3 %>% rownames_to_column("Study")) %>% flextable %>%
  merge_at(i=1, j=1:6) %>% merge_at(i=10, j=1:6) %>% merge_at(i=2, j=1:2) %>%
  merge_at(i=3, j=1:2) %>% merge_at(i=4, j=1:2) %>% merge_at(i=8, j=1:2) %>%
  merge_at(i=9, j=1:2) %>% merge_at(i=11, j=1:2) %>% merge_at(i=12, j=1:2)%>% 
  merge_at(i=13, j=1:2) %>% merge_at(i=17, j=1:2) %>%  merge_at(i=18, j=1:2)%>%
  bold(part = "header") %>% bold(i=c(1,10), j=1:6) %>% bold(i=1:18, j=1) %>% 
  italic(i=c(1,10), j=1:6) %>% italic(i=c(5:7, 14:16), j=2) %>%
  align(align="center", part="all") %>% autofit() -> tabs3

set_flextable_defaults(
  split=F, table_align="center", table.layout="autofit"); read_docx() %>%
  body_add_flextable(value=tabs3) %>% print(target = "Tables/TableS3.docx")


#Proportional hazards
models <- c(); for (i in c("Q", "A")){
  for (j in 1:2){models <- c(models, paste0("c", i, 0:5, ".", j))
  }}; as.list(models) %>%
  lapply(get) %>% sapply(function(x){
    modtest <- cox.zph(x); p <- modtest$table %>% nrow()
    modtest$table[c(1,p),c(3)]}) -> cox.zph.tab
#P-values for proportional hazards test
cox.zph.tab %>% t %>% `rownames<-`(models) %>% data.frame() %>%
  mutate("St"=rep(c(studies1[1:3],"MHAS",studies1[7:8]), 4),
         "Adj"=rep(c(rep("Simple",6), rep("Full",6)), 2),
         "Out"=c(rep("AccelQ4",12), rep("Accel",12))) %>% 
  arrange((GLOBAL)) %>% mutate("PR"=GLOBAL>(0.05))

#Visualize: AnthropoAgeAccel quartiles
plotlist <- list(); for(i in 1:6){
  surveys = c("Overall", "HRS", "ELSA", "MHAS", "CRELES", "CHARLS")
  pnew = survminer::ggcoxzph(
    (as.list(models[1:6]) %>% lapply(get) %>% lapply(cox.zph))[[i]],
    ggtheme = ggpubr::theme_pubclean(), point.col = "midnightblue", 
    point.alpha = 0.5, var = "accelQ4", se=F) + labs(subtitle=surveys[i])
  plotlist = append(plotlist, pnew)}
ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]],
          plotlist[[5]], plotlist[[6]], nrow=2, ncol=3)

#Visualize: AnthropoAgeAccel >=0
plotlist <- list(); for(i in 1:6){
  surveys = c("Overall", "HRS", "ELSA", "MHAS", "CRELES", "CHARLS")
  pnew = survminer::ggcoxzph(
    ((as.list(models[13:18]) %>% lapply(get) %>% lapply(cox.zph))[[i]]),
    ggtheme = ggpubr::theme_pubclean(), point.col = "midnightblue", 
    point.alpha = 0.5, var = "accel", se=F) + labs(subtitle=surveys[i])
  plotlist = append(plotlist, pnew)}
ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]],
          plotlist[[5]], plotlist[[6]], nrow=2, ncol=3)

####-------------------------------- FUNC/HEALTH DECLINE ----------#### ----####
#Correlation structures---------- ####
see.corrstr <- function(x){
  x %>% arrange(year) %>%
    pivot_wider(names_from=year, values_from=V, values_fill=NA) %>%
    select(-1) %>% cor(use="pairwise.complete.obs")}

### SRH ###
var1<-c(1,3,5,7); var2<-c(2,4,6)
(G2A_mega1 %>%transmute(ID,V=SRH,year) %>%see.corrstr)[var1,var1]
(G2A_mega1 %>%transmute(ID,V=SRH,year) %>%see.corrstr)[var2,var2]
G2A_mega2 %>%transmute(ID,V=SRH,year) %>% see.corrstr
G2A_mega3 %>%transmute(ID,V=SRH,year) %>% see.corrstr
G2A_mega4 %>%transmute(ID,V=SRH,year) %>% see.corrstr

### ADL ###
(G2A_mega1 %>%transmute(ID,V=ADL,year) %>%see.corrstr)[var1,var1]
(G2A_mega1 %>%transmute(ID,V=ADL,year) %>%see.corrstr)[var2,var2]
G2A_mega2 %>%transmute(ID,V=ADL,year) %>% see.corrstr
G2A_mega3 %>%transmute(ID,V=ADL,year) %>% see.corrstr
G2A_mega4 %>%transmute(ID,V=ADL,year) %>% see.corrstr
G2A_mega5 %>%transmute(ID,V=ADL,year) %>% see.corrstr

### IADL ###
(G2A_mega1 %>%transmute(ID,V=IADL,year) %>%see.corrstr)[var1,var1]
(G2A_mega1 %>%transmute(ID,V=IADL,year) %>%see.corrstr)[var2,var2]
G2A_mega2 %>%transmute(ID,V=IADL,year) %>% see.corrstr
G2A_mega3 %>%transmute(ID,V=IADL,year) %>% see.corrstr
G2A_mega4 %>%transmute(ID,V=IADL,year) %>% see.corrstr
G2A_mega5 %>%transmute(ID,V=IADL,year) %>% see.corrstr

### COMORBIDITIES ###
(G2A_mega1 %>%transmute(ID,V=n_comorb,year) %>%see.corrstr)[var1,var1]
(G2A_mega1 %>%transmute(ID,V=n_comorb,year) %>%see.corrstr)[var2,var2]
G2A_mega2 %>%transmute(ID,V=n_comorb,year) %>% see.corrstr
G2A_mega3 %>%transmute(ID,V=n_comorb,year) %>% see.corrstr
G2A_mega4 %>%transmute(ID,V=n_comorb,year) %>% see.corrstr
G2A_mega5 %>%transmute(ID,V=n_comorb,year) %>% see.corrstr


#Poisson GEE: SRH---------------- ####
## Self-reported health ##
###________________________________ SRH ________________________________###
#New data frame
G2A_SRH <- G2A_mega %>% mutate("Out"=as.numeric(SRH)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.SRH0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_SRH)
gee.SRH1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_SRH)
gee.SRH2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_SRH)
gee.SRH3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_SRH)
#QIC comparison
QIC.SRH <- QIC(gee.SRH1)-QIC(gee.SRH0)#SRH
#At risk and incident cases
risk.SRH <- nrow(G2A_SRH %>% filter(!duplicated(ID)))


#Poisson GEE: ADL, IADL---------- ####
###________________________________ ADL ________________________________###
#New data frame
ids_ADL <- (pooled %>% filter(ADL>=1))$ID
G2A_ADL <- G2A_mega %>% filter(!(ID %in% ids_ADL)) %>%
  mutate("Out"=as.numeric(ADL)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.ADL0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_ADL)
gee.ADL1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_ADL)
gee.ADL2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_ADL)
gee.ADL3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_ADL)
#QIC comparison
QIC.ADL <- QIC(gee.ADL1)-QIC(gee.ADL0)#ADL
#At risk and incident cases
risk.ADL <- table(G2A_ADL$ID) %>% length
case.ADL <- sum(table(G2A_ADL$ID, G2A_ADL$Out>0)[,2]>0)


###________________________________ IADL ________________________________###
#New data frame
ids_IAD <- (pooled %>% filter(IADL>=1))$ID
G2A_IAD <- G2A_mega %>% filter(!(ID %in% ids_IAD)) %>%
  mutate("Out"=as.numeric(IADL)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.IAD0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_IAD)
gee.IAD1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_IAD)
gee.IAD2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_IAD)
gee.IAD3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_IAD)
#QIC comparison
QIC.IAD <- QIC(gee.IAD1)-QIC(gee.IAD0)#IAD
#At risk and incident cases
risk.IAD <- table(G2A_IAD$ID) %>% length
case.IAD <- sum(table(G2A_IAD$ID, G2A_IAD$Out>0)[,2]>0)


#Poisson GEE: Diseases----------- ####
###________________________________ T2D ________________________________###
#New data frame
ids_T2D <- (pooled %>% filter(T2D==1))$ID
G2A_T2D <- G2A_mega %>% filter(!(ID %in% ids_T2D)) %>%
  mutate("Out"=as.numeric(T2D)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.T2D0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_T2D)
gee.T2D1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_T2D)
gee.T2D2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_T2D)
gee.T2D3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_T2D)
#QIC comparison
QIC.T2D <- QIC(gee.T2D1)-QIC(gee.T2D0)#T2D
#At risk and incident cases
risk.T2D <- table(G2A_T2D$ID) %>% length
case.T2D <- sum(table(G2A_T2D$ID, G2A_T2D$Out>0)[,2]>0)


###________________________________ HBP ________________________________###
#New data frame
ids_HBP <- (pooled %>% filter(HBP>=1))$ID
G2A_HBP <- G2A_mega %>% filter(!(ID %in% ids_HBP)) %>%
  mutate("Out"=as.numeric(HBP)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.HBP0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_HBP)
gee.HBP1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_HBP)
gee.HBP2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_HBP)
gee.HBP3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_HBP)
#QIC comparison
QIC.HBP <- QIC(gee.HBP1)-QIC(gee.HBP0)#HBP
#At risk and incident cases
risk.HBP <- table(G2A_HBP$ID) %>% length
case.HBP <- sum(table(G2A_HBP$ID, G2A_HBP$Out>0)[,2]>0)


###________________________________ CAN ________________________________###
#New data frame
ids_CAN <- (pooled %>% filter(CAN>=1))$ID
G2A_CAN <- G2A_mega %>% filter(!(ID %in% ids_CAN)) %>%
  mutate("Out"=as.numeric(CAN)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.CAN0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_CAN)
gee.CAN1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_CAN)
gee.CAN2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_CAN)
gee.CAN3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_CAN)
#QIC comparison
QIC.CAN <- QIC(gee.CAN1)-QIC(gee.CAN0)#CAN
#At risk and incident cases
risk.CAN <- table(G2A_CAN$ID) %>% length
case.CAN <- sum(table(G2A_CAN$ID, G2A_CAN$Out>0)[,2]>0)


###________________________________ LUD ________________________________###
#New data frame
ids_LUD <- (pooled %>% filter(LUD>=1))$ID
G2A_LUD <- G2A_mega %>% filter(!(ID %in% ids_LUD)) %>%
  mutate("Out"=as.numeric(LUD)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.LUD0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_LUD)
gee.LUD1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_LUD)
gee.LUD2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_LUD)
gee.LUD3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_LUD)
#QIC comparison
QIC.LUD <- QIC(gee.LUD1)-QIC(gee.LUD0)#LUD
#At risk and incident cases
risk.LUD <- table(G2A_LUD$ID) %>% length
case.LUD <- sum(table(G2A_LUD$ID, G2A_LUD$Out>0)[,2]>0)


###________________________________ AMI ________________________________###
#New data frame
ids_AMI <- (pooled %>% filter(AMI>=1))$ID
G2A_AMI <- G2A_mega %>% filter(!(ID %in% ids_AMI)) %>%
  mutate("Out"=as.numeric(AMI)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.AMI0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_AMI)
gee.AMI1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_AMI)
gee.AMI2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_AMI)
gee.AMI3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_AMI)
#QIC comparison
QIC.AMI <- QIC(gee.AMI1)-QIC(gee.AMI0)#AMI
#At risk and incident cases
risk.AMI <- table(G2A_AMI$ID) %>% length
case.AMI <- sum(table(G2A_AMI$ID, G2A_AMI$Out>0)[,2]>0)


###________________________________ EVC ________________________________###
#New data frame
ids_EVC <- (pooled %>% filter(EVC>=1))$ID
G2A_EVC <- G2A_mega %>% filter(!(ID %in% ids_EVC)) %>%
  mutate("Out"=as.numeric(EVC)) %>% filter(Out>=0) %>%
  select(ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex,
         Ethnicity3, coh, Education, Smoking, Drinking) %>% na.omit

#GEE models with Poisson variance function
gee.EVC0 <- geeglm( #Without AnthropoAgeAccel
  Out ~ Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_EVC)
gee.EVC1 <- geeglm( #AnthropoAgeAccel
  Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_EVC)
gee.EVC2 <- geeglm( #AnthropoAgeAccel>=0
  Out ~ (AnthropoAgeAccel>0) + Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_EVC)
gee.EVC3 <- geeglm( #AnthropoAgeAccel*Age
  Out ~ scale(AnthropoAgeAccel)*Age + Sex + Ethnicity3 +
    coh + Education + Smoking + Drinking, offset=log(TTDY+1),
  id=ID, wave=year, std.err="san.se", corstr = "ar1",
  family="poisson", weights=std_weights, data=G2A_EVC)
#QIC comparison
QIC.EVC <- QIC(gee.EVC1)-QIC(gee.EVC0)#EVC
#At risk and incident cases
risk.EVC <- table(G2A_EVC$ID) %>% length
case.EVC <- sum(table(G2A_EVC$ID, G2A_EVC$Out>0)[,2]>0)



#Tab 1: Functional/health decline ####
getRR_gee <- function(x, y){
  if(y%in%1:2){
    a <- summary(x)$coefficients[2,][1] %>% as.numeric
    b <- summary(x)$coefficients[2,][2] %>% as.numeric}
  else if(y==3){
    z <- nrow(summary(x)$coefficients)
    a <- summary(x)$coefficients[z,][1] %>% as.numeric
    b <- summary(x)$coefficients[z,][2] %>% as.numeric}
  c <- a %>% exp %>% sprintf(fmt="%#.3f")
  d <- (a + c(-1,1)*qnorm(.975)*b) %>% exp %>% sprintf(fmt="%#.3f")
  e <- paste0(c, " (", d[1], "-", d[2],")")
  e}

outs <- c("SRH", "ADL", "IAD", "T2D", "HBP", "CAN", "LUD", "AMI", "EVC")
at_risk <- as.list(paste0("risk.", outs)) %>% sapply(get)
events <- c("-", as.list(paste0("case.", outs[-1])) %>% sapply(get))
dQIC <- as.list(paste0("QIC.", outs)) %>% sapply(function(x){get(x)[1]})
RR1 <- as.list(paste0("gee.",outs,1)) %>% lapply(get) %>% sapply(getRR_gee,1)
RR2 <- as.list(paste0("gee.",outs,2)) %>% lapply(get) %>% sapply(getRR_gee,2)
RR3 <- as.list(paste0("gee.",outs,3)) %>% lapply(get) %>% sapply(getRR_gee,3)
names <- c("Change in self-\nreported health", paste0("New onset\n", c(
  "ADL deficit","IADL deficit","diabetes","hypertension","cancer diagnosis",
  "chronic lung disease", "myocardial infarction","stroke")))

data.frame(names, at_risk, events, RR1, RR2, RR3, dQIC) %>% `colnames<-`(c(
  "Outcome", "Population at\nrisk (baseline)", "Events",
  "RR.AnthropoAgeAccel (1-SD)", "RR.AnthropoAgeAccel ≥0",
  "RR.AnthropoAgeAccel*Age", "Delta QIC")) %>% flextable() %>% 
  separate_header(opts = c("span-top", "bottom-vspan", "center-hspan")) %>%
  bold(part="header") %>% italic(j=1, part="body") %>%
  align(align="center", part="all") -> tab4

set_flextable_defaults(
  split=F, table_align="center", table.layout="autofit"); read_docx() %>%
  body_add_flextable(value=tab4) %>% print(target = "Tables/Table1.docx")



#SFig9: Sensitivity: non incident ####
gee_get0 <- function(x,y){
  outs_gee = c("SRH", "ADL", "IADL", "n_comorb")
  outcome = (x %>% select(outs_gee[y]))[[1]]
  x2 = x %>% mutate("Out"=outcome) %>% select(
    ID, year, std_weights, Out, TTDY, AnthropoAgeAccel,
    Age, Sex, Ethnicity3) %>% na.omit
  geeglm(formula = Out ~ scale(AnthropoAgeAccel) + Age + Ethnicity3 + Sex,
         offset=log(TTDY+1), id=ID, wave=year, std.err="san.se",
         corstr = "ar1", family="poisson", weights=std_weights, data=x2)}
gee_get1 <- function(x,y){
  outs_gee = c("SRH", "ADL", "IADL", "n_comorb")
  outcome = (x %>% select(outs_gee[y]))[[1]]
  x2 = x %>% mutate("Out"=outcome) %>% select(
    ID, year, std_weights, Out, TTDY, AnthropoAgeAccel,
    Age, Sex, Ethnicity) %>% na.omit
  geeglm(formula = Out ~ scale(AnthropoAgeAccel) + Age + Sex + Ethnicity,
         offset=log(TTDY+1), id=ID, wave=year, std.err="san.se",
         corstr = "ar1", family="poisson", weights=std_weights, data=x2)}
gee_get2 <- function(x,y){
  outs_gee = c("SRH", "ADL", "IADL", "n_comorb")
  outcome = (x %>% select(outs_gee[y]))[[1]]
  x2 = x %>% mutate("Out"=outcome) %>% select(
    ID, year, std_weights, Out, TTDY, AnthropoAgeAccel, Age, Sex) %>% na.omit
  geeglm(formula = Out ~ scale(AnthropoAgeAccel) + Age + Sex,
         offset=log(TTDY+1), id=ID, wave=year, std.err="san.se",
         corstr = "ar1", family="poisson", weights=std_weights, data=x2)}
getRR_gee2 <- function(x, y=1){
  if(y%in%1:2){
    a <- summary(x)$coefficients[2,][1] %>% as.numeric
    b <- summary(x)$coefficients[2,][2] %>% as.numeric}
  else if(y==3){
    z <- nrow(summary(x)$coefficients)
    a <- summary(x)$coefficients[z,][1] %>% as.numeric
    b <- summary(x)$coefficients[z,][2] %>% as.numeric}
  c <- a %>% exp
  d <- (a + c(-1,1)*qnorm(.975)*b) %>% exp
  e <- c(c, d[1], d[2])
  e}

G2A_GEE <- G2A_mega %>% mutate(
  "SRH"=as.numeric(SRH), "ADL"=as.numeric(ADL),
  "IADL"=as.numeric(IADL), "n_comorb"=as.numeric(n_comorb)) %>%
  filter(SRH>=0, ADL>=0, IADL>=0, n_comorb>=0); G2A_GEE0 <- G2A_GEE
G2A_GEE1 <- G2A_GEE %>% filter(G2ASTUDY=="HRS")
G2A_GEE2 <- G2A_GEE %>% filter(G2ASTUDY=="ELSA")
G2A_GEE3 <- G2A_GEE %>% filter(G2ASTUDY=="MHAS")
G2A_GEE4 <- G2A_GEE %>% filter(G2ASTUDY=="CRELES")
G2A_GEE5 <- G2A_GEE %>% filter(G2ASTUDY=="CHARLS")

#Overall
GEE0.1 <- gee_get0(G2A_GEE0, 1); GEE0.2 <- gee_get0(G2A_GEE0, 2)
GEE0.3 <- gee_get0(G2A_GEE0, 3); GEE0.4 <- gee_get0(G2A_GEE0, 4)
#HRS
GEE1.1 <- gee_get1(G2A_GEE1, 1); GEE1.2 <- gee_get1(G2A_GEE1, 2)
GEE1.3 <- gee_get1(G2A_GEE1, 3); GEE1.4 <- gee_get1(G2A_GEE1, 4)
#ELSA
GEE2.1 <- gee_get2(G2A_GEE2, 1); GEE2.2 <- gee_get2(G2A_GEE2, 2)
GEE2.3 <- gee_get2(G2A_GEE2, 3); GEE2.4 <- gee_get2(G2A_GEE2, 4)
#MHAS
GEE3.1 <- gee_get2(G2A_GEE3, 1); GEE3.2 <- gee_get2(G2A_GEE3, 2)
GEE3.3 <- gee_get2(G2A_GEE3, 3); GEE3.4 <- gee_get2(G2A_GEE3, 4)
#CRELES
GEE4.1 <- gee_get2(G2A_GEE4, 1); GEE4.2 <- gee_get2(G2A_GEE4, 2)
GEE4.3 <- gee_get2(G2A_GEE4, 3); GEE4.4 <- gee_get2(G2A_GEE4, 4)
#CHARLES
GEE5.1 <- gee_get2(G2A_GEE5, 1); GEE5.2 <- gee_get2(G2A_GEE5, 2)
GEE5.3 <- gee_get2(G2A_GEE5, 3); GEE5.4 <- gee_get2(G2A_GEE5, 4)

EST0.1 <- getRR_gee2(GEE0.1); EST0.2 <- getRR_gee2(GEE0.2)
EST0.3 <- getRR_gee2(GEE0.3); EST0.4 <- getRR_gee2(GEE0.4)

EST1.1 <- getRR_gee2(GEE1.1); EST1.2 <- getRR_gee2(GEE1.2)
EST1.3 <- getRR_gee2(GEE1.3); EST1.4 <- getRR_gee2(GEE1.4)

EST2.1 <- getRR_gee2(GEE2.1); EST2.2 <- getRR_gee2(GEE2.2)
EST2.3 <- getRR_gee2(GEE2.3); EST2.4 <- getRR_gee2(GEE2.4)

EST3.1 <- getRR_gee2(GEE3.1); EST3.2 <- getRR_gee2(GEE3.2)
EST3.3 <- getRR_gee2(GEE3.3); EST3.4 <- getRR_gee2(GEE3.4)

EST4.1 <- getRR_gee2(GEE4.1); EST4.2 <- getRR_gee2(GEE4.2)
EST4.3 <- getRR_gee2(GEE4.3); EST4.4 <- getRR_gee2(GEE4.4)

EST5.1 <- getRR_gee2(GEE5.1); EST5.2 <- getRR_gee2(GEE5.2)
EST5.3 <- getRR_gee2(GEE5.3); EST5.4 <- getRR_gee2(GEE5.4)

studies <- c(
  "Overall", "HRS", "ELSA", "MHAS", "CRELES", "CHARLS"); outcomes <- c(
    "Changes in self-\nreported health", "IADL deficits", "ADL deficits",
    "Number of \ncomorbidities"); rbind(
  as.list(paste0("EST", 0:5, ".1")) %>% sapply(get) %>% t,
  as.list(paste0("EST", 0:5, ".2")) %>% sapply(get) %>% t,
  as.list(paste0("EST", 0:5, ".3")) %>% sapply(get) %>% t,
  as.list(paste0("EST", 0:5, ".4")) %>% sapply(get) %>% t) %>% 
  data.frame %>% `colnames<-`(c("mean", "lower", "upper")) %>% mutate(
    "Study"=factor(rep(0:5,4), levels=0:5, labels=studies),
    "Out"=factor(c(rep(1,6),rep(2,6),rep(3,6),rep(4,6)), 1:4, outcomes)) %>%
  
  ggplot(aes(x=mean, y=Out, color=Out, shape=Out, fill=Out)) +
  geom_vline(xintercept=1, linetype=2) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  labs(y=NULL, x="\nAdjusted rate ratio for AnthropoAgeAccel (95%CI)", 
       color="Outcome", shape="Outcome", fill="Outcome") +
  facet_wrap(~Study, nrow=3, ncol=2) + theme_pubclean() +
  scale_shape_manual(values=c(19,24,23,15),
                     guide=guide_legend(reverse = T,  byrow=T, ncol=2)) +
  scale_color_manual(values=c("#9683ec","#66439E","#8369CD","#35034F"),
                     guide=guide_legend(reverse = T,  byrow=T, ncol=2)) +
  scale_fill_manual(values=c("#9683ec","#66439E","#8369CD","#35034F"),
                    guide=guide_legend(reverse = T,  byrow=T, ncol=2)) +
  theme(
    legend.position = "bottom",  legend.key.size = unit(1,"cm"),
    legend.title = element_text(hjust=0.5, face = "bold", size=12),
    legend.text = element_text(hjust=0.5, size=11),
    strip.text = element_text(face = "bold.italic", size=11),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(vjust=1.5), strip.background = element_blank(),
    axis.title.x = element_text(face = "italic"),
    panel.spacing.y = unit(2,"lines"),
    panel.spacing.x = unit(4,"lines")) -> FigS9
  
ggsave(FigS9, file="Figures/FigureS9.jpg", bg="transparent",
       width=15, height=15, units=c("cm"), dpi=600, limitsize = FALSE)


####----DRAFTS-----------------------------------------------------#### ----####
#Sankey plots >0-------- ####
## CHARLS ##
sp.H <- G2A_mega %>% filter(G2ASTUDY=="CHARLS", participation==7, year%in%c(2006,2010,2018)) %>% select(ID,accel,year) %>%
  arrange(-desc(year)) %>% pivot_wider(names_from = year, values_from = accel) %>% transmute(
    `2006`=factor(`2006`, 0:1, c("Non-accelerated","Accelerated")),
    `2010`=factor(`2010`, 0:1, c("Non-accelerated","Accelerated")),
    `2018`=factor(`2018`, 0:1, c("Non-accelerated","Accelerated"))) %>% na.omit %>% make_long(`2006`,`2010`,`2018`); sp.H %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(flow.alpha = .5, node.color = "gray30", width=0.15) + labs(title="CHARLS", x=NULL) +
  geom_sankey_text(size = 4.5, angle=90, fontface="bold.italic", color="white") +
  scale_fill_manual(values=c("#35034f", "#9683ec")) + theme_sankey(base_size = 17) +
  theme(legend.position = "none", axis.text.x = element_text(hjust=0.5, vjust=5, face="bold.italic"),
        plot.title = element_text(size=17, hjust=0.5, vjust=-2, face="bold.italic")) -> Fsankey_CHARLS.1

## MHAS ##
sp.M <- G2A_mega %>% filter(G2ASTUDY=="MHAS", participation==3) %>% select(ID,accel,year) %>% arrange(-desc(year)) %>%
  pivot_wider(names_from = year, values_from = accel) %>% transmute(
    `2001`=factor(`2001`, 0:1, c("Non-accelerated","Accelerated")),
    `2003`=factor(`2003`, 0:1, c("Non-accelerated","Accelerated")),
    `2012`=factor(`2012`, 0:1, c("Non-accelerated","Accelerated"))) %>% na.omit %>% make_long(`2001`,`2003`,`2012`); sp.M %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(flow.alpha = .5, node.color = "gray30", width=0.15) + labs(title="MHAS", x=NULL) +
  geom_sankey_text(size = 4.5, angle=90, fontface="bold.italic", color="white") +
  scale_fill_manual(values=c("#35034f", "#9683ec")) + theme_sankey(base_size = 17) +
  theme(legend.position = "none", axis.text.x = element_text(hjust=0.5, vjust=5, face="bold.italic"),
        plot.title = element_text(size=17, hjust=0.5, vjust=-2, face="bold.italic")) -> Fsankey_MHAS.1

## ELSA ##
sp.E <- G2A_mega %>% filter(G2ASTUDY=="ELSA", participation==3) %>% select(ID,accel,year) %>% arrange(-desc(year)) %>%
  pivot_wider(names_from = year, values_from = accel) %>% transmute(
    `2004`=factor(`2004`, 0:1, c("Non-accelerated","Accelerated")),
    `2008`=factor(`2008`, 0:1, c("Non-accelerated","Accelerated")),
    `2012`=factor(`2012`, 0:1, c("Non-accelerated","Accelerated"))) %>% na.omit %>% make_long(`2004`,`2008`,`2012`); sp.E %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(flow.alpha = .5, node.color = "gray30", width=0.15) + labs(title="ELSA", x=NULL) +
  geom_sankey_text(size = 4.5, angle=90, fontface="bold.italic", color="white") +
  scale_fill_manual(values=c("#35034f", "#9683ec")) + theme_sankey(base_size = 17) +
  theme(legend.position = "none", axis.text.x = element_text(hjust=0.5, vjust=5, face="bold.italic"),
        plot.title = element_text(size=17, hjust=0.5, vjust=-2, face="bold.italic")) -> Fsankey_ELSA.1

## CHARLS ##
sp.C <- G2A_mega %>% filter(G2ASTUDY=="CHARLS", participation==3) %>% select(ID,accel,year) %>% arrange(-desc(year)) %>%
  pivot_wider(names_from = year, values_from = accel) %>% transmute(
    `2011`=factor(`2011`, 0:1, c("Non-accelerated","Accelerated")),
    `2013`=factor(`2013`, 0:1, c("Non-accelerated","Accelerated")),
    `2015`=factor(`2015`, 0:1, c("Non-accelerated","Accelerated"))) %>% na.omit %>% make_long(`2011`,`2013`,`2015`); sp.C %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(flow.alpha = .5, node.color = "gray30", width=0.15) + labs(title="CHARLS", x=NULL) +
  geom_sankey_text(size = 4.5, angle=90, fontface="bold.italic", color="white") +
  scale_fill_manual(values=c("#35034f", "#9683ec")) + theme_sankey(base_size = 17) +
  theme(legend.position = "none", axis.text.x = element_text(hjust=0.5, vjust=5, face="bold.italic"),
        plot.title = element_text(size=17, hjust=0.5, vjust=-2, face="bold.italic")) -> Fsankey_CHARLS.1

## Combine + labels ##
ggarrange(Fsankey_CHARLS.1,Fsankey_MHAS.1,Fsankey_ELSA.1,Fsankey_CHARLS.1, nrow=1, ncol=4) -> Fsankey1
sp.H %>% group_by(x) %>% summarise(Acc=(table(node) %>% prop.table)[1]*100, NAcc=(table(node) %>% prop.table)[2]*100)
sp.M %>% group_by(x) %>% summarise(Acc=(table(node) %>% prop.table)[1]*100, NAcc=(table(node) %>% prop.table)[2]*100)
sp.E %>% group_by(x) %>% summarise(Acc=(table(node) %>% prop.table)[1]*100, NAcc=(table(node) %>% prop.table)[2]*100)
sp.C %>% group_by(x) %>% summarise(Acc=(table(node) %>% prop.table)[1]*100, NAcc=(table(node) %>% prop.table)[2]*100)
#Regression models------ ####
health_var <- c("SRH", "ADL", "IADL", "HBP", "T2D", "CAN", "LUD", "AMI", "EVC", "ART")
G2A_mega$n_comorb <- G2A_mega %>% select(all_of(health_var[4:10])) %>% apply(1,sum, na.rm=T)
G2A_mega <- G2A_mega %>% filter(n_comorb>=0, ADL>=0, IADL>=0) 

M.AR1 <- geeglm(n_comorb~(AnthropoAgeAccel)*Country2+Sex+Age,offset=log(years+1),id=ID , wave = wave, 
                 std.err="san.se",   data=G2A_mega, weights=Weights, family = "poisson");summary(M.AR1)
M.AR2 <- geeglm(ADL~(AnthropoAgeAccel)*Country2+Sex+Age,offset=log(years+1),id=ID , wave = wave, 
                 std.err="san.se",   data=G2A_mega, weights=Weights, family = "poisson"); summary(M.AR2)
M.AR3 <- geeglm(IADL~(AnthropoAgeAccel)*Country2+Sex+Age,offset=log(years+1),id=ID , wave = wave, 
                 std.err="san.se",   data=G2A_mega, weights=Weights, family = "poisson"); summary(M.AR3)
CI.AR1 <- confint(M.AR1); CI.AR2 <- confint(M.AR2); CI.AR3 <- confint(M.AR3)

G.AR1 <- geeglm(n_comorb~(AnthropoAgeAccel)+Sex+Age,offset=log(years+1),id=ID , wave = wave, 
                 std.err="san.se",   data=G2A_mega, weights=Weights, family = "poisson")
G.AR2 <- geeglm(ADL~(AnthropoAgeAccel)+Sex+Age,offset=log(years+1),id=ID , wave = wave, 
                 std.err="san.se",   data=G2A_mega, weights=Weights, family = "poisson")
G.AR3 <- geeglm(IADL~(AnthropoAgeAccel)+Sex+Age,offset=log(years+1),id=ID , wave = wave, 
                 std.err="san.se",   data=G2A_mega, weights=Weights, family = "poisson")

FB.AR1 <- coefficients(summary(M.AR1))[c(2,9:12),1]; FC1 <- CI.AR1[c(2,9:12),]
FB.AR2 <- coefficients(summary(M.AR2))[c(2,9:12),1]; FC2 <- CI.AR2[c(2,9:12),]
FB.AR3 <- coefficients(summary(M.AR3))[c(2,9:12),1]; FC3 <- CI.AR3[c(2,9:12),]

GB.AR1 <- coefficients(summary(G.AR1))[c(2),1]; GC1 <- confint(G.AR1)[2,]
GB.AR2 <- coefficients(summary(G.AR2))[c(2),1]; GC2 <- confint(G.AR2)[2,]
GB.AR3 <- coefficients(summary(G.AR3))[c(2),1]; GC3 <- confint(G.AR3)[2,]


#Years and deaths------- ####

#Year of interview
#HRS
(HRS_alt %>% select(all_of(paste0("R",8:15,"IWENDY"))) %>% as.list %>%
    lapply(table)) %>% lapply(function(x){names(x)[1:2]}) %>%
  data.frame %>% t %>% data.frame
#MHAS
(MHAS %>% select(all_of(paste0("r",1:5,"iwy"))) %>% as.list %>%
    lapply(table)) %>% lapply(names) %>% data.frame %>%
  t %>% data.frame
#CHARLS
(CHARLS %>% select(all_of(paste0("r",1:4,"iwy"))) %>% as.list %>%
    lapply(table)) %>% lapply(names) %>% data.frame %>%
  t %>% data.frame
#ELSA
(ELSA %>% select(all_of(paste0("r",2:9,"iwindy"))) %>% as.list %>%
    lapply(table)) %>% lapply(names) %>% data.frame %>%
  t %>% data.frame
#CRELES
(CRELES %>% select(all_of(paste0("r",1:3,"iwy"))) %>% as.list %>%
    lapply(table)) %>% lapply(function(x){names(x)[2]}) %>%
  data.frame %>% t %>% data.frame

#Year of dearh
HRS_alt$RADYEAR %>% table #2020
MHAS$radyear %>% table #2018
CHARLS$radyear %>% table #2013
ELSA$radyear %>% table #2012
CRELES$radyear %>% table #2009
