
# Source debuged resid function for dlm_resid
source('dlm_resid_fixed.R')

# Function to find ages to forecast given the data supplied. ####
# There are no siblings to forecast jacks, and jacks don't affect mgmt decisions 
# Generally counts of the oldest age class are small, they contribute little to total forecast

find_ages_to_forecast <- function(df){
age_cols <- colnames(df)[str_detect(colnames(df),"Age")]
# age_cols[-c(1,length(age_cols))] %>% parse_number()
age_cols[-c(1)] %>% parse_number()
}

# Convert a Return Table to a Brood Table ####
return_to_brood <- function(rt){
  
# Find the min age column in the return table (rt), then make "sym" for tidy eval below
 min_age <- colnames(rt)[str_detect(colnames(rt),"Age")] %>% min() %>% rlang::sym()

# Reshape the data, filter rows with complete min_age data
  rt %>% 
  group_by(Stock) %>% 
  gather(key="AgeName", value="Return", contains("Age")) %>% 
  mutate(Age=parse_number(AgeName)) %>% 
  mutate(BroodYear=ReturnYear - Age) %>% 
  select(-ReturnYear, -Age) %>% 
  spread(AgeName, Return) %>%
    # Unquo the min_age thing
    filter(!is.na(!!min_age)) %>%
    ungroup() %>% 
    as.data.frame() %>% 
    arrange(Stock,BroodYear)
  
}

# Convert a brood table to a return table
brood_to_return <- function(bt){
  bt %>% group_by(Stock) %>% 
    pivot_longer(cols=contains("Age"), names_to="AgeNames", values_to="Return") %>% 
    mutate(Age=parse_number(AgeNames),
           ReturnYear=BroodYear+Age) %>% 
    filter(!is.na(Return)) %>% 
    select(Stock, ReturnYear, AgeNames, Return) %>% 
    pivot_wider(names_from=AgeNames, values_from=Return) 
}

# Prepare x and y vectors for models ####
prep_xy <- function(bt, ages) {
 
  # Make a list of dataframes for each age, 
  # with columns Stock, yt=log(age X returns) and xt=log(age X-1 returns)
  o <- lapply(ages, function(x){
    
        y_age <- paste0("Age", x)
        x_age <- paste0("Age", x - 1)

  bt %>% group_by(Stock) %>% 
        select(Stock, BroodYear, y_age, x_age) %>% 
        mutate_at(vars(y_age, x_age), ~log(.)) %>% 
        rename(yt=y_age, xt=x_age) %>% 
        ungroup() %>%
        mutate(model_Age=y_age) %>%
        select(Stock, BroodYear, model_Age, yt, xt) %>% 
        filter(!is.na(xt)) %>% 
        as.data.frame()

})
  
do.call(rbind,o)

}

# Create a list of functions to build DLM models ####
dlm_build_fns<- function(x.mat){
  
#   # Constant Intercept-only model
  constIntOnly <- function(parm,x.mat){
	parm <- exp(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], dW=0, addInt=FALSE))
    }

  # Random walk- time-varying intercept only
  tvIntOnly <- function(parm,x.mat){
	parm <- exp(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], dW=parm[2], addInt=FALSE))
  }

  # Time-varying Intercept and Slope
  tvIntSlope <- function(parm, x.mat){
	parm <- exp(parm)
	return(dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2], parm[3] )))
}

  # Time-varying slope
  tvSlope <- function(parm, x.mat){
	parm <- exp(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], dW=c(0, parm[2] )))
  }
  
  # Linear regression constant slope/intercept
  constLM <- function(parm, x.mat){
	parm <- exp(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], dW=c(0, 0)))
  }
  
  # Time-varying Cohort ratio (i.e., zero intercept model)
  tvCRzeroInt <- function(parm, x.mat){
	parm <- exp(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2]), addInt=FALSE))
  }
  
  # Constant cohort ratio (i.e., zero intercept model)
  constCRzeroInt <- function(parm, x.mat){
	parm <- exp(parm)
  return(dlmModReg(X=x.mat, dV=parm[1],dW=c(0), addInt=FALSE))
  }
  
  # Time-varying intercept, constant slope
  tvInt <- function(parm, x.mat){
    parm <- exp(parm)
    return(dlmModReg(X=x.mat, dV=parm[1],dW=c(parm[2], 0)))
  }
  
list(constIntOnly=constIntOnly, 
     tvIntOnly=tvIntOnly, 
     tvIntSlope=tvIntSlope, 
     tvSlope=tvSlope,
     tvInt=tvInt,
     constLM=constLM,
     tvCRzeroInt=tvCRzeroInt,
     constCRzeroInt=constCRzeroInt 
     )

}

# Estimate model params by MLE ####
est_parms <- function(npar, yt, xt, build_fns){
    
  fits <- lapply(1:length(npar), function(x){ 
          parm <- rep(0, npar[x])
          
        if(grepl(names(build_fns[x]),pattern="IntOnly")){
          xt <- rep(1,length(xt))}
          dlmMLE(y=yt, x.mat=xt, parm=parm, build=build_fns[[x]], hessian=T)
          
          })
     lapply(fits,"[[","par")
}

# Construct estimated models with the build functions and the MLE estimates ####
build_mods <- function(build_fns, fit_pars, xt){
            m <- lapply(1:length(build_fns), function(x){
            
              if(grepl(names(build_fns[x]),pattern="IntOnly")){
                xt <- rep(1,length(xt))}

                 build_fns[[x]](parm=fit_pars[[x]], x.mat=xt)}
            )
  names(m) <-  names(build_fns) 
  m
}

# Filter estimated models ####
filter_mods <- function(mods, yt){
      lapply(mods, function(x) { dlmFilter(y=yt, mod=x)})
}

# Smooth estimated models ####
smooth_mods <- function(mods, yt){
      lapply(mods, function(x) { dlmSmooth(y=yt, mod=x)})
}

# Calculate IC from models and data. ####
IC_table<-function(yt, mods, k_par){
  
 # Calc ICs for models
ICs <-   lapply(1:length(mods), function(x){
    
      # n data points
      n <- length(which(!is.na(yt)))
      
      # Negative log-likelihoods from fitted mods
      nLL <- c()
      nLL[x] <- dlmLL(y=yt, mods[[x]])
      
  # Calculate Information Criteria (ICs)
   AIC <- 2*k_par[x] +  2*nLL[x]
   
   AICc <- AIC + (2*k_par[x]^2 + 2*k_par[x]) / (n - k_par[x]-1)
   
   BIC  <- log(n) * k_par[x] + 2*nLL[x]
   
# Output a dataframe of the ICs for model x
  data.frame(AIC=AIC, AICc=AICc, BIC=BIC)
      }
    ) 

  # Name list elements by model names
  names(ICs) <- names(mods)
  
  #Rbind list elements intoa data frame, 
  do.call(rbind.data.frame, ICs) %>% 
    
    # Add model names and n parameters dataframe
    mutate(Model=rownames(.), npar=k_par) %>% 
    
    select(Model, npar, AIC, AICc, BIC) %>%
    
    #calculate deltaICs and model wts
    mutate(deltaAICc= AICc - min(AICc), 
           deltaBIC= BIC - min(BIC),
           AICc_wt= exp(-.5*deltaAICc) / sum(exp(-.5*deltaAICc)),
           BIC_wt= exp(-.5*deltaBIC) / sum(exp(-.5*deltaBIC))
           ) %>% 
    # Sort by lowest BIC
    arrange(AICc)
}

## Function to calculate intervals for forecasts ####
CIs <- function(fm,BY){
 
  # MAKE SURE TO SOURCE THIS from 'dlm_resid_fixed.R' FOR THE DEBUGGED VERSION
  sds <- residuals.dlmFiltered(fm, type="raw")$sd

# Forecast
  fore <- fm$f

  pred_mu <- exp(fore + 0.5 * sds^2)
  pred_sd <- sqrt(exp(2*fore+sds^2)*(exp(sds^2)-1))
  actual <- exp(fm$y)
  
  alpha <- c(0.2, 0.1, 0.05)
  
  interval_names <- paste0(rep(c("lcl_","ucl_"), each=length(alpha)), rep((1 - alpha)*100, times=2))
  
  qs <- alpha/2
  lls <- 0+qs
  uls <- 1-qs
  q <- c(lls,uls)
  
  lapply(seq_along(fore), function(x){
    #cls <- qlnorm(q, mean=fore[x], sd=sds[x])
    cls <- exp(qnorm(q, mean=fore[x], sd=sds[x]))
    names(cls) <- interval_names
    cls
    }) %>% do.call(rbind,.) %>% 
           cbind.data.frame(log_pred=fore, log_sd=sds, Pred_mu=pred_mu, Pred_med=exp(fore),actual,pred_sd=pred_sd,.) %>% 
           select(log_pred,log_sd, Pred_mu, Pred_med, pred_sd, actual, starts_with("lcl_"), starts_with("ucl_")) %>%
           mutate(APE=abs(Pred_med-actual)/actual,
                  SQE=(Pred_med-actual)^2,
                  Q=log(Pred_med/actual),
                  BroodYear=BY) %>% 
    # Could make this an argument to choose training set.
          tail(-10)
}

# Calculate error metrics- MdAPE, RMSE, MdLQ, MSA ####
err.metrix <- function(filtered_mods, BY,age){ 
  lapply(filtered_mods, CIs, BY) %>% 
  do.call(rbind.data.frame,.) %>% 
  mutate(Mod=rownames(.)) %>% 
  separate(Mod,into=c("Model", NA)) %>% 
  mutate(ReturnYear=BroodYear+age) %>% 
  group_by(Model) %>% 
    filter(is.finite(Pred_mu), !is.na(Pred_mu),ReturnYear>=2013) %>%  
  summarise(MdAPE=median(APE, na.rm=TRUE), 
            MAPE=mean(APE, na.rm=TRUE)*100, 
            RMSE=sqrt(mean(SQE,na.rm=TRUE)),
            MdLQ=median(Q, na.rm=TRUE),
            MSA=100*(exp(median(abs(Q),na.rm=TRUE))-1))}

# Wrapper to do forecasts and output list of results ####
do_forecasts <- function(df, table_type=c("brood", "return")) {

  table_type <- match.arg(table_type)
  if(table_type=="return"){bt <- return_to_brood(df)}
  else{bt <- df}

  stocks <- unique(bt$Stock)
  ages <- find_ages_to_forecast(bt)

lut <- ages %>% merge(stocks, stringsAsFactors=F) %>% select(stock=y,age=x)

out <- list()

# Vars specific to models ####
# Number of variance parameters in each model (Residual variance, Intercept and/or slope variance
  vpar <- c(1, 2, 3, 2, 2, 1, 2, 1)
  
# Number of coeffients for constant effects plus and initial states for time varying parms in each model
  nbetas <- c(1, 1, 2, 2, 2, 2, 1, 1)
  
  # Total parms estimated for each model to calc AIC
  k_par <- vpar-1 + nbetas
  
  # Prep xy data for models. Returns a list for each age class and stock
  o <- prep_xy(bt, ages=ages)
  
# Loop to do forecasts for all stocks/ages  
  for (i in 1:nrow(lut)){

  modAge <- paste0("Age", lut$age[i])
  stock <-as.character(lut$stock[i])

  df <- o %>% filter(Stock==lut$stock[i], model_Age==paste0("Age", lut$age[i]))
  BY <- df %>% pull(BroodYear)
  build_fns <- dlm_build_fns(x.mat=df$xt)
  fit_pars <- est_parms(npar=vpar, yt=df$yt, xt=df$xt, build_fns=build_fns)
  mods <- build_mods(build_fns=build_fns, fit_pars=fit_pars, xt=df$xt)
  filtered_mods <- filter_mods(mods, yt=df$yt)
  smoothed_mods <- smooth_mods(mods, yt=df$yt)

out[[i]] <- 
  list(
    
modSelTable=IC_table(df$yt, mods, k_par) %>% 
    left_join(err.metrix(filtered_mods,BY,age=lut$age[i]), by="Model") %>% 
              mutate(Stock=stock, Age=modAge) %>% 
  separate(Age, c(NA, "Age"), -1) %>% 
  mutate(Age=as.numeric(Age)),

Preds=lapply(filtered_mods, CIs, BY) %>% 
  do.call(rbind,.) %>% 
  mutate(Mod=rownames(.)) %>% 
  separate(Mod,into=c("Model",NA)) %>% 
  mutate(Stock=stock, Age=modAge) %>% 
  separate(Age,c(NA,"Age"),-1) %>% 
  mutate(Age=as.numeric(Age)) %>%  
  select(-APE,-SQE))

}

modsel <- lapply(out,"[[",1) %>% do.call(rbind,.)
preds <- lapply(out,"[[",2) %>% do.call(rbind,.) 

 o <- list(ModelSelectionResults=modsel, Predictions=preds)
 o
}
# TODO CLEAN UP THE DECIMAL PT BULLSHIT, CHANGE FORECASTS TO INTEGER.
# Make a model selection table for a stock and age class ####
make_mod_sel_tbl <- function(o, stock, age){

  o$Predictions %>% 
    filter(Stock==stock,Age==age,is.na(actual)) %>% 
    rename(Forecast=Pred_med) %>% #lapply(class)
  left_join(o$ModelSelectionResults,by=c("Model"="Model", "Stock"="Stock","Age"="Age" ))%>%
    group_by(Stock,Age) %>% 
    arrange(Age,AICc) %>% 
    select(-BIC,-BIC_wt,-MdAPE) %>% 
    mutate(AICc_wt=round(AICc_wt, 2)) %>% 
    mutate_at(vars("Forecast","lcl_95","ucl_95"),~round(.x,2)) %>% 
  flextable::flextable(col_keys=c("Model","npar","AICc","AICc_wt","MdLQ",
                                  
                                  "Forecast","MAPE","RMSE","MSA","lcl_95","ucl_95")) %>% 
  flextable::colformat_int(j=c("npar")) %>% 
  flextable::bold(j=1) %>% 
    flextable::bold(j=c(1,6),part="header") %>% 
  flextable::colformat_int(j="npar") %>% 
    flextable::width(j=6,width=.6) %>%
    flextable::width(j=1,width=1) %>% 
    flextable::width(j=2,width=0.5) %>% 
  flextable::colformat_double(j=c("AICc","Forecast","RMSE",
                                  "lcl_95","ucl_95"),digits=1) %>%
    #flextable::colformat_num(col_keys=c("Forecast","lcl_95","ucl_95"),digits=1) %>%  
  flextable::colformat_double(j=c("AICc_wt","MdLQ","MAPE","MSA"),digits=3) %>%  
    #flextable::colformat_num(col_keys=c("MdLQ"),digits=2) %>% 
  flextable::fontsize(j=2,size=8) %>%
  flextable::fontsize(size=8,part="body") %>% 
    flextable::fontsize(part="header",size=8) %>% 
    flextable::fontsize(j=7, size=8) %>% 
    flextable::padding(j=c(1,2), padding.right=5,part="header") %>% 
    flextable::padding(j=c(1,2), padding.right=5,part="body") %>% 
    flextable::padding(j=-c(1,2), padding.right=5,part="all") %>%
    flextable::padding(padding.top=6,part="all") %>% 
    flextable::padding(padding.bottom=6,part="all") %>% 
    flextable::align(align="center",part="all") %>% 
    flextable::autofit()
}

# Functionto calculate total forecast ####
mod_avg_totals <- function(o,stock){
 
o$Predictions %>% mutate(ReturnYear=BroodYear+Age) %>%  filter(Stock==stock) %>% 
  # Join to the Model selection 
  left_join(o$ModelSelectionResults, by=c("Model", "Stock", "Age")) %>% 
  select(Stock, Age,ReturnYear, Model, log_pred, log_sd, AICc_wt, BIC_wt,actual) %>% 
  group_by(Stock, Age, ReturnYear,actual) %>% 
  summarise(wtd_log_pred=sum(AICc_wt * log_pred),
            wtd_log_sd=sum(AICc_wt * sqrt(log_sd^2 + (log_pred - wtd_log_pred)^2))) %>%
  
        mutate(Pred_mu=exp(wtd_log_pred + 0.5*wtd_log_sd^2),
               Pred_med=exp(wtd_log_pred),
               lcl_95=qlnorm(.025,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               ucl_95=qlnorm(.975,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               lcl_80=qlnorm(0.10,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               ucl_80=qlnorm(0.90,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               lcl_50=qlnorm(0.25,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               ucl_50=qlnorm(0.75,meanlog=wtd_log_pred, sdlog=wtd_log_sd)
               
               ) %>% 
  
   select(Stock, Age, Pred_mu, Pred_med,ReturnYear,actual, contains("cl_")) %>% 
  ungroup %>%
  group_by(Stock,ReturnYear) %>%
    summarise (  Pred_med=sum(Pred_med),
                 actual=sum(actual)) %>% 
    ungroup %>% 
    mutate(log_resid=log(Pred_med)-log(actual),sd=lag(zoo::rollapplyr(log_resid,10,function(x){sqrt(sum(x^2)/(length(x)-1))},fill=NA)),

           lcl_95=exp(qnorm(.025,log(Pred_med),sd)),
           ucl_95=exp(qnorm(.975,log(Pred_med),sd)),
           lcl_80=exp(qnorm(.1,log(Pred_med),sd)),
           ucl_80=exp(qnorm(.9,log(Pred_med),sd)),
           lcl_50=exp(qnorm(.25,log(Pred_med),sd)),
           ucl_50=exp(qnorm(.75,log(Pred_med),sd))) %>% 
    select(Stock,ReturnYear,Forecast=Pred_med,lcl_95:ucl_50) %>% 

  mutate_if(is.numeric, round, 3) %>% as.data.frame
}

# Function to make model selection tables ####
# o is ouptut from fit/predict for-loop
mod_avg_totals_ft <- function(o, stock){

age_forecasts <- o$Predictions %>% 
  mutate(ReturnYear=BroodYear+Age) %>% 
  filter(Stock==stock) %>% #, is.na(actual)) %>%
  left_join(o$ModelSelectionResults, by=c("Model", "Stock", "Age")) %>%
  select(Stock,ReturnYear, actual,Age, Model, log_pred, log_sd, AICc_wt, BIC_wt) %>%
  group_by(Stock, Age,ReturnYear,actual) %>%
  summarise(wtd_log_pred=sum(AICc_wt * log_pred),
            wtd_log_sd=sum(AICc_wt * sqrt(log_sd^2 + (log_pred - wtd_log_pred)^2))) %>%

        mutate(Pred_mu=exp(wtd_log_pred + 0.5*wtd_log_sd^2),
               Pred_med=exp(wtd_log_pred),
               lcl_95=qlnorm(.025,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               ucl_95=qlnorm(.975,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               lcl_80=qlnorm(0.10,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               ucl_80=qlnorm(0.90,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               lcl_50=qlnorm(0.25,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               ucl_50=qlnorm(0.75,meanlog=wtd_log_pred, sdlog=wtd_log_sd),
               Age=as.character(Age)

               ) %>%
  
  bind_rows( group_by(.,Stock,ReturnYear) %>% 
              summarise (  Pred_med=sum(Pred_med),
                           actual=sum(actual)) %>% 
                ungroup %>% 
                mutate(
                  log_resid=log(Pred_med)-log(actual),
                  sd=lag(zoo::rollapplyr(log_resid,10,function(x){sqrt(sum(x^2)/(length(x)-1))},fill=NA)),

              lcl_80=exp(qnorm(.1,log(Pred_med),sd)),
              ucl_80=exp(qnorm(.9,log(Pred_med),sd)),
              lcl_95=exp(qnorm(.025,log(Pred_med),sd)),
              ucl_95=exp(qnorm(.975,log(Pred_med),sd)),
              # lcl_60=exp(qnorm(.2,log(Pred_med),sd)),
              # ucl_60=exp(qnorm(.8,log(Pred_med),sd)),
              lcl_50=exp(qnorm(.25,log(Pred_med),sd)),
              ucl_50=exp(qnorm(.75,log(Pred_med),sd))) %>% 
                           # Pred_mu=sum(Pred_med),
                           # lcl_80=sum(lcl_80),
                           # ucl_80=sum(ucl_80),
                           # lcl_95=sum(lcl_95),
                           # ucl_95=sum(ucl_95),
                           # lcl_60=sum(lcl_60),
                           # ucl_60=sum(ucl_60)) 
              mutate(Age="Total") ) %>% 
            
  group_by(Stock, Age,ReturnYear,actual) %>%
  mutate(SQE=(Pred_med-actual)^2,
         
         APE=abs(Pred_med-actual)/actual,
         SA=abs(log(actual/Pred_med))) %>% 
group_by(Stock, Age) %>% 
         mutate(
           across(c(SQE,APE,SA),function(x)ifelse(ReturnYear>=2013,x,NA)),
         RMSE=sqrt(mean(SQE,na.rm=TRUE)),
         MAPE=mean(APE, na.rm=TRUE)*100,
         MSA = 100*(exp(mean(SA, na.rm=TRUE))-1)
         ) %>% 
  ungroup %>% 
    filter(is.na(actual))%>%
  select(Stock, Age, Pred_mu, Pred_med,actual,RMSE,MAPE,MSA,ReturnYear, contains("cl_"))

combined <- age_forecasts %>%
  # group_by(Stock) %>%
  # summarise (  Forecast=sum(Pred_med),
  #              MAPE=sum(MAPE),
  #              RMSE=sum(RMSE),
  #                     lcl_80=sum(lcl_80),
  #                     ucl_80=sum(ucl_80),
  #                     lcl_95=sum(lcl_95),
  #                     ucl_95=sum(ucl_95),
  #                     lcl_60=sum(lcl_60),
  #                     ucl_60=sum(ucl_60)) %>% 
  # mutate(Age="Total") %>% 
  # bind_rows(age_forecasts %>% 
              # mutate(Age=as.character(Age)) %>%  
              select(Stock,Forecast=Pred_med,Age,MAPE,RMSE,MSA,contains("cl_")) %>% 
  arrange(Age) %>% 
 # mutate_if(is.numeric, round, 3) %>%
  as.data.frame() 
  
combined %>%  select(-Stock) %>% 
  flextable::flextable(col_keys=c("Age","Forecast","MAPE","RMSE","MSA", 
                                  "lcl_95","ucl_95","lcl_50","ucl_50")) %>% 
  merge_v(j=1) %>% 
    #flextable::colformat_int(j=c("lcl_95","ucl_95","lcl_80","ucl_80")) %>% 
  flextable::colformat_double(j=c("Forecast","RMSE",
                                  "lcl_95","ucl_95","lcl_50","ucl_50"),digits=1) %>%
  #flextable::colformat_num(col_keys=c("Forecast","lcl_95","ucl_95"),digits=1) %>%  
  flextable::colformat_double(j=c("MAPE","MSA"),digits=3) %>%  
  border(i=nrow(combined),border.top=officer::fp_border(width=2),border.bottom=officer::fp_border(width=2)) %>% 
    flextable::autofit() 
}

# Plot returns by age class for a stock ####
# Args: o= the output list from model fitting (analysis_script), stock will be looped over based on unique stocks in bt
plot_returns <- function(o, stock){
 dat <-  o$Predictions %>%
    select(Stock, BroodYear, Age, actual) %>% 
    distinct() %>% mutate(ReturnYear=BroodYear+Age) %>% 
    select(-BroodYear) %>% 
    filter(!is.na(actual),Stock==stock) 
 min_yr<-min(dat$ReturnYear)
 
 dat<-dat%>% 
   bind_rows(df %>% filter(Stock==stock,ReturnYear>=min_yr) %>% select(ReturnYear,actual=Age3,Stock) %>% mutate(Age=3))
 
 # Plot variables depending on dat
 
  # Facet plot with a row for each age class
  facet_rows <- length(unique(dat$Age))
  
  # TAC brood tables are N, Fall brood tables are thousands(N/1000)
  ylabel <- ifelse(any(dat$actual>2000,na.rm=TRUE), "Run size", "Run size (thousands)")
  
dat %>%   
  ggplot(aes(x=ReturnYear,y=actual,col=as.factor(Age)))+
    geom_line(lwd=1.2,show.legend=FALSE)+
    geom_point(col="black")+
     facet_wrap(vars(Age),nrow=facet_rows, strip.position="top", scales="free_y")+
  scale_y_continuous(labels=function(x)format(x, scientific=FALSE, big.mark=","))+
  theme(plot.title=element_text(face="bold",size=16,hjust=.5),
        plot.subtitle=element_text(hjust=0.5),
        strip.text=element_text(face="bold",size=12))+
        ylab(ylabel)+
        xlab("Return year")
}


# Plot predictions with size proportional to AIC wt. ####

plot_preds <- function(o, stock, age){
  
actuals <- o$Predictions %>% 
  filter(Stock==stock,Model=="constIntOnly") %>% 
      mutate(ReturnYear=BroodYear+Age) %>% 
      group_by(Stock, Model, ReturnYear,Age) %>% 
      summarise(Return=sum(actual)) %>% 
      filter(Age==age)

preds <- o$Predictions %>% filter(Stock==stock, is.na(actual)) %>%  
  left_join(o$ModelSelectionResults,by=c("Stock","Model","Age")) %>%#filter(AICc_wt>0) %>% 
  select(Model,Age,AICc_wt,starts_with("lcl"),starts_with("ucl"),Pred_med,BroodYear) %>%
  mutate(ReturnYear=BroodYear+Age) %>% 
  filter(Age==age) %>% 
  mutate(plot_helper=factor(x=AICc_wt, levels=sort(unique(AICc_wt),decreasing=FALSE),ordered=TRUE),
         scaled_AIC_wt=AICc_wt/max(AICc_wt)*6)

# Variables for plot that depend on data
# labels for models in the legend
modlabs <- preds %>% select(Model,AICc_wt) %>% arrange(AICc_wt) %>% pull(Model)
dotsizes <- preds$scaled_AIC_wt %>% sort()

# Y axis label
ylabel <- ifelse(any(actuals$Return>2000,na.rm=TRUE),"Run size","Run size (thousands)")

ggplot(data=actuals)+
      geom_line(aes(x=ReturnYear,y=Return,linetype=""),lwd=1,show.legend=NA)+
      facet_wrap(vars(Age),nrow=1,strip.position="top",scales="free_y")+
      ggtitle(stock)+
      theme(plot.title=element_text(face="bold",size=16,hjust=.5),
        plot.subtitle=element_text(hjust=0.5),
        strip.text=element_text(face="bold",size=12))+
        ylab(ylabel)+
        xlab("Return year")+
  geom_point(data=preds ,aes(x=ReturnYear, y=Pred_med, color=plot_helper, size=plot_helper),alpha=.5,show.legend=NA)+
  facet_wrap(vars(Age),nrow=1,strip.position="top",scales="free_y")+
  scale_size_manual(name="Model forecast", labels=modlabs, values=1:8)+
  scale_shape_manual(values=c(rep(23,8),NULL))+
  scale_color_discrete(name="Model forecast",labels=modlabs)+
  scale_linetype_manual(name="",labels="Observed",values="solid")+
  scale_y_continuous(labels=function(x)format(x,scientific=FALSE,big.mark=","))+
  theme(legend.key=element_blank())+
  #scale_x_continuous(limits=c(2000, NA))+
  coord_cartesian(xlim=c(2004, 2023))+
   guides(size = guide_legend(order=2,reverse=TRUE),
          color=guide_legend(order=2,reverse=TRUE),
          linetype=guide_legend(order=1,label.position="right"))
          
}

# Plot Total forecast ####
plot_tots <- function(o, stock){

 actuals_total <- o$Predictions %>% filter(Stock==stock,Model=="constIntOnly") %>% 
      mutate(ReturnYear=BroodYear+Age) %>% 
      group_by(Stock, Model, ReturnYear) %>% 
      summarise(Return=sum(actual)) %>% ungroup %>%  as.data.frame

 forecast_yr <- max(actuals_total$ReturnYear)
 
ylabel <- ifelse(any(actuals_total>2000,na.rm=TRUE),"Run size", "Run size (thousands)")

tot_mod_avg <-  mod_avg_totals(o,stock)

  actuals_total %>% 
    ggplot(data=.,aes(x=ReturnYear,y=Return,col="Actual")) +
    geom_point(size=2,show.legend=NA)+geom_line(aes(col="Actual"),lwd=.5,lty=2,show.legend=TRUE)+
  
    # MAKE YEAR max(ReturnYear +1)
    geom_point(data=tot_mod_avg %>% filter(ReturnYear<forecast_yr),aes(x=ReturnYear, y=Forecast,col="Predicted"),size=2,inherit.aes=F,show.legend=TRUE)+
    scale_color_manual(name="",values=c("black","red","red"),guide=guide_legend(override.aes=list(linetype=c("dashed","blank"))))+
    theme(legend.key=element_blank())+
   geom_point(data=tot_mod_avg %>% filter(ReturnYear==forecast_yr),aes(x=ReturnYear, y=Forecast),size=3.5,inherit.aes=FALSE,col="red")+
    geom_errorbar(data=tot_mod_avg %>% filter(ReturnYear==forecast_yr), aes(x=ReturnYear,ymin=lcl_95,ymax=ucl_95),width=.75,size=.75,col="red",inherit.aes=F,show.legend=F)+
    geom_errorbar(data=tot_mod_avg %>% filter(ReturnYear==forecast_yr), aes(x=ReturnYear,ymin=lcl_50,ymax=ucl_50),width=.5,size=.5,col="red",inherit.aes=F,show.legend=F)+
    
  ylab(ylabel)+
    scale_y_continuous(labels=function(x)format(x,scientific=FALSE,big.mark=",")) + 
    scale_x_continuous(minor_breaks=seq(1900,2022,1))+
  xlab("Return year") +
  ggtitle(label=paste0(stock, " total adults"),subtitle= paste0("Ages ", paste(unique(o$Predictions$Age),collapse=",")))+
      theme(plot.title=element_text(face="bold",size=16,hjust=.5),
        plot.subtitle=element_text(hjust=0.5),
        strip.text=element_text(face="bold",size=12))
}
