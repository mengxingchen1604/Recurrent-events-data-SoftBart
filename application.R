library(MASS)
library(msm)
library(glmnet)
library(SoftBart)
library(diversitree)
library(extraDistr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)
library(writexl)
#competing methods
library(recforest)
source("Rcpp proportional intensity function.R")
source("Martingale Residual for RecSBART.R")
source('cumulative test for application RecSBART.R')
source('RecSBART bayesian margin.R')
#dataset
library(frailtypack)
data("readmission")

#set up dataset
readmission$gender<-as.numeric(readmission$sex=="Male")
readmission$group<-as.numeric(readmission$chemo=="Treated")
readmission$dukes1<-as.numeric(readmission$dukes)

readmission$charlson_base <- as.numeric(ave(readmission$charlson, 
                                 readmission$id, 
                                 FUN = function(x) rep(x[1], length(x))))
#1-2 = 1.5
readmission$charlson_base<-sapply(readmission$charlson_base,function(x)if(x==2){x=1.5} else{x})
#A-B = 0.5
readmission$dukes1<-sapply(readmission$dukes1,function(x)if(x==1){x=0.5} else{x})

#convert dataset

#X is total for bayesian methods
X <- matrix(NA, nrow = 403, ncol = 4)
for (i in 1:403) {
  m <- subset(readmission, id == i)
  X[i, ] <- as.numeric(m[1, 12:15])
}

ID<-sort(unique(readmission$id))
readmission$terminal<-ave(readmission$t.stop,readmission$id,FUN = max)
t<- subset(readmission, t.stop==terminal)
terminal<-as.numeric(t$terminal)


recurrent_event <- lapply(split(readmission, readmission$id), function(x) {
  as.numeric(subset(x, event == 1)$t.stop)
})

recurrent_event<-lapply(recurrent_event,function(i) i/max(terminal))
terminal<- terminal/max(terminal)

#for recforest method
forest_train<-readmission[,c(1,3,4,6,12:15)]

recforest_martingale_residual<-train_forest(data = forest_train,
                                            id_var = 'id',
                                            covariates = c('gender','group','dukes1','charlson_base'),
                                            event = 'event',
                                            time_vars = c('t.start','t.stop'),
                                            death_var = NULL,
                                            n_trees = 50,
                                            n_bootstrap = 100,
                                            mtry = 4,
                                            minsplit = 10,
                                            nodesize = 100,
                                            method = "NAa",
                                            min_score = 100,
                                            max_nodes = 100,
                                            seed = 111,
                                            parallel = FALSE,
                                            verbose = FALSE)

expected_number_recforest<-predict(recforest_martingale_residual,
                                   newdata=forest_train,
                                   id_var = "id",
                                   covariates = c('gender','group','dukes1','charlson_base'),
                                   event = 'event',
                                   time_vars = c('t.start','t.stop'))
forest_train$cum<-expected_number_recforest
forest_train$cumsum<-ave(forest_train$event,forest_train$id,FUN = function(x) cumsum(x))
forest_train$martingale<- forest_train$cumsum-forest_train$cum

#RecSBART

#martingale residual
readmission_softbart<-fit_RecSBART_martingale(X_train = X,
                                            X_test = X,
                                            recurrent_train = recurrent_event,
                                            recurrent_test = recurrent_event,
                                            terminal_train = terminal,
                                            terminal_test = terminal,
                                            num_burn = 2500,
                                            num_thin = 1,
                                            num_save = 2500)

true_event<-lapply(recurrent_event,function(i) {if(length(i)==0){0} else {as.numeric(c(1:length(i),length(i)))}})

martingale_estimate<-readmission_softbart$martingale_mean

martingale_residual <- mapply(function(x, y) {
  x-y
}, true_event, martingale_estimate)

#subgroup martingale
softbart_readmission<-as.data.frame(forest_train)
softbart_readmission$m<-unlist(martingale_residual)

softbart_readmission<-read.csv("softbart_readmission_martingale.csv")


#proportional intensity model
Intensity_prop<- fit_proportional_intensity(X_train = X,
                                            X_test = X,
                                            recurrent_train = recurrent_event,
                                            recurrent_test = recurrent_event,
                                            terminal_train = terminal,
                                            terminal_test = terminal,
                                            num_burn = 2500,
                                            num_thin = 1,
                                            num_save = 2500)

beta_estimate<-colMeans(Intensity_prop$beta)
W_estimate<-colMeans(Intensity_prop$W)
eta_estimate<-mean(Intensity_prop$eta)
Lambda_estimate<-colMeans(Intensity_prop$Lambda0)

time_point <- c(0, unique(sort(unlist(recurrent_event))))
cum_Lambda <- cumsum(Lambda_estimate)
last_cum <- tail(cum_Lambda, 1)

theta <- as.vector(exp(X %*% beta_estimate))
Wv <- W_estimate

n <- length(recurrent_event)
M <- vector("list", n)

for (i in seq_len(n)) {
  rec <- recurrent_event[[i]]
  nEvent <- length(rec)
  
  if (nEvent == 0L) {
    M[[i]] <- - theta[i] * last_cum * Wv[i]
  } else {
    idx <- findInterval(rec, time_point, left.open = TRUE)
    cum_vals <- ifelse(idx == 0L, 0, cum_Lambda[idx])
    
    event_resid <- seq_len(nEvent) - (theta[i] * Wv[i]) * cum_vals
    terminal_resid <- nEvent - theta[i] * last_cum * Wv[i]
    
    M[[i]] <- c(event_resid, terminal_resid)
  }
}

prop_martingale<- readmission[,c(1,3,4,6,12:15)]
prop_martingale$cumnum <- ave(prop_martingale$event,prop_martingale$id,FUN = function(x) cumsum(x))
prop_martingale$m <- unlist(M)

