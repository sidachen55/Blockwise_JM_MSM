##### update sim2 motivated by cprd last updated 2023.2.4 (sbp as biomarker)

library(simsurv)
library(mvtnorm)

n<-1000
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

set.seed(task_id)
# define transition specific hazard functions
haz_01 <- function(t, x, betas, ...) {
  betas[["delta"]]*(t^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*t))
}

haz_02 <- function(t, x, betas, ...) {
  betas[["delta"]]*(t^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*t))
}

haz_04 <- function(t, x, betas, ...) {
  betas[["delta"]]*(t^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*t))
}

# t_e, time entered the current state
haz_13 <- function(t, x, betas, t_e, ...) {
  betas[["delta"]]*((t)^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*(t+t_e)))
}

haz_14 <- function(t, x, betas, t_e, ...) {
  betas[["delta"]]*((t)^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*(t+t_e)))
}

haz_23 <- function(t, x, betas, t_e, ...) {
  betas[["delta"]]*((t)^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*(t+t_e)))
}

haz_24 <- function(t, x, betas, t_e, ...) {
  betas[["delta"]]*((t)^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*(t+t_e)))
}

haz_34 <- function(t, x, betas, t_e, ...) {
  betas[["delta"]]*((t)^(betas[["delta"]]-1))*exp(betas[["gamma_0"]]+betas[["gamma_1"]]*x[["w1"]]+betas[["alpha"]]*(betas[["beta_0i"]]+betas[["beta_1i"]]*(t+t_e)))
}

####################################################################
#### basic simulation settings 
# generate random effects (motivated from cprd hba1c)
rho <- -0.5
b_corrmat <- matrix(c(1, rho, rho, 1), 2, 2)
b_sds <- c(0.7, 0.1)
b_means <- rep(0, 2)
b_z <- MASS::mvrnorm(n = n, mu = b_means, Sigma = b_corrmat)
b <- sapply(1:length(b_sds), FUN = function(x) b_sds[x] * b_z[,x])
# LMM parameters (motivated from cprd hba1c)
beta_0<- 0 # LMM intercept (in our control!)
beta_1<- -0.01 # LMM slope
beta_0i<-rep(beta_0, n) + b[, 1] # subject-specific intercept
beta_1i<-rep(beta_1, n) + b[, 2] # subject-specific slope 
# parameters in MSM transitions (alpha in our control; other pars motivated from MSM fitted to CPRD: CV, T2, MH, D)
betas_01 <- data.frame(delta = rep(1.19, n), gamma_0 = rep(log(0.011), n),
                       gamma_1 = rep(0.63, n), alpha = rep(0.13, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_02 <- data.frame(delta = rep(0.95, n), gamma_0 = rep(log(0.017), n),
                       gamma_1 = rep(-0.36, n), alpha = rep(0.07, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_04 <- data.frame(delta = rep(1.63, n), gamma_0 = rep(log(0.0025), n),
                       gamma_1 = rep(1.29, n), alpha = rep(-0.64, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_13 <- data.frame(delta = rep(0.87, n), gamma_0 = rep(log(0.02), n),
                       gamma_1 = rep(-0.35, n), alpha = rep(0.01, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_14 <- data.frame(delta = rep(1.48, n), gamma_0 = rep(log(0.014), n),
                       gamma_1 = rep(0.96, n), alpha = rep(-0.61, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_23 <- data.frame(delta = rep(1.21, n), gamma_0 = rep(log(0.01), n),
                       gamma_1 = rep(0.45, n), alpha = rep(0.08, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_24 <- data.frame(delta = rep(1.51, n), gamma_0 = rep(log(0.004), n),
                       gamma_1 = rep(1.4, n), alpha = rep(-0.85, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)

betas_34 <- data.frame(delta = rep(1.36, n), gamma_0 = rep(log(0.018), n),
                       gamma_1 = rep(1.07, n), alpha = rep(-0.65, n),
                       beta_0i = beta_0i, beta_1i = beta_1i)



# baseline covariates in MSM 
age   <- c(runif(n*0.55, min = 18, max = 65), runif(n*0.3, min = 65, max = 80), runif(n*0.15, min = 80, max = 90))
age <- (age - mean(age))/sd(age)
covdat <- data.frame(w1 = age)



#### generate event times
# initiate vectors to save true event times (clock-reset)
trueT01 <- trueT02 <- trueT04 <- trueT13 <- trueT14 <- trueT23 <- trueT24 <- trueT34 <- numeric(n)

# sample censoring times (non-informative)
mean.Cens <- 12
Ctimes <- runif(n, 13, 2 * mean.Cens)

# simulate time-to-event data (for each transition)
for (i in 1:n) {
  # Transition 0->1
  trueT01[i] <- simsurv(hazard = haz_01, x = covdat[i,,drop=FALSE], betas = betas_01[i,], maxt = 100)$eventtime
  # Transition 0->2
  trueT02[i] <- simsurv(hazard = haz_02, x = covdat[i,,drop=FALSE], betas = betas_02[i,], maxt = 100)$eventtime
  
  trueT04[i] <- simsurv(hazard = haz_04, x = covdat[i,,drop=FALSE], betas = betas_04[i,], maxt = 100)$eventtime
  # 
  if (Ctimes[i]<min(trueT01[i],trueT02[i],trueT04[i])) {
    trueT13[i]<-500; trueT14[i]<-500; trueT23[i]<-500; trueT24[i]<-500; trueT34[i]<-500
  } else if (trueT04[i]<min(trueT01[i],trueT02[i])) {
    trueT13[i]<-500; trueT14[i]<-500; trueT23[i]<-500; trueT24[i]<-500; trueT34[i]<-500
  } else {
    if (trueT01[i]<trueT02[i]) {
      trueT13[i] <- simsurv(hazard = haz_13, x = covdat[i,,drop=FALSE], betas = betas_13[i,],t_e=trueT01[i], maxt = 100)$eventtime
      trueT14[i] <- simsurv(hazard = haz_14, x = covdat[i,,drop=FALSE], betas = betas_14[i,],t_e=trueT01[i], maxt = 100)$eventtime
      if ((trueT13[i]<trueT14[i]) & (trueT01[i]+trueT13[i]<Ctimes[i])) {
        trueT34[i] <- simsurv(hazard = haz_34, x = covdat[i,,drop=FALSE], betas = betas_34[i,],t_e=trueT01[i]+trueT13[i], maxt = 100)$eventtime
        trueT23[i]<-500; trueT24[i]<-500
      } else {
        trueT23[i]<-500; trueT24[i]<-500; trueT34[i]<-500
      }
    } else {
      trueT23[i] <- simsurv(hazard = haz_23, x = covdat[i,,drop=FALSE], betas = betas_23[i,],t_e=trueT02[i], maxt = 100)$eventtime
      trueT24[i] <- simsurv(hazard = haz_24, x = covdat[i,,drop=FALSE], betas = betas_24[i,],t_e=trueT02[i], maxt = 100)$eventtime
      if ((trueT23[i]<trueT24[i]) & (trueT02[i]+trueT23[i]<Ctimes[i])) {
        trueT34[i] <- simsurv(hazard = haz_34, x = covdat[i,,drop=FALSE], betas = betas_34[i,],t_e=trueT02[i]+trueT23[i], maxt = 100)$eventtime
        trueT13[i]<-500; trueT14[i]<-500
      } else {
        trueT13[i]<-500; trueT14[i]<-500; trueT34[i]<-500
      }
    }
  }
  print(i)
}

data_mstate <- data.frame('id' = 1:n, 'trueT01' = trueT01, 'trueT02' = trueT02,　'trueT04' = trueT04, 'trueT13' = trueT13, 'trueT14' = trueT14,'trueT23' = trueT23,'trueT24' = trueT24, 'trueT34' = trueT34,
                          'Ctimes' = Ctimes, 'X1' = covdat$w1)
data_mstate_split.by.id <- split(data_mstate, data_mstate$id)

ms_arrange <- function (x) {
  if (x$Ctimes < min(x$trueT01, x$trueT02, x$trueT04)) {
    x_new <- data.frame('id' = rep(x$id, 3), 'from_state' = rep(0, 3), 'to_state' = c(1,2,4), 
                        'transition' = 1:3, 'Tstart' = rep(0, 3), 'Tstop' = x$Ctimes, 'status' = rep(0, 3), 
                        'X1' = x$X1)
  } else if (x$trueT04 < min(x$trueT01, x$trueT02)) {
    x_new <- data.frame('id' = rep(x$id, 3), 'from_state' = rep(0, 3), 'to_state' = c(1,2,4), 
                        'transition' = 1:3, 'Tstart' = rep(0, 3), 'Tstop' = rep(x$trueT04,3), 'status' = c(0,0,1), 
                        'X1' = x$X1)
  } else {
    if (x$trueT01 < x$trueT02) {
      if(x$Ctimes<min(x$trueT01+x$trueT13,x$trueT01+x$trueT14)){
        x_new <- data.frame('id' = rep(x$id, 5), 'from_state' = c(0, 0, 0, 1,1), 'to_state' = c(1,2,4,3,4), 
                            'transition' = 1:5, 'Tstart' = c(rep(0, 3), rep(x$trueT01,2)), 
                            'Tstop' = c(rep(x$trueT01, 3), rep(x$Ctimes,2)), 'status' = c(1, 0, 0,0,0), 
                            'X1' = x$X1)
      } else if (x$trueT14<x$trueT13) {
        x_new <- data.frame('id' = rep(x$id, 5), 'from_state' = c(0, 0, 0, 1,1), 'to_state' = c(1,2,4,3,4), 
                            'transition' = 1:5, 'Tstart' = c(rep(0, 3), rep(x$trueT01,2)), 
                            'Tstop' = c(rep(x$trueT01, 3), rep(x$trueT01+x$trueT14,2)), 'status' = c(1, 0, 0,0,1), 
                            'X1' = x$X1)
      } else {
        if (x$Ctimes < x$trueT01+x$trueT13+x$trueT34) {
          x_new <- data.frame('id' = rep(x$id, 6), 'from_state' = c(0, 0, 0, 1,1,3), 'to_state' = c(1,2,4,3,4,4), 
                              'transition' = c(1:5,8), 'Tstart' = c(rep(0, 3), rep(x$trueT01,2), x$trueT01+x$trueT13), 
                              'Tstop' = c(rep(x$trueT01, 3), rep(x$trueT01+x$trueT13,2), x$Ctimes), 'status' = c(1, 0, 0, 1,0,0), 
                              'X1' = x$X1)
        } else {
          x_new <- data.frame('id' = rep(x$id, 6), 'from_state' = c(0, 0, 0, 1,1,3), 'to_state' = c(1,2,4,3,4,4), 
                              'transition' = c(1:5,8), 'Tstart' = c(rep(0, 3), rep(x$trueT01,2), x$trueT01+x$trueT13), 
                              'Tstop' = c(rep(x$trueT01, 3), rep(x$trueT01+x$trueT13,2), x$trueT01+x$trueT13+x$trueT34), 'status' = c(1, 0, 0, 1,0,1), 
                              'X1' = x$X1)
        }
      }
    } else {
      if (x$Ctimes<min(x$trueT02+x$trueT23,x$trueT02+x$trueT24)) {
        x_new <- data.frame('id' = rep(x$id, 5), 'from_state' = c(0, 0, 0, 2,2), 'to_state' = c(1,2,4,3,4), 
                            'transition' = c(1,2,3,6,7), 'Tstart' = c(rep(0, 3), rep(x$trueT02,2)), 
                            'Tstop' = c(rep(x$trueT02, 3), rep(x$Ctimes,2)), 'status' = c(0, 1, 0,0,0), 
                            'X1' = x$X1)
      } else if (x$trueT24<x$trueT23) {
        x_new <- data.frame('id' = rep(x$id, 5), 'from_state' = c(0, 0, 0, 2,2), 'to_state' = c(1,2,4,3,4), 
                            'transition' = c(1,2,3,6,7), 'Tstart' = c(rep(0, 3), rep(x$trueT02,2)), 
                            'Tstop' = c(rep(x$trueT02, 3), rep(x$trueT02+x$trueT24,2)), 'status' = c(0, 1, 0,0,1), 
                            'X1' = x$X1)
      } else {
        if (x$Ctimes < x$trueT02+x$trueT23+x$trueT34) {
          x_new <- data.frame('id' = rep(x$id, 6), 'from_state' = c(0, 0, 0, 2,2,3), 'to_state' = c(1,2,4,3,4,4), 
                              'transition' = c(1,2,3,6,7,8), 'Tstart' = c(rep(0, 3), rep(x$trueT02,2),x$trueT02+x$trueT23), 
                              'Tstop' = c(rep(x$trueT02, 3), rep(x$trueT02+x$trueT23,2),x$Ctimes), 'status' = c(0, 1, 0, 1,0,0), 
                              'X1' = x$X1)
        } else {
          x_new <- data.frame('id' = rep(x$id, 6), 'from_state' = c(0, 0, 0, 2,2,3), 'to_state' = c(1,2,4,3,4,4), 
                              'transition' = c(1,2,3,6,7,8), 'Tstart' = c(rep(0, 3), rep(x$trueT02,2),x$trueT02+x$trueT23), 
                              'Tstop' = c(rep(x$trueT02, 3), rep(x$trueT02+x$trueT23,2),x$trueT02+x$trueT23+x$trueT34), 'status' = c(0, 1, 0, 1,0,1), 
                              'X1' = x$X1)
        }
      }
    }
  }
}

data_mstate_split.by.id <- lapply(data_mstate_split.by.id, ms_arrange)
data_mstate <- do.call(rbind, data_mstate_split.by.id)
data_mstate$transition <- factor(data_mstate$transition)     
data_mstate$years<-data_mstate$Tstop-data_mstate$Tstart

#### simulate longitudinal data (increasing visiting freq)
DBRO_1 <- 2.6 # distance between repeated observations for block 1
DBRO_2 <- 2 # distance between repeated observations after #1 trans
DBRO_3 <- 1.2
sigma.re = 0.75 # measurement error sd;  
beta.tilde = c(beta_0,beta_1)
visits.out <- vector()
times<-tapply(data_mstate$Tstop, data_mstate$id, max)
ID <- longit.out <- vector()
for(i in 1:n){
  if (max(data_mstate$Tstart[which(data_mstate$id==i)])==0) {
    visits <- seq(0,times[i], by = DBRO_1) # visit time t_ij
  } else if (length(data_mstate$id[which(data_mstate$id==i)])==5) {
    visits <- c(seq(0,data_mstate$Tstop[which(data_mstate$id==i)[1]], by = DBRO_1),seq(data_mstate$Tstop[which(data_mstate$id==i)[1]],times[i], by = DBRO_2)) ### implication: 
  } else {
    visits <- c(seq(0,data_mstate$Tstop[which(data_mstate$id==i)[1]], by = DBRO_1),seq(data_mstate$Tstop[which(data_mstate$id==i)[1]],data_mstate$Tstop[which(data_mstate$id==i)[4]], by = DBRO_2),seq(data_mstate$Tstop[which(data_mstate$id==i)[4]],times[i], by = DBRO_3))
  }
  yt <- beta.tilde[1] + beta.tilde[2]*visits +b[i,1] + b[i,2]*visits + rnorm(length(visits),0,sigma.re)
  longit.out <- c(longit.out,yt)
  ID <- c(ID,rep(i,length(visits)))
  visits.out <- c(visits.out,visits)
}
long.proc <- as.matrix(cbind(longit.out,ID)) # Longitudinal process
obj <- list(long.proc,visits.out)
names(obj) <- c("longitudinal","visits")

save.image(paste0("data_",task_id,".RData"))