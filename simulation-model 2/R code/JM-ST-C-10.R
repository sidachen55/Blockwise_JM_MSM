###################################################
# JM-ST-C for Ferrer's sim model 
# Trans=9
###################################################
library(simsurv)
library(mvtnorm)
library(splines)
library(MASS)
library(rstan)

####################
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

set.seed(task_id)

input_file <- paste0("data_", task_id, ".RData")
load(input_file)
####################
###################
# B-splines for modelling the baseline transition
T_max <- 24
inknots <- c(7.127)

B_1 <- function(t) {
  result <- spline.des(knots = c(rep(0,4), inknots, rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,1]
  return(result)
}
B_2 <- function(t) {
  result <- spline.des(knots = c(rep(0,4), inknots, rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,2]
  return(result)
}
B_3 <- function(t) {
  result <- spline.des(knots = c(rep(0,4), inknots, rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,3]
  return(result)
}
B_4 <- function(t) {
  result <- spline.des(knots = c(rep(0,4), inknots, rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,4]
  return(result)
}
B_5 <- function(t) {
  result <- spline.des(knots = c(rep(0,4), inknots, rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,5]
  return(result)
}

f1 <- function (t) {
  return ((1+t)^(-1.2)-1)
}

f2 <- function (t) {return (t)}


################################################################################
# data preparation 
################################################################################
cr_id <- data_mstate$id[which(data_mstate$transition==10)]
data_UJ <- data_mstate[data_mstate$transition==10,]
L_data <- data.frame(ID=ID,long.data=longit.out,visits=visits.out)
L_data_cr<-L_data[L_data$ID %in% cr_id,]
L_data_cr$Tstop <- rep(data_UJ$Tstop,as.vector(table(L_data_cr$ID)))
L_data_cr$Tstart <- rep(data_UJ$Tstart,as.vector(table(L_data_cr$ID)))
L_data_split.by.id <- split(L_data_cr, L_data_cr$ID)
cr_extract <- function (x) {
  if (any(x$visits >= x$Tstart & x$visits < x$Tstop)) {
    x_new <- x[(x$visits >= x$Tstart) & (x$visits < x$Tstop), ]
  } else {
    x_new <- x[x$visits < x$Tstop, ]
  }
}
L_data_split.by.id <- lapply(L_data_split.by.id, cr_extract)
data_lg <- do.call(rbind, L_data_split.by.id)


  # relabel
  data_lg$ID<-rep(1:length(cr_id),as.vector(table(data_lg$ID)))
  data_UJ$id<-rep(1:length(cr_id),each=1)
  data_UJ$transition<-as.numeric(data_UJ$transition)
  data_UJ$transition<-replace(data_UJ$transition,data_UJ$transition==c(10),c(1))
  
  ############################################
  # Required quantities for model fitting
  
  n <- length(cr_id)              # total number of subjects
  y <- data_lg$long.data          # longitudinal outcomes
  ID <- data_lg$ID                # patient IDs
  visits <- data_lg$visits         # visit times for repeated observations
  N <- length(y)                     # total number of longitudinal outcomes
  
  
  library(statmod)
  glq <- gauss.quad(15, kind = "legendre")
  xk <- glq$nodes   # nodes
  wk <- glq$weights # weights
  
  Nevents<-sum(data_UJ$status==1)
  Qn<-length(data_UJ$status)*15
  nrow_q<-Nevents+Qn
  ID_q<-c(data_UJ$id[which(data_UJ$status==1)],rep(data_UJ$id,each=15))
  
  # t_event_q: updated version when using clock forward timescale
  t_event_q<-vector()
  t_event_q[1:Nevents]<-data_UJ$Tstop[which(data_UJ$status==1)]
  for (i in 1:length(data_UJ$status)) {
    t_event_q<-c(t_event_q,(data_UJ$Tstop[i]-data_UJ$Tstart[i])*xk/2 + (data_UJ$Tstop[i]+data_UJ$Tstart[i])/2)
  }
  
  B_1_event_q <- B_1(t_event_q)
  B_2_event_q <- B_2(t_event_q)
  B_3_event_q <- B_3(t_event_q)
  B_4_event_q <- B_4(t_event_q)
  B_5_event_q <- B_5(t_event_q)
  
  f1_1 <- f1(visits)
  f2_1 <- f2(visits)
  f1_2 <- f1(t_event_q)
  f2_2 <- f2(t_event_q)
  
  qwts <- vector()
  for (i in 1:length(data_UJ$status)) {
    qwts<-c(qwts,(data_UJ$Tstop[i]-data_UJ$Tstart[i])*wk/2)
  }
  
  covdat_cr<-data_UJ$X1
  
  fitLN1 <- stan(file = "JM-ST-H_cf.stan", 
                 data = list(y=y,N=N,n=n,ID=ID,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,qwts=qwts,X=as.matrix(covdat_cr),
                             B_1_event_q=B_1_event_q, B_2_event_q=B_2_event_q, B_3_event_q=B_3_event_q, B_4_event_q=B_4_event_q, B_5_event_q=B_5_event_q, 
                             f1_1=f1_1, f2_1=f2_1, f1_2=f1_2, f2_2=f2_2),        
                 warmup = 400,  # for sample size n=1000 this seems sufficient               
                 iter = 1400,
                 chains = 1, 
                 seed = 2024,
                 init = 0,
                 cores = 1)
  



mean <- c(task_id, summary(fitLN1, pars = c("gamma", "alpha"), probs = c(0.025, 0.975))$summary[, c("mean")])
sd <- c(task_id,summary(fitLN1, pars = c("gamma", "alpha"), probs = c(0.025, 0.975))$summary[, c("sd")])
p95_L <- c(task_id,summary(fitLN1, pars = c("gamma", "alpha"), probs = c(0.025, 0.975))$summary[, c("2.5%")])
p95_U <- c(task_id,summary(fitLN1, pars = c("gamma", "alpha"), probs = c(0.025, 0.975))$summary[, c("97.5%")])
neff <- c(task_id,summary(fitLN1, pars = c("gamma", "alpha"), probs = c(0.025, 0.975))$summary[, c("n_eff")])
Rhat <- c(task_id,summary(fitLN1, pars = c("gamma", "alpha"), probs = c(0.025, 0.975))$summary[, c("Rhat")])
x <- list(mean=mean, sd=sd, p95_L=p95_L, p95_U=p95_U, neff=neff, Rhat=Rhat, runtime=max(rowSums(unname(get_elapsed_time(fitLN1)))))

save(x, file=paste0("JM-ST-C-10_",task_id,".rdata"))
