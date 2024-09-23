###################################################
# JM-MSM for Ferrer's sim model 
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
############################################
# Required quantities for model fitting

y <- obj$longitudinal[,1]          # longitudinal outcomes
ID <- obj$longitudinal[,2]         # patient IDs
visits <- obj$visits               # visit times for repeated observations
N <- length(y)                     # total number of longitudinal outcomes


library(statmod)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights

Nevents<-sum(data_mstate$status==1)
Qn<-length(data_mstate$status)*15
nrow_q<-Nevents+Qn
ID_q<-c(data_mstate$id[which(data_mstate$status==1)],rep(data_mstate$id,each=15))
trans_q<-c(as.numeric(data_mstate$transition)[which(data_mstate$status==1)],rep(as.numeric(data_mstate$transition),each=15))

# t_event_q: updated version when using clock forward timescale
t_event_q<-vector()
t_event_q[1:Nevents]<-data_mstate$Tstop[which(data_mstate$status==1)]
for (i in 1:length(data_mstate$status)) {
  t_event_q<-c(t_event_q,(data_mstate$Tstop[i]-data_mstate$Tstart[i])*xk/2 + (data_mstate$Tstop[i]+data_mstate$Tstart[i])/2)
}

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

B_1_event_q <- B_1(t_event_q)
B_2_event_q <- B_2(t_event_q)
B_3_event_q <- B_3(t_event_q)
B_4_event_q <- B_4(t_event_q)
B_5_event_q <- B_5(t_event_q)
#
f1 <- function (t) {
  return ((1+t)^(-1.2)-1)
}

f2 <- function (t) {return (t)}

f1_1 <- f1(visits)
f2_1 <- f2(visits)
f1_2 <- f1(t_event_q)
f2_2 <- f2(t_event_q)

#
qwts <- vector()
for (i in 1:length(data_mstate$status)) {
  qwts<-c(qwts,(data_mstate$Tstop[i]-data_mstate$Tstart[i])*wk/2)
}

fitLN1 <- stan(file = "JM-MSM_cf.stan", 
               data = list(y=y,N=N,n=n,ID=ID,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,trans_q=trans_q,qwts=qwts,X=as.matrix(covdat),
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

save(x, file=paste0("JM-MSM_",task_id,".rdata"))

pdf(paste0("JM-MSM_tr_",task_id,".pdf"))
traceplot(fitLN1, pars = c("gamma", "alpha"),inc_warmup = T) 
dev.off()

