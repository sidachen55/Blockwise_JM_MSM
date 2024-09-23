#### JM_MSM simulation model 2
library(simsurv)
library(mvtnorm)
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
nid <- length(unique(ID))          # number of patients
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
t_start_q<-c(data_mstate$Tstart[which(data_mstate$status==1)],rep(data_mstate$Tstart,each=15))

t_event_q<-vector()
t_event_q[1:Nevents]<-data_mstate$years[which(data_mstate$status==1)]
for (i in 1:length(data_mstate$status)) {
  t_event_q<-c(t_event_q,data_mstate$years[i]*(1+xk)/2)
}

t_cf_event_q<-t_event_q+t_start_q

qwts<-vector()
for (i in 1:length(data_mstate$status)) {
  qwts<-c(qwts,data_mstate$years[i]*wk/2)
}


nchain<-1
fitLN1 <- stan(file = "stan_JM_MSM.stan", 
               data = list(y=y,N=N,n=n,ID=ID,visits=visits,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,t_cf_event_q=t_cf_event_q,trans_q=trans_q,qwts=qwts,X=as.matrix(covdat)),        
               warmup = 500,                 
               iter = 1500,
               chains = nchain,
               seed = 2022,
               init = 0,
               cores = 1)

mean <- c(task_id, summary(fitLN1, pars = c("gamma", "beta_tilde","Var_e","alpha","delta","lambda"), probs = c(0.025, 0.975))$summary[, c("mean")],sum(unname(get_elapsed_time(fitLN1))))
sd <- c(task_id,summary(fitLN1, pars = c("gamma", "beta_tilde","Var_e","alpha","delta","lambda"), probs = c(0.025, 0.975))$summary[, c("sd")],sum(unname(get_elapsed_time(fitLN1))))
p95_L <- c(task_id,summary(fitLN1, pars = c("gamma", "beta_tilde","Var_e","alpha","delta","lambda"), probs = c(0.025, 0.975))$summary[, c("2.5%")],sum(unname(get_elapsed_time(fitLN1))))
p95_U <- c(task_id,summary(fitLN1, pars = c("gamma", "beta_tilde","Var_e","alpha","delta","lambda"), probs = c(0.025, 0.975))$summary[, c("97.5%")],sum(unname(get_elapsed_time(fitLN1))))
neff <- c(task_id,summary(fitLN1, pars = c("gamma", "beta_tilde","Var_e","alpha","delta","lambda"), probs = c(0.025, 0.975))$summary[, c("n_eff")],sum(unname(get_elapsed_time(fitLN1))))
Rhat <- c(task_id,summary(fitLN1, pars = c("gamma", "beta_tilde","Var_e","alpha","delta","lambda"), probs = c(0.025, 0.975))$summary[, c("Rhat")],sum(unname(get_elapsed_time(fitLN1))))
x <- list(mean=mean, sd=sd, p95_L=p95_L, p95_U=p95_U, neff=neff, Rhat=Rhat)
save(x, file=paste0("MSM_",task_id,".rdata"))