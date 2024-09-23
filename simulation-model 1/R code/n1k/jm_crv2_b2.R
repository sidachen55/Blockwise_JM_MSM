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

###### fitting Joint competing risk model (later blocks e.g. block 2)
#### version 2: using all historical longitudinal data when fitting the longitudinal model 
## data extraction (after simulating the joint MSM data set)
#
cr_id<-data_mstate$id[which(data_mstate$transition==4)]
length(cr_id)
#
data_cr<-data_mstate[sort(c(which(data_mstate$transition==4),which(data_mstate$transition==5))),]
summary(unique(data_cr$Tstart))
summary(unique(data_cr$years))
table(data_cr$status)

#
L_data<-data.frame(ID=ID,long.data=longit.out,visits=visits.out)
L_data_cr<-L_data[L_data$ID %in% cr_id,]
L_data_cr$Tstop<- rep(unique(data_cr$Tstop),as.vector(table(L_data_cr$ID)))   # different from version 1
L_data_split.by.id <- split(L_data_cr, L_data_cr$ID)
cr_extract <- function (x) {
  x_new<-x[x$visits<x$Tstop,]
}
L_data_split.by.id <- lapply(L_data_split.by.id, cr_extract)
data_lg <- do.call(rbind, L_data_split.by.id)

cr_id2<-unique(data_lg$ID) # excluding subjects in data_cr that have no longitudinal measurements in cr block
length(cr_id2)==length(cr_id) # should be T

# relabel
as.vector(table(data_lg$ID))
data_lg$ID<-rep(1:length(cr_id2),as.vector(table(data_lg$ID)))

data_cr$id<-rep(1:length(cr_id2),each=2)
data_cr$transition<-as.numeric(data_cr$transition)
data_cr$transition<-replace(data_cr$transition,data_cr$transition==c(4,5),c(1,2))

# check
table(data_lg$ID)

############################################
# Required quantities for model fitting

y <- data_lg$long.data          # longitudinal outcomes
ID <- data_lg$ID                # patient IDs
nid <- length(unique(ID))          # number of patients
visits <- data_lg$visits         # visit times for repeated observations
N <- length(y)                     # total number of longitudinal outcomes
n<-length(cr_id2) # total number of subjects involved in the competing risk block

library(statmod)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights

Nevents<-sum(data_cr$status==1)
Qn<-length(data_cr$status)*15
nrow_q<-Nevents+Qn
ID_q<-c(data_cr$id[which(data_cr$status==1)],rep(data_cr$id,each=15))
trans_q<-c(data_cr$transition[which(data_cr$status==1)],rep(data_cr$transition,each=15))
t_start_q<-c(data_cr$Tstart[which(data_cr$status==1)],rep(data_cr$Tstart,each=15)) # time at entry initial state in B_v

t_event_q<-vector()
t_event_q[1:Nevents]<-data_cr$years[which(data_cr$status==1)]
for (i in 1:length(data_cr$status)) {
  t_event_q<-c(t_event_q,data_cr$years[i]*(1+xk)/2)
}

t_cf_event_q<-t_event_q+t_start_q

qwts<-vector()
for (i in 1:length(data_cr$status)) {
  qwts<-c(qwts,data_cr$years[i]*wk/2)
}

covdat_cr<-covdat[cr_id2,]

nchain<-1
fitLN1 <- stan(file = "stan_JM_cr.stan", 
               data = list(y=y,N=N,n=n,ID=ID,visits=visits,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,t_cf_event_q=t_cf_event_q,trans_q=trans_q,qwts=qwts,X=as.matrix(covdat_cr)),        
               warmup = 300,                 
               iter = 1300,
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
save(x, file=paste0("crv2_b2_",task_id,".rdata"))