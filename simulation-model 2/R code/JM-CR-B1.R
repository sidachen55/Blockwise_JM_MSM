###################################################
# JM-CR for Ferrer's sim model 
# block 1 (JM-CR-H and JM-CR-C are the same)
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

################################################################################
# data preparation 
################################################################################
cr_id <- data_mstate$id[which(data_mstate$transition==1)]
data_cr <- data_mstate[sort(c(which(data_mstate$transition==1),which(data_mstate$transition==2),which(data_mstate$transition==3),which(data_mstate$transition==4))),]
L_data <- data.frame(ID=ID,long.data=longit.out,visits=visits.out)
L_data_cr<-L_data[L_data$ID %in% cr_id,]
#L_data_cr$Tstop<- rep(unique(data_cr$Tstop),as.vector(table(L_data_cr$ID))) this only work if different id has different Tstop
L_data_cr$Tstop <- data_cr[data_cr$id %in% L_data_cr$ID, "Tstop"][match(L_data_cr$ID, data_cr$id)] # `match(L_data_cr$ID, data_cr$id)` returns the position of the first occurrence of each `ID` from `L_data_cr` in `data_cr$id`.

L_data_split.by.id <- split(L_data_cr, L_data_cr$ID)
cr_extract <- function (x) {
  x_new<-x[x$visits<x$Tstop,]
}
L_data_split.by.id <- lapply(L_data_split.by.id, cr_extract)
data_lg <- do.call(rbind, L_data_split.by.id)

# relabel
data_lg$ID<-rep(1:length(cr_id),as.vector(table(data_lg$ID)))
data_cr$id<-rep(1:length(cr_id),each=4)
data_cr$transition<-as.numeric(data_cr$transition)
data_cr$transition<-replace(data_cr$transition,data_cr$transition==c(1:4),c(1:4))

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

Nevents<-sum(data_cr$status==1)
Qn<-length(data_cr$status)*15
nrow_q<-Nevents+Qn
ID_q<-c(data_cr$id[which(data_cr$status==1)],rep(data_cr$id,each=15))
trans_q<-c(data_cr$transition[which(data_cr$status==1)],rep(data_cr$transition,each=15))

# t_event_q: updated version when using clock forward timescale
t_event_q<-vector()
t_event_q[1:Nevents]<-data_cr$Tstop[which(data_cr$status==1)]
for (i in 1:length(data_cr$status)) {
  t_event_q<-c(t_event_q,(data_cr$Tstop[i]-data_cr$Tstart[i])*xk/2 + (data_cr$Tstop[i]+data_cr$Tstart[i])/2)
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
for (i in 1:length(data_cr$status)) {
  qwts<-c(qwts,(data_cr$Tstop[i]-data_cr$Tstart[i])*wk/2)
}

covdat_cr<-covdat[cr_id,]

fitLN1 <- stan(file = "JM-CR_B1.stan", 
               data = list(y=y,N=N,n=n,ID=ID,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,trans_q=trans_q,qwts=qwts,X=as.matrix(covdat_cr),
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

save(x, file=paste0("JM-CR_B1_",task_id,".rdata"))


