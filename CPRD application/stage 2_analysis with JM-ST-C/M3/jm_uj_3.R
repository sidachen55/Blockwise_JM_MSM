################################################################################
# stage 1: JM-uj-3
################################################################################
library(bit64)
library(ggplot2)
library(rstan)
library(grDevices)
library(splines)
################
load("MSMdata.RData")
set.seed(1)

#######################################
data_cr<-data_mstate[sort(c(which(data_mstate$transition==1),which(data_mstate$transition==2),which(data_mstate$transition==3))),]


data_UJ<-data_cr[sort(c(which((data_cr$transition==3) & (data_cr$status==1)),which((data_cr$transition==3) & (data_cr$status==0)))),]


#
L_data_cr<-long
L_data_cr$Tstop<- rep(data_UJ$Tstop,as.vector(table(L_data_cr$ID)))
L_data_split.by.id <- split(L_data_cr, L_data_cr$ID)
cr_extract <- function (x) {
  x_new<-x[x$visits<x$Tstop,]
}
L_data_split.by.id <- lapply(L_data_split.by.id, cr_extract)
data_lg <- do.call(rbind, L_data_split.by.id)

cr_id2<-unique(data_lg$ID) 
data_UJ<-data_UJ[data_UJ$id %in% cr_id2,]

# relabel
data_lg$ID<-rep(1:length(cr_id2),as.vector(table(data_lg$ID)))

data_UJ$id<-rep(1:length(cr_id2),each=1)
data_UJ$transition<-as.numeric(data_UJ$transition)
data_UJ$transition<-replace(data_UJ$transition,data_UJ$transition==c(3),c(1))

############################################
# Required quantities for model fitting

y <- data_lg$long.data          # longitudinal outcomes
ID <- data_lg$ID                # patient IDs
visits <- data_lg$visits         # visit times for repeated observations
N <- length(y)                     # total number of longitudinal outcomes
n<-length(cr_id2) # total number of subjects involved in the competing risk block

id_lg_min <- match(1:n,ID)
id_lg_n <- as.vector(table(ID))


library(statmod)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights

Nevents<-sum(data_UJ$status==1)
Qn<-length(data_UJ$status)*15
nrow_q<-Nevents+Qn
ID_q<-c(data_UJ$id[which(data_UJ$status==1)],rep(data_UJ$id,each=15))

t_start_q<-rep(0,nrow_q) # clock-reset, Tstart =0
X_q <- c(which(data_UJ$status==1),rep(1:length(data_UJ$id), each=15)) # for cov
t_event_q<-vector()
t_event_q[1:Nevents]<-data_UJ$years[which(data_UJ$status==1)]
for (i in 1:length(data_UJ$status)) {
  t_event_q<-c(t_event_q,data_UJ$years[i]*(1+xk)/2)
}

t_cf_event_q<-t_event_q+t_start_q

T_e <- data_UJ$years
T_status <- data_UJ$status
#################################################################################
# construct B-spline basis for modelling the baseline transition
# 5 basis functions (1 internal knot located at median observed transition time)
#################################################################################
obs_median_transtime <- summary(data_mstate[data_mstate$transition == 3 & data_mstate$status == 1, ]$years)["Median"]

T_max <- as.numeric(ceiling(summary(data_mstate$years)["Max."]))
inknots <- obs_median_transtime

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
# 
B_1_event_q <- B_1(t_event_q)
B_2_event_q <- B_2(t_event_q)
B_3_event_q <- B_3(t_event_q)
B_4_event_q <- B_4(t_event_q)
B_5_event_q <- B_5(t_event_q)

#
B_1_gq_q <- B_2_gq_q <- B_3_gq_q <- B_4_gq_q <- B_5_gq_q <- matrix(NA,n,15)
for (i in 1:n) {
  for (j in 1:15) {
    B_1_gq_q[i,j] <- B_1(T_e[i]*(1+xk[j])/2)
    B_2_gq_q[i,j] <- B_2(T_e[i]*(1+xk[j])/2)
    B_3_gq_q[i,j] <- B_3(T_e[i]*(1+xk[j])/2)
    B_4_gq_q[i,j] <- B_4(T_e[i]*(1+xk[j])/2)
    B_5_gq_q[i,j] <- B_5(T_e[i]*(1+xk[j])/2)
  }
}
#
B_1_gq_haz <- B_2_gq_haz <- B_3_gq_haz <- B_4_gq_haz <- B_5_gq_haz <- rep(NA,n)
for (i in 1:n) {
  B_1_gq_haz[i] <- B_1(T_e[i])
  B_2_gq_haz[i] <- B_2(T_e[i])
  B_3_gq_haz[i] <- B_3(T_e[i])
  B_4_gq_haz[i] <- B_4(T_e[i])
  B_5_gq_haz[i] <- B_5(T_e[i])
}
################################################################################


qwts<-vector()
for (i in 1:length(data_UJ$status)) {
  qwts<-c(qwts,data_UJ$years[i]*wk/2)
}

covdat <- cbind(data_UJ$X1, data_UJ$X2, data_UJ$X3, data_UJ$X4) 

nchain<-1
fitLN1 <- stan(file = "stan_JM_ST_as3.stan", 
               data = list(y=y,N=N,n=n,ID=ID,T_e=T_e, T_status=T_status,id_lg_min=id_lg_min, id_lg_n=id_lg_n, xk=xk, wk=wk,visits=visits,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,t_cf_event_q=t_cf_event_q,qwts=qwts,X=covdat,
                           B_1_event_q=B_1_event_q, B_2_event_q=B_2_event_q, B_3_event_q=B_3_event_q, B_4_event_q=B_4_event_q, B_5_event_q=B_5_event_q, X_q=X_q, nrow_msmdata=length(data_UJ$status),
                           B_1_gq_q=B_1_gq_q, B_2_gq_q=B_2_gq_q, B_3_gq_q=B_3_gq_q, B_4_gq_q=B_4_gq_q, B_5_gq_q=B_5_gq_q, B_1_gq_haz=B_1_gq_haz, B_2_gq_haz=B_2_gq_haz, B_3_gq_haz=B_3_gq_haz, B_4_gq_haz=B_4_gq_haz, B_5_gq_haz=B_5_gq_haz),
               warmup = 500,                 
               iter = 1500,
               chains = nchain,
               seed = 2024,
               init = 0,
               cores = 1)

mean <- c(summary(fitLN1, pars = c("gamma1","gamma2","gamma3","gamma4", "beta_tilde","Var_e","alpha","rho"), probs = c(0.025, 0.975))$summary[, c("mean")],sum(unname(get_elapsed_time(fitLN1))))
sd <- c(summary(fitLN1, pars = c("gamma1","gamma2","gamma3","gamma4", "beta_tilde","Var_e","alpha","rho"), probs = c(0.025, 0.975))$summary[, c("sd")],sum(unname(get_elapsed_time(fitLN1))))
p95_L <- c(summary(fitLN1, pars = c("gamma1","gamma2","gamma3","gamma4", "beta_tilde","Var_e","alpha","rho"), probs = c(0.025, 0.975))$summary[, c("2.5%")],sum(unname(get_elapsed_time(fitLN1))))
p95_U <- c(summary(fitLN1, pars = c("gamma1","gamma2","gamma3","gamma4", "beta_tilde","Var_e","alpha","rho"), probs = c(0.025, 0.975))$summary[, c("97.5%")],sum(unname(get_elapsed_time(fitLN1))))
neff <- c(summary(fitLN1, pars = c("gamma1","gamma2","gamma3","gamma4", "beta_tilde","Var_e","alpha","rho"), probs = c(0.025, 0.975))$summary[, c("n_eff")],sum(unname(get_elapsed_time(fitLN1))))
Rhat <- c(summary(fitLN1, pars = c("gamma1","gamma2","gamma3","gamma4", "beta_tilde","Var_e","alpha","rho"), probs = c(0.025, 0.975))$summary[, c("Rhat")],sum(unname(get_elapsed_time(fitLN1))))

detach("package:rstan", unload = TRUE)
library(loo)
log_lik_1 <- extract_log_lik(fitLN1)
xraw <- list(mean=mean, sd=sd, p95_L=p95_L, p95_U=p95_U, neff=neff, Rhat=Rhat,log_lik=log_lik_1)

save(xraw, file="uj_3_fitobj.rdata")

pdf("trace_uj3.pdf")
traceplot(fitLN1, pars = c("alpha","gamma1","gamma2","gamma3","gamma4","rho","Var_e"),inc_warmup = T) 
dev.off()
