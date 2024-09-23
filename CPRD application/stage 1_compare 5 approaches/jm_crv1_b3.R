################################################################################
# stage 1: JM-crv1b3
################################################################################
library(bit64)
library(ggplot2)
library(rstan)
library(grDevices)
library(splines)
################

load("MSMdata.RData")
set.seed(1)

cr_id<-data_mstate$id[which(data_mstate$transition==6)]
data_cr<-data_mstate[sort(c(which(data_mstate$transition==6),which(data_mstate$transition==7))),]

L_data_cr<-long[long$ID %in% cr_id,]
L_data_cr$Tstart<- rep(data_cr[!duplicated(data_cr$id), "Tstart"],as.vector(table(L_data_cr$ID)))
L_data_cr$Tstop<- rep(data_cr[!duplicated(data_cr$id), "Tstop"],as.vector(table(L_data_cr$ID)))
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

cr_id2<-unique(data_lg$ID) 
data_cr<-data_cr[data_cr$id %in% cr_id2,]

data_lg$ID<-rep(1:length(cr_id2),as.vector(table(data_lg$ID)))
data_cr$id<-rep(1:length(cr_id2),each=2)
data_cr$transition<-as.numeric(data_cr$transition)
data_cr$transition<-replace(data_cr$transition,data_cr$transition==c(6,7),c(1,2))

############################################
# Required quantities for model fitting

y <- data_lg$long.data          # longitudinal outcomes
ID <- data_lg$ID                # patient IDs
visits <- data_lg$visits         
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
t_start_q<-c(data_cr$Tstart[which(data_cr$status==1)],rep(data_cr$Tstart,each=15))

X_q <- c(which(data_cr$status==1),rep(1:length(data_cr$id), each=15)) # for cov


t_event_q<-vector()
t_event_q[1:Nevents]<-data_cr$years[which(data_cr$status==1)]
for (i in 1:length(data_cr$status)) {
  t_event_q<-c(t_event_q,data_cr$years[i]*(1+xk)/2)
}

t_cf_event_q<-t_event_q+t_start_q

#################################################################################
# construct B-spline basis for modelling the baseline transition
# 5 basis functions (1 internal knot located at median observed transition time)
#################################################################################
summary(data_mstate$years)
obs_median_transtime <- rep(NA,2)
trans_id <- c(6,7)
for (i in 1:2) {obs_median_transtime[i] <- summary(data_mstate[data_mstate$transition == trans_id[i] & data_mstate$status == 1, ]$years)["Median"]}

T_max <- as.numeric(ceiling(summary(data_mstate$years)["Max."]))
inknots <- obs_median_transtime

B_1 <- function(t, trans) {
  result <- spline.des(knots = c(rep(0,4), inknots[trans], rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,1]
  return(result)
}
B_2 <- function(t, trans) {
  result <- spline.des(knots = c(rep(0,4), inknots[trans], rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,2]
  return(result)
}
B_3 <- function(t, trans) {
  result <- spline.des(knots = c(rep(0,4), inknots[trans], rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,3]
  return(result)
}
B_4 <- function(t, trans) {
  result <- spline.des(knots = c(rep(0,4), inknots[trans], rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,4]
  return(result)
}
B_5 <- function(t, trans) {
  result <- spline.des(knots = c(rep(0,4), inknots[trans], rep(T_max,4)), x = t, ord = 4, outer.ok = TRUE)$design[,5]
  return(result)
}
# Vectorizing 
B_1_vectorized <- Vectorize(B_1, vectorize.args = c("t", "trans"))
B_2_vectorized <- Vectorize(B_2, vectorize.args = c("t", "trans"))
B_3_vectorized <- Vectorize(B_3, vectorize.args = c("t", "trans"))
B_4_vectorized <- Vectorize(B_4, vectorize.args = c("t", "trans"))
B_5_vectorized <- Vectorize(B_5, vectorize.args = c("t", "trans"))


B_1_event_q <- B_1_vectorized(t_event_q, trans_q)
B_2_event_q <- B_2_vectorized(t_event_q, trans_q)
B_3_event_q <- B_3_vectorized(t_event_q, trans_q)
B_4_event_q <- B_4_vectorized(t_event_q, trans_q)
B_5_event_q <- B_5_vectorized(t_event_q, trans_q)
################################################################################

qwts<-vector()
for (i in 1:length(data_cr$status)) {
  qwts<-c(qwts,data_cr$years[i]*wk/2)
}

covdat <- cbind(data_cr$X1, data_cr$X2, data_cr$X3, data_cr$X4) 

nchain<-1
fitLN1 <- stan(file = "stan_JM_cr.stan", 
               data = list(y=y,N=N,n=n,ID=ID,visits=visits,Nevents=Nevents,Qn=Qn,nrow_q=nrow_q,ID_q=ID_q,t_event_q=t_event_q,t_cf_event_q=t_cf_event_q,trans_q=trans_q,qwts=qwts,X=covdat,
                           B_1_event_q=B_1_event_q, B_2_event_q=B_2_event_q, B_3_event_q=B_3_event_q, B_4_event_q=B_4_event_q, B_5_event_q=B_5_event_q, X_q=X_q, nrow_msmdata=length(data_cr$status)),
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
x <- list(mean=mean, sd=sd, p95_L=p95_L, p95_U=p95_U, neff=neff, Rhat=Rhat)
save(x, file=paste0("crv1b3.rdata"))

pdf("trace_crv1b3.pdf")
traceplot(fitLN1, pars = c("alpha","gamma1","gamma2","gamma3","gamma4","rho","Var_e"),inc_warmup = T) 
dev.off()