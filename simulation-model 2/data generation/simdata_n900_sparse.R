################################################################################
# new simulation model - motivated from Ferrer et al (2016) 
################################################################################
library(simsurv)
library(mvtnorm)
library(splines)
library(MASS)
n<-900
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

set.seed(task_id)

################################################################################
# useful functions 
################################################################################
# longitudinal submodel:
f1 <- function (t) {
  return ((1+t)^(-1.2)-1)
}

f2 <- function (t) {return (t)}

mu_long <- function (t, x, betas) {
  return(betas[["beta_0i"]]+betas[["beta_0x"]]*x[["w1"]]+(betas[["beta_1i"]]+betas[["beta_1x"]]*x[["w1"]])*f1(t) + (betas[["beta_2i"]]+betas[["beta_2x"]]*x[["w1"]])*f2(t))
}

#
T_max <- 24
inknots <- c(7.127)
theta_B_surv_1 <- c(-6.990, -3.428, -5.578, -6.847, -6.067)
theta_B_surv_2 <- c(-7.537, -5.525, -7.191, -4.692, -6.420)
theta_B_surv_3 <- c(-7.322, -4.345, -3.770, -1.496, -1.261)

# baseline hazards
#Bs_surv <- function(t, theta_B) {
#  result <- ifelse(t < 0.026 | t > 18.201, 0, 
#                   exp(as.vector(spline.des(knots = c(rep(0.026,4), inknots, rep(18.201,4)), 
#                                            x = t, ord = 4, outer.ok = TRUE)$design %*% theta_B)))
#  return(result)
# }

Bs_surv <- function(t, theta_B) {
  result <- exp(as.vector(spline.des(knots = c(rep(0,4), inknots, rep(T_max,4)), 
                                     x = t, ord = 4, outer.ok = TRUE)$design %*% theta_B))
  return(result)
}

# transition specific hazard functions [clock-forward]  
haz_01 <- haz_03 <- haz_04 <- function(t, x, betas, ...) {
  Bs_surv(t, c(betas[["B_1"]],betas[["B_2"]],betas[["B_3"]],betas[["B_4"]],betas[["B_5"]]))*exp(betas[["gamma_1"]]*x[["w1"]] + betas[["alpha"]]*mu_long(t,x,betas))
}

haz_02 <- haz_12 <- function(t, x, betas) {
  exp(betas[["zeta"]])*haz_01(t, x, betas)
}


haz_13 <- haz_23 <- function(t, x, betas) {
  exp(betas[["zeta"]])*haz_03(t, x, betas)
}

haz_14 <- haz_24 <- haz_34 <- function(t, x, betas) {
  exp(betas[["zeta"]])*haz_04(t, x, betas)
}

# transition specific cumulative hazard functions [cf]
cum_haz_01 <- function (t, u, x, betas) {
  integrate(haz_01, lower = 0, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_02 <- function (t, u, x, betas) {
  integrate(haz_02, lower = 0, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_03 <- function (t, u, x, betas) {
  integrate(haz_03, lower = 0, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_04 <- function (t, u, x, betas) {
  integrate(haz_04, lower = 0, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_12 <- function (t, left_trunc, u, x, betas) {
  integrate(haz_12, lower = left_trunc, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_13 <- function (t, left_trunc, u, x, betas) {
  integrate(haz_13, lower = left_trunc, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_14 <- function (t, left_trunc, u, x, betas) {
  integrate(haz_14, lower = left_trunc, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_23 <- function (t, left_trunc, u, x, betas) {
  integrate(haz_23, lower = left_trunc, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_24 <- function (t, left_trunc, u, x, betas) {
  integrate(haz_24, lower = left_trunc, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}

cum_haz_34 <- function (t, left_trunc, u, x, betas) {
  integrate(haz_34, lower = left_trunc, upper = t, x, betas, subdivisions = 10000L)$value + log(u)
}


################################################################################
# step 1: generate random effects
################################################################################
Sigma_b <- matrix(c(0.37, 0.35, 0.011,
                    0.35, 1.7, 0.31,
                    0.011, 0.31, 0.17), nrow=3, byrow=TRUE)
b <- mvrnorm(n = n, mu = rep(0,3), Sigma = Sigma_b)

################################################################################
# step 2: simulate multistate data
################################################################################
# 
beta_0 <- -0.26
beta_1 <- 0.95
beta_2 <- -0.09
beta_0i<-rep(beta_0, n) + b[, 1] 
beta_1i<-rep(beta_1, n) + b[, 2] 
beta_2i<-rep(beta_2, n) + b[, 3]


betas_01 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_1[1], n),  B_2 = rep(theta_B_surv_1[2], n),  B_3 = rep(theta_B_surv_1[3], n),  B_4 = rep(theta_B_surv_1[4], n),  B_5 = rep(theta_B_surv_1[5], n),
                       gamma_1 = rep(0.047, n), alpha = rep(0.36, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_03 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_2[1], n),  B_2 = rep(theta_B_surv_2[2], n),  B_3 = rep(theta_B_surv_2[3], n),  B_4 = rep(theta_B_surv_2[4], n),  B_5 = rep(theta_B_surv_2[5], n),
                       gamma_1 = rep(0.05, n), alpha = rep(0.4, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_04 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_3[1], n),  B_2 = rep(theta_B_surv_3[2], n),  B_3 = rep(theta_B_surv_3[3], n),  B_4 = rep(theta_B_surv_3[4], n),  B_5 = rep(theta_B_surv_3[5], n),
                       gamma_1 = rep(0.06, n), alpha = rep(-0.15, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_02 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_1[1], n),  B_2 = rep(theta_B_surv_1[2], n),  B_3 = rep(theta_B_surv_1[3], n),  B_4 = rep(theta_B_surv_1[4], n),  B_5 = rep(theta_B_surv_1[5], n),
                       zeta= rep(-0.9,n), gamma_1 = rep(0.489, n), alpha = rep(0.5, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_12 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_1[1], n),  B_2 = rep(theta_B_surv_1[2], n),  B_3 = rep(theta_B_surv_1[3], n),  B_4 = rep(theta_B_surv_1[4], n),  B_5 = rep(theta_B_surv_1[5], n),
                       zeta= rep(2,n), gamma_1 = rep(0.916, n), alpha = rep(-0.2, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_13 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_2[1], n),  B_2 = rep(theta_B_surv_2[2], n),  B_3 = rep(theta_B_surv_2[3], n),  B_4 = rep(theta_B_surv_2[4], n),  B_5 = rep(theta_B_surv_2[5], n),
                       zeta= rep(4,n), gamma_1 = rep(0.351, n), alpha = rep(-0.43, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_23 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_2[1], n),  B_2 = rep(theta_B_surv_2[2], n),  B_3 = rep(theta_B_surv_2[3], n),  B_4 = rep(theta_B_surv_2[4], n),  B_5 = rep(theta_B_surv_2[5], n),
                       zeta= rep(3,n), gamma_1 = rep(-0.21, n), alpha = rep(-0.17, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_14 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_3[1], n),  B_2 = rep(theta_B_surv_3[2], n),  B_3 = rep(theta_B_surv_3[3], n),  B_4 = rep(theta_B_surv_3[4], n),  B_5 = rep(theta_B_surv_3[5], n),
                       zeta= rep(2,n), gamma_1 = rep(-0.385, n), alpha = rep(0.1, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_24 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_3[1], n),  B_2 = rep(theta_B_surv_3[2], n),  B_3 = rep(theta_B_surv_3[3], n),  B_4 = rep(theta_B_surv_3[4], n),  B_5 = rep(theta_B_surv_3[5], n),
                       zeta= rep(0.5,n), gamma_1 = rep(0.007, n), alpha = rep(0.05, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)

betas_34 <- data.frame(beta_0x = rep(0.799, n), beta_1x = rep(0.905, n), beta_2x = rep(0.207, n),
                       B_1 = rep(theta_B_surv_3[1], n),  B_2 = rep(theta_B_surv_3[2], n),  B_3 = rep(theta_B_surv_3[3], n),  B_4 = rep(theta_B_surv_3[4], n),  B_5 = rep(theta_B_surv_3[5], n),
                       zeta= rep(1.7,n), gamma_1 = rep(0.262, n), alpha = rep(0.05, n),
                       beta_0i = beta_0i, beta_1i = beta_1i, beta_2i = beta_2i)


# baseline covariates
w1 <- rnorm(n, mean = 2.04, sd = sqrt(0.5))
covdat <- data.frame(w1 = w1)


# initiate vectors to save true event times (in a clock-forward framework)
trueT01 <- trueT02 <- trueT03 <- trueT04 <- trueT12 <- trueT13 <- trueT14 <- trueT23 <- trueT24 <- trueT34 <- numeric(n)

# generate u[0,1] r.v. to be used for the inversion method
u01 <- runif(n, 0, 1)
u02 <- runif(n, 0, 1)
u03 <- runif(n, 0, 1)
u04 <- runif(n, 0, 1)
u12 <- runif(n, 0, 1)
u13 <- runif(n, 0, 1)
u14 <- runif(n, 0, 1)
u23 <- runif(n, 0, 1)
u24 <- runif(n, 0, 1)
u34 <- runif(n, 0, 1)
# sample censoring times
Ctimes <- runif(n, 4, T_max)

Up <- 200
# simulate time-to-event data (for each transition)
for (i in 1:n) {
  Root01 <- Root02 <- Root03 <- Root04 <- Root12 <- Root13 <- Root14 <- Root23 <- Root24 <- Root34 <-  NULL
  # Transition 0->1
  Root01 <- try(uniroot(cum_haz_01, interval = c(1e-05, Up), u = u01[i], x = covdat[i,,drop=FALSE], betas = betas_01[i,])$root, TRUE)
  trueT01[i] <- if (!inherits(Root01, "try-error")) Root01 else 500
  # Transition 0->2
  Root02 <- try(uniroot(cum_haz_02, interval = c(1e-05, Up), u = u02[i], x = covdat[i,,drop=FALSE], betas = betas_02[i,])$root, TRUE)
  trueT02[i] <- if (!inherits(Root02, "try-error")) Root02 else 500
  # Transition 0->3
  Root03 <- try(uniroot(cum_haz_03, interval = c(1e-05, Up), u = u03[i], x = covdat[i,,drop=FALSE], betas = betas_03[i,])$root, TRUE)
  trueT03[i] <- if (!inherits(Root03, "try-error")) Root03 else 500
  # Transition 0->4
  Root04 <- try(uniroot(cum_haz_04, interval = c(1e-05, Up), u = u04[i], x = covdat[i,,drop=FALSE], betas = betas_04[i,])$root, TRUE)
  trueT04[i] <- if (!inherits(Root04, "try-error")) Root04 else 500
  
  # censored within block 1
  if (Ctimes[i]<min(trueT01[i], trueT02[i], trueT03[i], trueT04[i])) {
    trueT12[i]<-500; trueT13[i]<-500; trueT14[i]<-500; trueT23[i]<-500; trueT24[i]<-500; trueT34[i]<-500
  } else if (trueT04[i]<min(trueT01[i], trueT02[i], trueT03[i])) {
    trueT12[i]<-500; trueT13[i]<-500; trueT14[i]<-500; trueT23[i]<-500; trueT24[i]<-500; trueT34[i]<-500
  } else {
    if (trueT01[i]<min(trueT02[i], trueT03[i])) {
      Root12 <- try(uniroot(cum_haz_12, interval = c(trueT01[i], Up), left_trunc=trueT01[i], u = u12[i], x = covdat[i,,drop=FALSE], betas = betas_12[i,])$root, TRUE)
      trueT12[i] <- if (!inherits(Root12, "try-error")) {Root12-trueT01[i]} else {500}
      Root13 <- try(uniroot(cum_haz_13, interval = c(trueT01[i], Up), left_trunc=trueT01[i], u = u13[i], x = covdat[i,,drop=FALSE], betas = betas_13[i,])$root, TRUE)
      trueT13[i] <- if (!inherits(Root13, "try-error")) {Root13-trueT01[i]} else {500}
      Root14 <- try(uniroot(cum_haz_14, interval = c(trueT01[i], Up), left_trunc=trueT01[i], u = u14[i], x = covdat[i,,drop=FALSE], betas = betas_14[i,])$root, TRUE)
      trueT14[i] <- if (!inherits(Root14, "try-error")) {Root14-trueT01[i]} else {500}
      
      if ((trueT12[i]<min(trueT13[i],trueT14[i])) & (trueT01[i]+trueT12[i]<Ctimes[i])) {
        Root23 <- try(uniroot(cum_haz_23, interval = c(trueT01[i]+trueT12[i], Up), left_trunc=trueT01[i]+trueT12[i], u = u23[i], x = covdat[i,,drop=FALSE], betas = betas_23[i,])$root, TRUE)
        trueT23[i] <- if (!inherits(Root23, "try-error")) {Root23-(trueT01[i]+trueT12[i])} else {500}
        Root24 <- try(uniroot(cum_haz_24, interval = c(trueT01[i]+trueT12[i], Up), left_trunc=trueT01[i]+trueT12[i], u = u24[i], x = covdat[i,,drop=FALSE], betas = betas_24[i,])$root, TRUE)
        trueT24[i] <- if (!inherits(Root24, "try-error")) {Root24-(trueT01[i]+trueT12[i])} else {500}
        
        if ((trueT23[i]<trueT24[i]) & (trueT01[i]+trueT12[i]+trueT23[i]<Ctimes[i])) {
          Root34 <- try(uniroot(cum_haz_34, interval = c(trueT01[i]+trueT12[i]+trueT23[i], Up), left_trunc=trueT01[i]+trueT12[i]+trueT23[i], u = u34[i], x = covdat[i,,drop=FALSE], betas = betas_34[i,])$root, TRUE)
          trueT34[i] <- if (!inherits(Root34, "try-error")) {Root34-(trueT01[i]+trueT12[i]+trueT23[i])} else {500}
          
        } else {
          trueT34[i]<-500
        }
      } else if ((trueT13[i]<min(trueT12[i],trueT14[i])) & (trueT01[i]+trueT13[i]<Ctimes[i])) {
        Root34 <- try(uniroot(cum_haz_34, interval = c(trueT01[i]+trueT13[i], Up), left_trunc=trueT01[i]+trueT13[i], u = u34[i], x = covdat[i,,drop=FALSE], betas = betas_34[i,])$root, TRUE)
        trueT34[i] <- if (!inherits(Root34, "try-error")) {Root34-(trueT01[i]+trueT13[i])} else {500}
        trueT23[i]<-500; trueT24[i]<-500
      } else {
        trueT23[i]<-500; trueT24[i]<-500; trueT34[i]<-500
      }
    } else if (trueT02[i]<min(trueT01[i], trueT03[i])) {
      Root23 <- try(uniroot(cum_haz_23, interval = c(trueT02[i], Up), left_trunc=trueT02[i], u = u23[i], x = covdat[i,,drop=FALSE], betas = betas_23[i,])$root, TRUE)
      trueT23[i] <- if (!inherits(Root23, "try-error")) {Root23-trueT02[i]} else {500}
      Root24 <- try(uniroot(cum_haz_24, interval = c(trueT02[i], Up), left_trunc=trueT02[i], u = u24[i], x = covdat[i,,drop=FALSE], betas = betas_24[i,])$root, TRUE)
      trueT24[i] <- if (!inherits(Root24, "try-error")) {Root24-trueT02[i]} else {500}
      if ((trueT23[i]<trueT24[i]) & (trueT02[i]+trueT23[i]<Ctimes[i])) {
        Root34 <- try(uniroot(cum_haz_34, interval = c(trueT02[i]+trueT23[i], Up), left_trunc=trueT02[i]+trueT23[i], u = u34[i], x = covdat[i,,drop=FALSE], betas = betas_34[i,])$root, TRUE)
        trueT34[i] <- if (!inherits(Root34, "try-error")) {Root34-(trueT02[i]+trueT23[i])} else {500}
        trueT13[i]<-500; trueT14[i]<-500
      } else {
        trueT13[i]<-500; trueT14[i]<-500; trueT34[i]<-500
      }
    } else if (trueT03[i]<min(trueT01[i], trueT02[i])) {
      Root34 <- try(uniroot(cum_haz_34, interval = c(trueT03[i], Up), left_trunc=trueT03[i], u = u34[i], x = covdat[i,,drop=FALSE], betas = betas_34[i,])$root, TRUE)
      trueT34[i] <- if (!inherits(Root34, "try-error")) {Root34-trueT03[i]} else {500}
      trueT12[i]<-500; trueT13[i]<-500; trueT14[i]<-500; trueT23[i]<-500; trueT24[i]<-500
    }
  }
  print(i)
}

data_mstate <- data.frame('id' = 1:n, 'trueT01' = trueT01, 'trueT02' = trueT02, 'trueT03' = trueT03, 'trueT04' = trueT04, 'trueT12' = trueT12, 'trueT13' = trueT13, 'trueT14' = trueT14,'trueT23' = trueT23, 'trueT24' = trueT24, 'trueT34' = trueT34,
                          'Ctimes' = Ctimes, 'X1' = covdat$w1)
data_mstate_split.by.id <- split(data_mstate, data_mstate$id)

ms_arrange <- function (x) {
  if (x$Ctimes < min(x$trueT01, x$trueT02, x$trueT03, x$trueT04)) {
    x_new <- data.frame('id' = rep(x$id, 4), 'from_state' = rep(0, 4), 'to_state' = 1:4, 
                        'transition' = 1:4, 'Tstart' = rep(0, 4), 'Tstop' = x$Ctimes, 'status' = rep(0, 4), 
                        'X1' = x$X1)
  } else if (x$trueT04 < min(x$trueT01, x$trueT02, x$trueT03)) {
    x_new <- data.frame('id' = rep(x$id, 4), 'from_state' = rep(0, 4), 'to_state' = 1:4, 
                        'transition' = 1:4, 'Tstart' = rep(0, 4), 'Tstop' = rep(x$trueT04,4), 'status' = c(0,0,0,1), 
                        'X1' = x$X1)
  } else {
    if (x$trueT01 < min(x$trueT02, x$trueT03)) {
      if(x$Ctimes < min(x$trueT01+x$trueT12,x$trueT01+x$trueT13,x$trueT01+x$trueT14)){
        x_new <- data.frame('id' = rep(x$id, 7), 'from_state' = c(0, 0, 0, 0, 1,1,1), 'to_state' = c(1,2,3,4,2,3,4), 
                            'transition' = 1:7, 'Tstart' = c(rep(0, 4), rep(x$trueT01,3)), 
                            'Tstop' = c(rep(x$trueT01, 4), rep(x$Ctimes,3)), 'status' = c(1, 0, 0,0,0,0,0), 
                            'X1' = x$X1)
      } else if (x$trueT14<min(x$trueT12, x$trueT13)) {
        x_new <- data.frame('id' = rep(x$id, 7), 'from_state' = c(0, 0, 0, 0, 1,1,1), 'to_state' = c(1,2,3,4,2,3,4), 
                            'transition' = 1:7, 'Tstart' = c(rep(0, 4), rep(x$trueT01,3)), 
                            'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT14,3)), 'status' = c(1, 0, 0,0,0,0,1), 
                            'X1' = x$X1)
      } else if (x$trueT12<min(x$trueT14, x$trueT13)) {
        if (x$Ctimes < min(x$trueT01+x$trueT12+x$trueT23, x$trueT01+x$trueT12+x$trueT24)) {
          x_new <- data.frame('id' = rep(x$id, 9), 'from_state' = c(0, 0, 0, 0, 1,1,1,2,2), 'to_state' = c(1,2,3,4,2,3,4,3,4), 
                              'transition' = 1:9, 'Tstart' = c(rep(0, 4), rep(x$trueT01,3), rep(x$trueT01+x$trueT12,2)), 
                              'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT12,3), rep(x$Ctimes,2)), 'status' = c(1, 0, 0, 0, 1,0,0,0,0), 
                              'X1' = x$X1)
        } else if (x$trueT24 < x$trueT23) {
          x_new <- data.frame('id' = rep(x$id, 9), 'from_state' = c(0, 0, 0, 0, 1,1,1,2,2), 'to_state' = c(1,2,3,4,2,3,4,3,4), 
                              'transition' = 1:9, 'Tstart' = c(rep(0, 4), rep(x$trueT01,3), rep(x$trueT01+x$trueT12,2)), 
                              'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT12,3), rep(x$trueT01+x$trueT12+x$trueT24,2)), 'status' = c(1, 0, 0, 0, 1,0,0,0,1), 
                              'X1' = x$X1)
        } else {
          if (x$Ctimes < x$trueT01+x$trueT12+x$trueT23+x$trueT34) {
            x_new <- data.frame('id' = rep(x$id, 10), 'from_state' = c(0, 0, 0, 0, 1,1,1,2,2,3), 'to_state' = c(1,2,3,4,2,3,4,3,4,4), 
                                'transition' = 1:10, 'Tstart' = c(rep(0, 4), rep(x$trueT01,3), rep(x$trueT01+x$trueT12,2),x$trueT01+x$trueT12+x$trueT23), 
                                'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT12,3), rep(x$trueT01+x$trueT12+x$trueT23,2),x$Ctimes), 'status' = c(1, 0, 0, 0, 1,0,0,1,0,0), 
                                'X1' = x$X1)
          } else {
            x_new <- data.frame('id' = rep(x$id, 10), 'from_state' = c(0, 0, 0, 0, 1,1,1,2,2,3), 'to_state' = c(1,2,3,4,2,3,4,3,4,4), 
                                'transition' = 1:10, 'Tstart' = c(rep(0, 4), rep(x$trueT01,3), rep(x$trueT01+x$trueT12,2),x$trueT01+x$trueT12+x$trueT23), 
                                'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT12,3), rep(x$trueT01+x$trueT12+x$trueT23,2),x$trueT01+x$trueT12+x$trueT23+x$trueT34), 'status' = c(1, 0, 0, 0, 1,0,0,1,0,1), 
                                'X1' = x$X1)
          }
        }
      } else {
        if (x$Ctimes < x$trueT01+x$trueT13+x$trueT34) {
          x_new <- data.frame('id' = rep(x$id, 8), 'from_state' = c(0, 0, 0,0, 1,1,1,3), 'to_state' = c(1,2,3,4,2,3,4,4), 
                              'transition' = c(1:7,10), 'Tstart' = c(rep(0, 4), rep(x$trueT01,3), x$trueT01+x$trueT13), 
                              'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT13,3), x$Ctimes), 'status' = c(1, 0, 0,0,0, 1,0,0), 
                              'X1' = x$X1)
        } else {
          x_new <- data.frame('id' = rep(x$id, 8), 'from_state' = c(0, 0, 0,0, 1,1,1,3), 'to_state' = c(1,2,3,4,2,3,4,4), 
                              'transition' = c(1:7,10), 'Tstart' = c(rep(0, 4), rep(x$trueT01,3), x$trueT01+x$trueT13), 
                              'Tstop' = c(rep(x$trueT01, 4), rep(x$trueT01+x$trueT13,3), x$trueT01+x$trueT13+x$trueT34), 'status' = c(1, 0, 0,0,0, 1,0,1), 
                              'X1' = x$X1)
        }
      }
    } else if (x$trueT02 < min(x$trueT01, x$trueT03)) {
      if (x$Ctimes<min(x$trueT02+x$trueT23,x$trueT02+x$trueT24)) {
        x_new <- data.frame('id' = rep(x$id, 6), 'from_state' = c(0, 0, 0, 0, 2,2), 'to_state' = c(1,2,3,4,3,4), 
                            'transition' = c(1,2,3,4,8,9), 'Tstart' = c(rep(0, 4), rep(x$trueT02,2)), 
                            'Tstop' = c(rep(x$trueT02, 4), rep(x$Ctimes,2)), 'status' = c(0, 1,0, 0,0,0), 
                            'X1' = x$X1)
      } else if(x$trueT24<x$trueT23) {
        x_new <- data.frame('id' = rep(x$id, 6), 'from_state' = c(0, 0, 0,0, 2,2), 'to_state' = c(1,2,3,4,3,4), 
                            'transition' = c(1,2,3,4,8,9), 'Tstart' = c(rep(0, 4), rep(x$trueT02,2)), 
                            'Tstop' = c(rep(x$trueT02, 4), rep(x$trueT02+x$trueT24,2)), 'status' = c(0, 1, 0,0,0,1), 
                            'X1' = x$X1)
      } else {
        if (x$Ctimes < x$trueT02+x$trueT23+x$trueT34) {
          x_new <- data.frame('id' = rep(x$id, 7), 'from_state' = c(0, 0, 0,0, 2,2,3), 'to_state' = c(1,2,3,4,3,4,4), 
                              'transition' = c(1,2,3,4,8,9,10), 'Tstart' = c(rep(0, 4), rep(x$trueT02,2),x$trueT02+x$trueT23), 
                              'Tstop' = c(rep(x$trueT02, 4), rep(x$trueT02+x$trueT23,2),x$Ctimes), 'status' = c(0, 1, 0,0, 1,0,0), 
                              'X1' = x$X1)
        } else {
          x_new <- data.frame('id' = rep(x$id, 7), 'from_state' = c(0, 0, 0,0, 2,2,3), 'to_state' = c(1,2,3,4,3,4,4), 
                              'transition' = c(1,2,3,4,8,9,10), 'Tstart' = c(rep(0, 4), rep(x$trueT02,2),x$trueT02+x$trueT23), 
                              'Tstop' = c(rep(x$trueT02, 4), rep(x$trueT02+x$trueT23,2),x$trueT02+x$trueT23+x$trueT34), 'status' = c(0, 1, 0,0, 1,0,1), 
                              'X1' = x$X1)
        }
      }
    } else {
      if (x$Ctimes < x$trueT03+x$trueT34) {
        x_new <- data.frame('id' = rep(x$id, 5), 'from_state' = c(0, 0, 0,0,3), 'to_state' = c(1,2,3,4,4), 
                            'transition' = c(1,2,3,4,10), 'Tstart' = c(rep(0, 4), x$trueT03), 
                            'Tstop' = c(rep(x$trueT03, 4), x$Ctimes), 'status' = c(0, 0, 1,0, 0), 
                            'X1' = x$X1)
      } else {
        x_new <- data.frame('id' = rep(x$id, 5), 'from_state' = c(0, 0, 0,0,3), 'to_state' = c(1,2,3,4,4), 
                            'transition' = c(1,2,3,4,10), 'Tstart' = c(rep(0, 4), x$trueT03), 
                            'Tstop' = c(rep(x$trueT03, 4), x$trueT03+x$trueT34), 'status' = c(0, 0, 1,0, 1), 
                            'X1' = x$X1)
      }
    }
  }
}

data_mstate_split.by.id <- lapply(data_mstate_split.by.id, ms_arrange)
data_mstate <- do.call(rbind, data_mstate_split.by.id)
data_mstate$transition <- factor(data_mstate$transition)     
data_mstate$years<-data_mstate$Tstop-data_mstate$Tstart

################################################################################
# step 3: simulated longitudinal data
################################################################################
# par setting
DBRO<-1
sigma.re = sqrt(0.074)
# 
visits.out <- vector()
times<-tapply(data_mstate$Tstop, data_mstate$id, max) # 
ID <- longit.out <- vector()

for(i in 1:n){
  visits <- seq(0, times[i], by = DBRO)
  yt <- mu_long(visits, x=covdat[i,,drop=FALSE], betas = betas_01[i,]) + rnorm(length(visits),0,sigma.re)
  longit.out <- c(longit.out,yt)
  ID <- c(ID,rep(i,length(visits)))
  visits.out <- c(visits.out,visits)
}

long.proc <- as.matrix(cbind(longit.out,ID)) # Longitudinal process
obj <- list(long.proc,visits.out)
names(obj) <- c("longitudinal","visits")

save.image(paste0("data_",task_id,".RData"))