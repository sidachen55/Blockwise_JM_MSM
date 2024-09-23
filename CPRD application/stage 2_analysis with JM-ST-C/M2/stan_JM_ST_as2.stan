// stan code fitting joint uni transition model (JM-ST approach) - for real data application 
// stage 2_JM-ST-C analysis
// current value + age-sbp interaction (as2)

functions{
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
    vector linear_predictor(vector visits, int[] ID, vector beta_tilde, matrix bi){ 
         int N = num_elements(visits);
         vector[N] out = beta_tilde[1] + beta_tilde[2]*visits + bi[ID,1] + bi[ID,2].*visits;
         return out;
    }
    
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE multistate SUBMODEL (time scale since origin)                
// ------------------------------------------------------ 
    vector linear_predictor_S(matrix X, vector t_cf_event_q, int[] ID_q, int[] X_q, vector beta_tilde, real gamma1, real gamma2, real gamma3, real gamma4, matrix bi, vector alpha){
      int N_q = num_elements(t_cf_event_q);
      vector[N_q] lp= beta_tilde[1] + beta_tilde[2]*t_cf_event_q + bi[ID_q,1] + bi[ID_q,2].*t_cf_event_q;
      vector[N_q] out = X[X_q,1]*gamma1 + X[X_q,2]*gamma2 + X[X_q,3]*gamma3 + X[X_q,4]*gamma4 + alpha[1]*lp + alpha[2]*(X[X_q,1].*lp);
      return out;
    }
    


// ----- baseline log-hazard for a specific transition
    vector baseline_loghaz(vector B_1_event_q, vector B_2_event_q, vector B_3_event_q, vector B_4_event_q, vector B_5_event_q, real eta_1, real eta_2, real eta_3, real eta_4, real eta_5){
      int N_q = num_elements(B_1_event_q);
      vector[N_q] out = eta_1*B_1_event_q + eta_2*B_2_event_q + eta_3*B_3_event_q + eta_4*B_4_event_q + eta_5*B_5_event_q;
      return out;
    }
}


data{
  int N; // total number of longitudinal outcomes
  int n; // total number of subjects
  int Nevents; // total number of observed transitions
  int Qn; 
  int nrow_msmdata;
  int nrow_q; // Nevents+Qn
  vector[N] y; // longitudinal outcome
  vector[nrow_q] t_cf_event_q; // clock forward version of t_event_q
  vector[Qn] qwts;
  int<lower=1,upper=n> ID[N]; // ID for longitudinal submodel
  int<lower=1,upper=n> ID_q[nrow_q]; // ID for MSM submodel
  vector[N] visits; // visit times 
  int<lower=1,upper=nrow_q> X_q[nrow_q]; // row index for cov X
  matrix[nrow_msmdata,4] X; //design matrix for covariates in each MSM transition
  vector[nrow_q] B_1_event_q;
  vector[nrow_q] B_2_event_q;
  vector[nrow_q] B_3_event_q;
  vector[nrow_q] B_4_event_q;
  vector[nrow_q] B_5_event_q;
  matrix[n, 15] B_1_gq_q; // for computing pointwise llk in gq block
  matrix[n, 15] B_2_gq_q;
  matrix[n, 15] B_3_gq_q;
  matrix[n, 15] B_4_gq_q;
  matrix[n, 15] B_5_gq_q;
  vector[n] B_1_gq_haz; // for computing pointwise llk in gq block
  vector[n] B_2_gq_haz;
  vector[n] B_3_gq_haz;
  vector[n] B_4_gq_haz;
  vector[n] B_5_gq_haz; 
  int<lower=1,upper=N> id_lg_min[n]; // for producing pointwise llk
  int<lower=1,upper=N> id_lg_n[n]; // for producing pointwise llk
  vector[n] T_e; // event times T_i 
  vector[15] xk; // quadrature 
  vector[15] wk; // quadrature 
  vector[n] T_status; //event indicator delta_i
}

parameters{
  vector[2] beta_tilde;
  real<lower=0> Var_b[2];
  real<lower=-1, upper=1> rho; 
  real<lower=0> Var_e;  
  matrix[n,2] bi;
  vector[2] alpha;
  real gamma1;
  real gamma2;
  real gamma3;
  real gamma4;
  real eta_1;
  real eta_2;
  real eta_3;
  real eta_4;
  real eta_5;
}

transformed parameters{
  cov_matrix[2] Sigma;
  Sigma[1,1] = Var_b[1];
  Sigma[2,2] = Var_b[2];
  Sigma[1,2] = sqrt(Var_b[1]*Var_b[2])*rho;
  Sigma[2,1] = Sigma[1,2];
}

model{
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR longitudinal (univariate) SUBMODEL                
// ------------------------------------------------------
{
   vector[N] linpred; 

   // Linear predictor
   linpred = linear_predictor(visits, ID, beta_tilde, bi);

   // Longitudinal Normal log-likelihood
   target += normal_lpdf(y | linpred, sqrt(Var_e));
}

// ------------------------------------------------------
//          LOG-LIKELIHOOD FOR multistate SUBMODEL  (Gauss quadrature)              
// ------------------------------------------------------
{
 vector[Nevents] log_haz_etimes;
 vector[Qn] log_haz_qtimes;
 vector[nrow_q] log_haz_all;
 vector[nrow_q] log_basehaz;
 vector[nrow_q] lp_s;
 
 log_basehaz = baseline_loghaz(B_1_event_q, B_2_event_q, B_3_event_q, B_4_event_q, B_5_event_q, eta_1, eta_2, eta_3, eta_4, eta_5);
 
 lp_s = linear_predictor_S(X, t_cf_event_q, ID_q, X_q, beta_tilde, gamma1, gamma2, gamma3, gamma4, bi, alpha);
 
 log_haz_all = log_basehaz + lp_s;
 
 log_haz_etimes = head(log_haz_all,Nevents);
 
 log_haz_qtimes = tail(log_haz_all,Qn);
 //
 target += sum(log_haz_etimes)-dot_product(qwts,exp(log_haz_qtimes));   
}

// --- log-priors-----------------------------------
  target += normal_lpdf(beta_tilde | 0, 10);
  target += normal_lpdf(gamma1 | 0, 10);
  target += normal_lpdf(gamma2 | 0, 10);
  target += normal_lpdf(gamma3 | 0, 10);
  target += normal_lpdf(gamma4 | 0, 10);
  target += normal_lpdf(alpha | 0, 10);
  target += normal_lpdf(eta_1 | 0, 10);
  target += normal_lpdf(eta_2 | 0, 10);
  target += normal_lpdf(eta_3 | 0, 10);
  target += normal_lpdf(eta_4 | 0, 10);
  target += normal_lpdf(eta_5 | 0, 10);
  
  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:2] | rep_vector(0.0,2), Sigma); }
  
  target += inv_gamma_lpdf(Var_b | 0.01, 0.01);
  target += beta_lpdf((rho+1)/2 | 0.5, 0.5); 
  target += inv_gamma_lpdf(Var_e | 0.01, 0.01);   
}  

generated quantities {
  vector[n] log_lik_s;
  vector[n] log_lik_l;
  vector[n] log_lik;
  vector[n] cumhaz;
  vector[15] cumhaz_q;
  vector[n] log_haz;
  vector[N] linpred;
  vector[N] llk_lg_all;
  for (i in 1:n) {
    for (j in 1:15) {cumhaz_q[j] = (T_e[i]/2)*wk[j]*(exp(eta_1*B_1_gq_q[i,j] + eta_2*B_2_gq_q[i,j] + eta_3*B_3_gq_q[i,j] + eta_4*B_4_gq_q[i,j] + eta_5*B_5_gq_q[i,j] + X[i,1]*gamma1 + X[i,2]*gamma2 + X[i,3]*gamma3 + X[i,4]*gamma4 + alpha[1]*(beta_tilde[1] + beta_tilde[2]*T_e[i]*(1+xk[j])/2 + bi[i,1] + bi[i,2]*T_e[i]*(1+xk[j])/2) + alpha[2]*X[i,1]*(beta_tilde[1] + beta_tilde[2]*T_e[i]*(1+xk[j])/2 + bi[i,1] + bi[i,2]*T_e[i]*(1+xk[j])/2)));}
    cumhaz[i] = sum(cumhaz_q);
    }
  for (i in 1:n) {log_haz[i] = T_status[i]*(eta_1*B_1_gq_haz[i] + eta_2*B_2_gq_haz[i] + eta_3*B_3_gq_haz[i] + eta_4*B_4_gq_haz[i] + eta_5*B_5_gq_haz[i] + X[i,1]*gamma1 + X[i,2]*gamma2 + X[i,3]*gamma3 + X[i,4]*gamma4 + alpha[1]*(beta_tilde[1] + beta_tilde[2]*T_e[i] + bi[i,1] + bi[i,2]*T_e[i]) + alpha[2]*X[i,1]*(beta_tilde[1] + beta_tilde[2]*T_e[i] + bi[i,1] + bi[i,2]*T_e[i]));}
  for (i in 1:n) {log_lik_s[i] = log_haz[i]-cumhaz[i];}
  for (i in 1:N) {linpred[i] = beta_tilde[1] + beta_tilde[2]*visits[i] + bi[ID[i],1] + bi[ID[i],2]*visits[i];}
  for (i in 1:N) {llk_lg_all[i] = normal_lpdf(y[i] | linpred[i], sqrt(Var_e));}
  for (i in 1:n) {log_lik_l[i] = sum(segment(llk_lg_all, id_lg_min[i], id_lg_n[i]));}
  for (i in 1:n) {log_lik[i] = log_lik_s[i]+log_lik_l[i];}
}
