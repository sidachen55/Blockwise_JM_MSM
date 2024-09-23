//------------------------------------------------------------------------------
// stan code fitting joint multistate model using the JM-CR approach (clock forward scale, updated 2024.3.5)
// block1 (trans=1,2,3,4)
//------------------------------------------------------------------------------
functions{
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
    vector linear_predictor(matrix X, vector f1_1, vector f2_1, int[] ID, vector beta, vector beta_x, matrix bi){ 
         int N = num_elements(f1_1);
         vector[N] out;

         out = beta[1] + beta_x[1]*X[ID,1] + bi[ID,1] + rows_dot_product(beta[2] + beta_x[2]*X[ID,1] + bi[ID,2],f1_1) + rows_dot_product(beta[3] + beta_x[3]*X[ID,1] + bi[ID,3],f2_1);

         return out;
    }
    
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE multistate SUBMODEL (time scale since origin)                
// ------------------------------------------------------ 
    vector linear_predictor_S(matrix X, vector f1_2, vector f2_2, int[] ID_q, int[] trans_q, vector beta, vector beta_x, vector gamma, matrix bi, vector alpha){
      int N_q = num_elements(f1_2);
      vector[N_q] out;
      vector[N_q] out0;
      out0 = beta[1] + beta_x[1]*X[ID_q,1] + bi[ID_q,1] + rows_dot_product(beta[2] + beta_x[2]*X[ID_q,1] + bi[ID_q,2],f1_2) + rows_dot_product(beta[3] + beta_x[3]*X[ID_q,1] + bi[ID_q,3],f2_2);
      out = rows_dot_product(X[ID_q,1],gamma[trans_q]) + rows_dot_product(alpha[trans_q],out0);
      
      return out;
    }
    


// ----- baseline log-hazard for a specific transition (here use weibull baseline)
    vector baseline_loghaz(vector B_1_event_q, vector B_2_event_q, vector B_3_event_q, vector B_4_event_q, vector B_5_event_q, int[] trans_q, vector eta_1, vector eta_2, vector eta_3, vector eta_4, vector eta_5){
      int N_q = num_elements(B_1_event_q);
      vector[N_q] out;
      
      out = rows_dot_product(eta_1[trans_q],B_1_event_q) + rows_dot_product(eta_2[trans_q],B_2_event_q) + rows_dot_product(eta_3[trans_q],B_3_event_q) + rows_dot_product(eta_4[trans_q],B_4_event_q) + rows_dot_product(eta_5[trans_q],B_5_event_q);
      
      return out;
    }
}


data{
  int N; // total number of longitudinal outcomes
  int n; // total number of subjects
  int Nevents; // total number of observed transitions
  int Qn; 
  int nrow_q; // Nevents+Qn
  vector[N] y; // longitudinal outcome
  vector[nrow_q] t_event_q; 
  vector[Qn] qwts;
  int<lower=1,upper=n> ID[N]; // ID for longitudinal submodel
  int<lower=1,upper=n> ID_q[nrow_q]; // ID for MSM submodel
  int<lower=1,upper=4> trans_q[nrow_q]; // transition type indicator
  matrix[n,1] X; //design matrix for time-independent covariates in each MSM transition
  vector[nrow_q] B_1_event_q;
  vector[nrow_q] B_2_event_q;
  vector[nrow_q] B_3_event_q;
  vector[nrow_q] B_4_event_q;
  vector[nrow_q] B_5_event_q;
  vector[N] f1_1;
  vector[N] f2_1;
  vector[nrow_q] f1_2;
  vector[nrow_q] f2_2;
}

parameters{
  vector[3] beta;
  vector[3] beta_x;
  real<lower=0> Var_e;  
  matrix[n,3] bi;
  vector[4] eta_1;
  vector[4] eta_2;
  vector[4] eta_3;
  vector[4] eta_4;
  vector[4] eta_5;
  vector[4] alpha;
  vector[4] gamma;
  cholesky_factor_corr[3] L_Omega; // Cholesky factor of a 3x3 correlation matrix for Sigma
  vector<lower=0>[3] tau;          // Standard deviations for Sigma
}

transformed parameters{
  matrix[3,3] Sigma;               // Covariance matrix for random effect
  matrix[3,3] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  Sigma = quad_form_diag(Omega, tau);
}

model{
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR longitudinal (univariate) SUBMODEL                
// ------------------------------------------------------
{
   vector[N] linpred; 

   // Linear predictor
   linpred = linear_predictor(X, f1_1, f2_1, ID, beta, beta_x, bi);

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
 
 log_basehaz = baseline_loghaz(B_1_event_q, B_2_event_q, B_3_event_q, B_4_event_q, B_5_event_q, trans_q, eta_1, eta_2, eta_3, eta_4, eta_5);
 
 lp_s = linear_predictor_S(X, f1_2, f2_2, ID_q, trans_q, beta, beta_x, gamma, bi, alpha);
 
 log_haz_all = log_basehaz + lp_s;
 
 log_haz_etimes = head(log_haz_all,Nevents);
 
 log_haz_qtimes = tail(log_haz_all,Qn);
 //
 target += sum(log_haz_etimes)-dot_product(qwts,exp(log_haz_qtimes));  
}

// --- log-priors-----------------------------------
  target += normal_lpdf(beta | 0, 100);
  target += normal_lpdf(beta_x | 0, 100);
  target += normal_lpdf(gamma | 0, 100);
  target += normal_lpdf(alpha | 0, 100); 
  target += normal_lpdf(eta_1 | 0, 10);
  target += normal_lpdf(eta_2 | 0, 10);
  target += normal_lpdf(eta_3 | 0, 10);
  target += normal_lpdf(eta_4 | 0, 10);
  target += normal_lpdf(eta_5 | 0, 10);
  target += inv_gamma_lpdf(Var_e | 0.01, 0.01);
  
  // Random-effects variance-covariance matrices
  L_Omega ~ lkj_corr_cholesky(2);  // LKJ prior on the Cholesky factor
  tau ~ cauchy(0, 2.5);            // Prior for standard deviations
  
  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:3] | rep_vector(0,3), Sigma); }
}  
