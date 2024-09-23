// stan code fitting joint uni transition model

functions{
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
    vector linear_predictor(vector visits, int[] ID, vector beta_tilde, matrix bi){ 
         int N = num_elements(visits);
         vector[N] out;

         out = beta_tilde[1] + beta_tilde[2]*visits + bi[ID,1] + rows_dot_product(bi[ID,2],visits);

         return out;
    }
    
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE single transition SUBMODEL (time scale since Tstart)                
// ------------------------------------------------------ 
    vector linear_predictor_S(matrix X, vector t_cf_event_q, int[] ID_q, vector beta_tilde, real gamma, matrix bi, real alpha){
      int N_q = num_elements(t_cf_event_q);
      vector[N_q] out;
      
      out = X[ID_q,1]*gamma + alpha*(beta_tilde[1] + beta_tilde[2]*t_cf_event_q + bi[ID_q,1] + rows_dot_product(bi[ID_q,2],t_cf_event_q));
      
      return out;
    }
    


// ----- baseline log-hazard for a specific transition (here use weibull baseline)
    vector baseline_loghaz(vector t_event_q, real delta, real lambda){
      int N_q = num_elements(t_event_q);
      vector[N_q] out;
      
      out = log(delta)+log(lambda)+(delta-1)*log(t_event_q);
      
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
  vector[nrow_q] t_cf_event_q; // clock forward version of t_event_q
  vector[Qn] qwts;
  int<lower=1,upper=n> ID[N]; // ID for longitudinal submodel
  int<lower=1,upper=n> ID_q[nrow_q]; // ID for MSM submodel
  vector[N] visits; // visit times 
  matrix[n,1] X; //design matrix for time-independent covariates in each MSM transition
}

parameters{
  vector[2] beta_tilde;
  real<lower=0> Var_e;  
  matrix[n,2] bi;
  real<lower=0> delta;
  real<lower=0> lambda;
  real alpha;
  real gamma;
  cholesky_factor_corr[2] L_Omega; // Cholesky factor of a 3x3 correlation matrix for Sigma
  vector<lower=0>[2] tau;          // Standard deviations for Sigma
}

transformed parameters{
  matrix[2,2] Sigma;               // Covariance matrix for random effect
  matrix[2,2] Omega;
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
 
 log_basehaz = baseline_loghaz(t_event_q, delta, lambda);
 
 lp_s = linear_predictor_S(X, t_cf_event_q, ID_q, beta_tilde, gamma, bi, alpha);
 
 log_haz_all = log_basehaz + lp_s;
 
 log_haz_etimes = head(log_haz_all,Nevents);
 
 log_haz_qtimes = tail(log_haz_all,Qn);
 //
 target += sum(log_haz_etimes)-dot_product(qwts,exp(log_haz_qtimes));   
}

// --- log-priors-----------------------------------
  target += normal_lpdf(beta_tilde | 0, 100);
  target += normal_lpdf(gamma | 0, 100);
  target += normal_lpdf(alpha | 0, 100);
  target += cauchy_lpdf(delta | 0, 1);
  target += cauchy_lpdf(lambda | 0, 1);
  target += inv_gamma_lpdf(Var_e | 0.01, 0.01); 

  // Random-effects variance-covariance matrices
  L_Omega ~ lkj_corr_cholesky(2);  // LKJ prior on the Cholesky factor
  tau ~ cauchy(0, 2.5);            // Prior for standard deviations

  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:2] | rep_vector(0,2), Sigma); }
}
