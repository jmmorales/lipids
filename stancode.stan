functions { 
  
  /* compute the kronecker product
  * Copied from brms: Paul-Christian Bürkner (2018). 
  * Advanced Bayesian Multilevel Modeling with the R Package brms. 
  * The R Journal, 10(1), 395-411. <doi:10.32614/RJ-2018-017>
  * Args: 
    *   A,B: matrices 
  * Returns: 
    *   kronecker product of A and B
  */ 
    matrix kronecker(matrix A, matrix B) { 
      matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
      for (i in 1:cols(A)) { 
        for (j in 1:rows(A)) { 
          kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
        } 
      } 
      return kron; 
    } 
  
  // copied from R. McElreath (2018). 
  // Statistical rethinking: A Bayesian course with examples in R and Stan. 
   matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }    
} 

data { 
  int<lower=1> N;                 // total number of observations 
  int Y[N];                       // response variable 
  int<lower=1> K;                 // number of population-level effects (2)
  int<lower=1> J;                 // num of groups (bird gens)
  int<lower=1> L;                 // num group level predictors
  int<lower=1,upper=J> jj[N];     // group id 
  matrix[N, K] X;                 // obs-level design matrix 
  matrix[J,L] TT;                 // group-level traits
  matrix[J,J] C;                  // phylogenetic correlation matrix
  vector[J] ones;                 // vector on 1s
  int<lower=1> N_1;               // num sites
  int<lower=1> M_1;               // one site r.e.
  int<lower=1> M_2;               // one plant r.e.
  int<lower=1,upper=N_1> J_1[N];  // ids for sites
  int<lower=1> J_2[N];            // ids for plants
  int<lower=1> N_2;               // num plant spp
  matrix[N_2, N_2] DistP;         // cophenetic distance among plants
  matrix[N_1,N_1] Dmat;           // sites distance matrix
}

parameters {
  corr_matrix[K] Omega;
  vector<lower=0>[K] tau;
  vector[J * K] beta;
  real<lower=0,upper=1> rho;
  vector[L * K] z;
  vector[N_1] r_1_1;          //  site-level effects
  vector[N_2] r_1_2;          // plqnt-level effects
  real<lower=0> etasq;
  real<lower=0> rhosq;
  real<lower=0> etasqp;
  real<lower=0> rhosqp;
  real<lower=0> delta;
  real<lower=0> deltap;
  real<lower=0> phi;
  real<lower=0,upper=1> theta;
}

transformed parameters { 
  matrix[K, K] Sigma = quad_form_diag(Omega, tau);
  matrix[J+J, J+J] S = kronecker(Sigma, rho * C + (1-rho) * diag_matrix(ones));
  matrix[L, K] Z = to_matrix(z, L, K);
  vector[J * K] m = to_vector(TT * Z);
  matrix[J, K] b_m = to_matrix(beta, J, K);
} 

model {
  matrix[N_1,N_1] SIGMA;
  matrix[N_2,N_2] SIGMAP;
  
  rhosq ~ exponential( 0.5 );
  etasq ~ exponential( 2 );
  delta ~ normal(0, 2.5);
  rhosqp ~ exponential( 0.5 );
  etasqp ~ exponential( 2 );
  deltap ~ normal(0, 2.5);
    
  Omega ~ lkj_corr(2);
  tau ~ student_t(3,0,10); // cauchy(0, 2.5);
  beta ~ multi_normal(m, S);
  //rho ~ beta(1,10);
  z ~ normal(0,1);
  
  phi ~ cauchy(0, 3);
  
  SIGMA = cov_GPL2(Dmat, etasq, rhosq, delta);
  r_1_1 ~ multi_normal( rep_vector(0,N_1) , SIGMA );
  
  SIGMAP = cov_GPL2(DistP, etasqp, rhosqp, deltap);
  r_1_2 ~ multi_normal( rep_vector(0,N_2) , SIGMAP );

 {
  vector[N] mu;

  for (n in 1:N){
     mu[n] = X[n,1] * b_m[jj[n],1] + X[n,2] * b_m[jj[n],2] + r_1_1[J_1[n]] + r_1_2[J_2[n]];
  }
  
  Y ~ poisson_log(mu);
  Y ~ neg_binomial_2(exp(mu), phi);
  target += log_sum_exp(log(theta) +  beta_lpdf(rho|1, 100), log(1-theta) +  beta_lpdf(rho|2, 2));
  }
}