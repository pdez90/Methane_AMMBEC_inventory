// inversion.stan — joint CH4 / C2H6 / CO source inversion for the AMMBEC flights -----------------
//
// Observations are segment means (ppb). For segment i on flight f(i), leg l(i):
//
//   CH4_i  = b_CH4[f] + tau[f] * sum_k alpha_k H_ik              + u_CH4[l]  + e
//   C2H6_i = b_C2H6[f] + tau[f] * sum_k alpha_k r_k H_ik         + u_C2H6[l] + e
//   CO_i   = b_CO[f]  + tau[f] * gamma * Hco_i                   + u_CO[l]   + e
//
//   alpha_k  scale factor on prior component k (lognormal prior)         <- the answer
//   r_k      ethane:methane ratio of component k (normal prior, or fixed at 0 for biogenic)
//   tau_f    flight-level transport/dilution factor shared by all three species (lognormal).
//            CH4 alone cannot separate tau_f from alpha; CO, whose inventory is well known
//            (gamma tightly around 1), pins tau_f down — the enhancement-ratio idea of Paper 1
//            inside the inversion.
//   b        flight backgrounds; u leg offsets (correlated error within a leg)
//   e        model-data mismatch, Student-t with nu degrees of freedom and scale
//            sqrt(s^2 + (phi * predicted enhancement_i)^2). nu = 4 (model v2) keeps unresolved point-source
//            plumes (the 10-h, 4-km model cannot place them) from inflating phi and pulling the signal
//            into the backgrounds; nu >= 1e6 is the Gaussian model v1.
//
// C2H6 and CO enter only through the index sets idx_c2h6 / idx_co, so passing N_c2h6 = 0 or
// N_co = 0 gives the CH4-only and CH4+C2H6 versions of the same model (identifiability tests).
// ------------------------------------------------------------------------------------------------
data {
  int<lower=1> N;  int<lower=1> K;  int<lower=1> F;  int<lower=1> L;
  array[N] int<lower=1, upper=F> flight;
  array[N] int<lower=1, upper=L> leg;
  matrix[N, K] H;                  // ppb per unit alpha
  vector[N] Hco;                   // ppb, prior CO
  vector[N] y_ch4;
  int<lower=0> N_c2h6; array[N_c2h6] int<lower=1, upper=N> idx_c2h6; vector[N_c2h6] y_c2h6;
  int<lower=0> N_co;   array[N_co] int<lower=1, upper=N> idx_co;     vector[N_co] y_co;

  vector[K] la_mu;  vector<lower=0>[K] la_sd;            // log(alpha) prior
  vector<lower=0>[K] r_mu;  vector<lower=0>[K] r_sd;     // ethane ratio prior (sd 0 = fixed)
  int<lower=0, upper=K> K_rfree; array[K_rfree] int<lower=1, upper=K> rfree;

  vector[F] b_ch4_mu;  vector[F] b_c2h6_mu;  vector[F] b_co_mu;
  real<lower=0> b_ch4_sd;  real<lower=0> b_c2h6_sd;  real<lower=0> b_co_sd;
  real<lower=0> tau_sd;  real<lower=0> gamma_sd;
  real<lower=0> s_ch4_scale;  real<lower=0> s_c2h6_scale;  real<lower=0> s_co_scale;
  real<lower=0> leg_ch4_scale;  real<lower=0> leg_c2h6_scale;  real<lower=0> leg_co_scale;
  real<lower=0> phi_scale;
  real<lower=1> nu;                         // Student-t degrees of freedom (1e6 = Gaussian)

  vector<lower=0>[K] E_box;                 // prior emissions in the Paper 1 box (t CH4 / h)
  array[K] int<lower=0, upper=1> is_fossil;
}
parameters {
  vector[K] la_z;
  vector<lower=0>[K_rfree] r_free;
  vector[F] b_ch4_z;  vector[F] b_c2h6_z;  vector[F] b_co_z;
  vector[F] ltau_z;  real lgamma_z;
  real<lower=0> s_ch4;  real<lower=0> s_c2h6;  real<lower=0> s_co;  real<lower=0> phi;
  vector[L] u_ch4_z;  vector[L] u_c2h6_z;  vector[L] u_co_z;
  real<lower=0> sl_ch4;  real<lower=0> sl_c2h6;  real<lower=0> sl_co;
}
transformed parameters {
  vector[K] alpha = exp(la_mu + la_sd .* la_z);
  vector[K] r = r_mu;
  vector[F] tau = exp(tau_sd * ltau_z);
  real gamma = exp(gamma_sd * lgamma_z);
  r[rfree] = r_free;
}
model {
  vector[N] enh_ch4 = tau[flight] .* (H * alpha);
  // priors
  la_z ~ std_normal();
  r_free ~ normal(r_mu[rfree], r_sd[rfree]);
  b_ch4_z ~ std_normal();  b_c2h6_z ~ std_normal();  b_co_z ~ std_normal();
  ltau_z ~ std_normal();  lgamma_z ~ std_normal();
  s_ch4 ~ normal(0, s_ch4_scale);  s_c2h6 ~ normal(0, s_c2h6_scale);  s_co ~ normal(0, s_co_scale);
  phi ~ normal(0, phi_scale);
  u_ch4_z ~ std_normal();  u_c2h6_z ~ std_normal();  u_co_z ~ std_normal();
  sl_ch4 ~ normal(0, leg_ch4_scale);  sl_c2h6 ~ normal(0, leg_c2h6_scale);  sl_co ~ normal(0, leg_co_scale);
  // CH4
  y_ch4 ~ student_t(nu, b_ch4_mu[flight] + b_ch4_sd * b_ch4_z[flight] + enh_ch4 + sl_ch4 * u_ch4_z[leg],
                 sqrt(square(s_ch4) + square(phi * enh_ch4)));
  // C2H6
  if (N_c2h6 > 0) {
    vector[N_c2h6] enh = tau[flight[idx_c2h6]] .* (H[idx_c2h6] * (alpha .* r));
    y_c2h6 ~ student_t(nu, b_c2h6_mu[flight[idx_c2h6]] + b_c2h6_sd * b_c2h6_z[flight[idx_c2h6]] + enh
                    + sl_c2h6 * u_c2h6_z[leg[idx_c2h6]], sqrt(square(s_c2h6) + square(phi * enh)));
  }
  // CO
  if (N_co > 0) {
    vector[N_co] enh = gamma * tau[flight[idx_co]] .* Hco[idx_co];
    y_co ~ student_t(nu, b_co_mu[flight[idx_co]] + b_co_sd * b_co_z[flight[idx_co]] + enh
                  + sl_co * u_co_z[leg[idx_co]], sqrt(square(s_co) + square(phi * enh)));
  }
}
generated quantities {
  vector[K] E_post = alpha .* E_box;                  // t CH4 / h in the Paper 1 box
  real E_total = sum(E_post);
  real fossil_share;
  {
    real fos = 0;
    for (k in 1:K) if (is_fossil[k] == 1) fos += E_post[k];
    fossil_share = E_total > 0 ? fos / E_total : 0;
  }
}
