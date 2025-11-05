rm(list = ls(all = TRUE))
set.seed(1973)

# --- Settings ---
n = 100
p = 15
q = 10
alpha = 0.05
xrho = 0.95          # rho in Eq (18)
tau = 1
iter = 1000

d = 0.1
delta0 = 0

Ip = diag(p)
eps = 1e-10
standardize_X = FALSE   # Sec. 8.1 Eq (18) typically uses raw X; set TRUE only if you want scaling

# --- Helper: draw X via Sec. 8.1, Eq. (18) ---
draw_X_eq18 = function(n, p, rho) {
  T1 = matrix(rnorm(n * p), n, p)        # T_ij ~ N(0,1)
  T2 = matrix(rnorm(n), n, 1)            # T_{i,p+1} ~ N(0,1), common across columns for row i
  X  = sqrt(1 - rho^2) * T1 + rho * T2 %*% matrix(1, 1, p)
  X
}

make_beta = function(delta) c(rep(1, p - q), rep(0, q - 1), delta)

# --- KL component-k helper (unchanged) ---
kl_component_k = function(alpha_i, lambda_i, s2) {
  if (lambda_i < eps) lambda_i = eps
  k_i = s2 / (2 * (alpha_i^2) + s2 / lambda_i)
  max(k_i, 0)
}

# --- Step 1: Choose k via KL heuristic (average across replications) using Eq (18) X ---
k_opt_vec = numeric(iter)
for (t in 1:iter) {
  X = draw_X_eq18(n, p, xrho)
  if (standardize_X) X = scale(X, center = TRUE, scale = TRUE)

  y = as.vector(X %*% make_beta(delta0) + rnorm(n, 0, tau))

  fit = lm(y ~ X - 1)
  s2  = as.numeric(summary(fit)$sigma)^2

  G   = crossprod(X)
  EV  = eigen(G, symmetric = TRUE)
  Q   = EV$vectors
  lam = pmax(EV$values, eps)

  Z         = X %*% Q
  alpha_hat = as.vector(crossprod(Z, y) / lam)

  k_i = mapply(kl_component_k, alpha_hat, lam, MoreArgs = list(s2 = s2))
  k_opt_vec[t] = min(k_i)
}
k = mean(k_opt_vec)
cat("Selected k (KL heuristic avg, Eq. 18 X) =", round(k, 6), "\n")

# --- Results table (baseline = KL) with parameters included ---
delta_grid = seq(0, 2, by = 0.1)
RE_table = data.frame(
  delta = delta_grid,
  RE_RKL_vsKL = NA_real_,
  RE_KLPT_vsKL = NA_real_,
  RE_KLS_vsKL = NA_real_,
  RE_KLPS_vsKL = NA_real_,
  RE_LSKL_vsKL = NA_real_,
  RE_SPTKL_vsKL = NA_real_,
  n = n, p = p, q = q, k = round(k, 5), d = d, rho = xrho
)

# --- Simulation over delta values (again using Eq (18) X) ---
for (dd in seq_along(delta_grid)) {
  delta = delta_grid[dd]
  thetatrue = make_beta(delta)

  mseKL    = numeric(iter)
  mseRKL   = numeric(iter)
  mseKLPT  = numeric(iter)
  mseKLS   = numeric(iter)
  mseKLPS  = numeric(iter)
  mseLSKL  = numeric(iter)
  mseSPTKL = numeric(iter)

  R = cbind(matrix(0, nrow = q, ncol = p - q), diag(q))
  lalpha = qf(1 - alpha, q, n - p)

  for (i in 1:iter) {
    X = draw_X_eq18(n, p, xrho)
    if (standardize_X) X = scale(X, center = TRUE, scale = TRUE)

    y = as.vector(X %*% thetatrue + rnorm(n, 0, tau))

    G    = crossprod(X)
    Ginv = solve(G + eps * Ip)

    fit = lm(y ~ X - 1)
    thetaOLS = coef(fit)
    se = drop(crossprod(y - X %*% thetaOLS)) / (n - p)

    # Restricted OLS under R theta = 0
    thetaR = thetaOLS - Ginv %*% t(R) %*% solve(R %*% Ginv %*% t(R)) %*% (R %*% thetaOLS)

    # KL & RKL
    thetaKL  = solve(G + k * Ip) %*% ((G - k * Ip) %*% thetaOLS)
    thetaRKL = solve(G + k * Ip) %*% ((G - k * Ip) %*% thetaR)

    # Pretest statistic
    phin = drop(t(R %*% thetaOLS) %*% solve(R %*% Ginv %*% t(R)) %*% (R %*% thetaOLS) / (q * se))

    # KLPT
    thetaKLPT = if (phin <= lalpha) thetaRKL else thetaKL

    # Stein-type KLS
    piconstan = (q - 2) * (n - p) / (q * (n - p + 2))
    thetaKLS  = thetaKL - (thetaKL - thetaRKL) * piconstan / phin

    # Positive-rule KLPS
    thetaKLPS = if (phin <= piconstan) {
      thetaKLS - (thetaKL - thetaRKL) * (1 - piconstan / phin)
    } else thetaKLS

    # LSKL
    thetaLSKL = thetaKL - d * (thetaKL - thetaRKL)

    # SPTKL (pretest + Liu-step)
    thetaSPTKL = if (phin <= lalpha) {
      thetaKL - d * (thetaKL - thetaRKL)
    } else thetaKL

    # Parameter-risk MSEs
    mseKL[i]    = sum((thetaKL   - thetatrue)^2)
    mseRKL[i]   = sum((thetaRKL  - thetatrue)^2)
    mseKLPT[i]  = sum((thetaKLPT - thetatrue)^2)
    mseKLS[i]   = sum((thetaKLS  - thetatrue)^2)
    mseKLPS[i]  = sum((thetaKLPS - thetatrue)^2)
    mseLSKL[i]  = sum((thetaLSKL - thetatrue)^2)
    mseSPTKL[i] = sum((thetaSPTKL - thetatrue)^2)
  }

  # Relative efficiencies vs KL (baseline = KL)
  S_KL = sum(mseKL)
  RE_table$RE_RKL_vsKL[dd]   = S_KL / sum(mseRKL)
  RE_table$RE_KLPT_vsKL[dd]  = S_KL / sum(mseKLPT)
  RE_table$RE_KLS_vsKL[dd]   = S_KL / sum(mseKLS)
  RE_table$RE_KLPS_vsKL[dd]  = S_KL / sum(mseKLPS)
  RE_table$RE_LSKL_vsKL[dd]  = S_KL / sum(mseLSKL)
  RE_table$RE_SPTKL_vsKL[dd] = S_KL / sum(mseSPTKL)

  cat("iter=", iter, " delta=", delta, " k=", round(k,4), " d=", d, " rho=", xrho,
      " n=", n, " p=", p, " q=", q,
      " | RE vs KL -> RKL:", round(RE_table$RE_RKL_vsKL[dd],3),
      " KLPT:", round(RE_table$RE_KLPT_vsKL[dd],3),
      " KLS:", round(RE_table$RE_KLS_vsKL[dd],3),
      " KLPS:", round(RE_table$RE_KLPS_vsKL[dd],3),
      " LSKL:", round(RE_table$RE_LSKL_vsKL[dd],3),
      " SPTKL:", round(RE_table$RE_SPTKL_vsKL[dd],3), "\n")
}

# Final table
RE_table
