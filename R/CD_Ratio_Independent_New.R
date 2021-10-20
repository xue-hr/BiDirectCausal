CD_Ratio_Independent_New <- function (sig_part, V_T1 = NULL, V_T2 = NULL)
{
  p = nrow(sig_part)
  SNP = matrix(rnorm(p * p * 100), p * 100)
  SNP = scale(SNP)
  N_T1 = sig_part[, 8]
  N_T2 = sig_part[, 11]
  T1_T = sig_part[, 6]/sig_part[, 7]
  T1_r = T1_T/sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[, 9]/sig_part[, 10]
  T2_r = T2_T/sqrt(N_T2 - 2 + T2_T^2)
  rho_T1 = matrix(0, ncol = (p + 1), nrow = (p + 1))
  rho_T1[1:p, 1:p] = cor(SNP)
  rho_T1[1:p, p + 1] = T1_r
  rho_T1[p + 1, 1:p] = T1_r
  rho_T1[p + 1, p + 1] = 1
  rho_T2 = matrix(0, ncol = (p + 1), nrow = (p + 1))
  rho_T2[1:p, 1:p] = cor(SNP)
  rho_T2[1:p, p + 1] = T2_r
  rho_T2[p + 1, 1:p] = T2_r
  rho_T2[p + 1, p + 1] = 1
  if (is.null(V_T1)) {
    V_T1 = calculate_asymptotic_variance(SNP, rho_T1)
  }
  if (is.null(V_T2)) {
    V_T2 = calculate_asymptotic_variance(SNP, rho_T2)
  }
  jacobian = cbind(diag(1/T1_r), -diag(T2_r/T1_r^2))
  combined_V = rbind(cbind(V_T2, matrix(0, ncol = p, nrow = p))/mean(N_T2),
                     cbind(matrix(0, ncol = p, nrow = p), V_T1)/mean(N_T1))
  V = jacobian %*% combined_V %*% t(jacobian)
  inv_V = solve(V, tol = 0)
  est_vec = T2_r/T1_r
  gls_est = sum(inv_V %*% est_vec)/sum(inv_V)
  gls_var = 1/sum(inv_V)
  T1toT2 = c(gls_est, sqrt(gls_var))
  Q_T1toT2 = (est_vec - gls_est) %*% inv_V %*% (est_vec - gls_est)
  #jacobian = cbind(diag(1/T2_r), -diag(T1_r/T2_r^2))
  #combined_V = rbind(cbind(V_T1, matrix(0, ncol = p, nrow = p))/mean(N_T1),
  #                   cbind(matrix(0, ncol = p, nrow = p), V_T2)/mean(N_T2))
  #V = jacobian %*% combined_V %*% t(jacobian)
  #inv_V = solve(V, tol = 0)
  #est_vec = T1_r/T2_r
  #gls_est = sum(inv_V %*% est_vec)/sum(inv_V)
  #gls_var = 1/sum(inv_V)
  #T2toT1 = c(gls_est, sqrt(gls_var))
  #Q_T2toT1 = (est_vec - gls_est) %*% inv_V %*% (est_vec - gls_est)
  names(T1toT2) = c("K", "se(K)")
  #names(T2toT1) = c("K", "se(K)")
  return(list(T1toT2 = T1toT2, #T2toT1 = T2toT1,
              Q_T1toT2 = Q_T1toT2
              #Q_T2toT1 = Q_T2toT1
  ))
}
