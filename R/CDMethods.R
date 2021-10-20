CDMethods <- function(b_exp,b_out,
                      se_exp,se_out,
                      n_exp,n_out)
{
  r_exp = b_exp / sqrt(b_exp^2 + (n_exp-2)*se_exp^2)
  r_out = b_out / sqrt(b_out^2 + (n_out-2)*se_out^2)

  ### calculate covariance matrices of r_exp and r_out

  p = length(r_exp)
  SNP = matrix(rnorm(p*p*100),p*100)
  SNP = scale(SNP)
  #
  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T1[1:p,1:p] = diag(p)
  rho_T1[1:p,p+1] = r_exp
  rho_T1[p+1,1:p] = r_exp
  rho_T1[p+1,p+1] = 1
  #
  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T2[1:p,1:p] = diag(p)
  rho_T2[1:p,p+1] = r_out
  rho_T2[p+1,1:p] = r_out
  rho_T2[p+1,p+1] = 1
  #
  #V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  #V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  V_T1 = diag((1-r_exp^2)^2)
  V_T2 = diag((1-r_out^2)^2)

  sig_part = data.frame(chr = 1, pos = 1, rsid = "a", A1 = "A", A2 = "G",
                        beta_T1 = b_exp, se_T1 = se_exp, N_T1 = n_exp,
                        beta_T2 = b_out, se_T2 = se_out, N_T2 = n_out,
                        loci = 1)

  CD_Ratio_result =
    CD_Ratio_Independent_New(sig_part,V_T1 = V_T1,V_T2 = V_T2)

  CD_Egger_result =
    CD_Egger_Independent(sig_part,num_iteration = 20,
                         V_T1 = V_T1,V_T2 = V_T2)


  return(list(CD_Ratio_result = CD_Ratio_result,
              CD_Egger_result = CD_Egger_result)
  )
}
