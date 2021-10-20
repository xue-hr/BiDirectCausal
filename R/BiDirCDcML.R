#' Bi-directional CDcML Method
#'
#' @param b_X Vector of estimated effect sizes for trait 1.
#' @param b_Y Vector of estimated effect sizes for trait 2.
#' @param se_X Vector of standard errors of b_X.
#' @param se_Y Vector of standard errors of b_Y.
#' @param n_X Sample size for trait 1.
#' @param n_Y Sample size for trait 2.
#' @param sig.cutoff Significant cutoff to choose IV, defaul is 5e-8.
#' @param num_pert Number of perturbations, default is 100.
#' @param random_start Number of starting points for solving non-convex optimization
#' in cML method, defaul is 0, i.e. no random starting point.
#' @param random.seed Random seed, default is 0, i.e. no random seed is set.
#'
#' @return List of estimated correlation ratio (est) and
#' corresponding standard error (se) obtained with CDcML, for
#' both directions. "NoS" indicates no screening, "S" indicates with screening.
#' "DP" indicates results with data perturbation, otherwise without data perturbation.
#' @export
#'
#' @examples
BiDirCDcML <- function(b_X,b_Y,se_X,se_Y,n_X,n_Y,
                       sig.cutoff = 5e-8, num_pert = 100,
                       random_start = 0, random.seed = 0
                       )
{
  pvalue.X = pnorm(-abs(b_X/se_X))*2
  pvalue.Y = pnorm(-abs(b_Y/se_Y))*2

  ind_X = which(pvalue.X<(sig.cutoff))
  ind_Y = which(pvalue.Y<(sig.cutoff))

  cor_X = b_X / sqrt(b_X^2 + (n_X-2)*se_X^2)
  cor_Y = b_Y / sqrt(b_Y^2 + (n_Y-2)*se_Y^2)

  se_corX = sqrt((1-cor_X^2)^2/n_X)
  se_corY = sqrt((1-cor_Y^2)^2/n_Y)

  if(random.seed)
  {
    set.seed(random.seed)
  }

  # Without Screening
  CDcML_XtoY_NoS = mr_cML_DP(b_exp = cor_X[ind_X],
                             b_out = cor_Y[ind_X],
                             se_exp = se_corX[ind_X],
                             se_out = se_corY[ind_X],
                             n = min(n_X,n_Y),random_start = random_start,
                             random_start_pert = random_start,
                             num_pert = num_pert)
  set.seed(1)
  CDcML_YtoX_NoS = mr_cML_DP(b_exp = cor_Y[ind_Y],
                             b_out = cor_X[ind_Y],
                             se_exp = se_corY[ind_Y],
                             se_out = se_corX[ind_Y],
                             n = min(n_X,n_Y),random_start = random_start,
                             random_start_pert = random_start,
                             num_pert = num_pert)

  # With Screening
  intersect.ind.X.Y = intersect(ind_X,ind_Y)
  ind_X_new = setdiff(ind_X,
                      intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                          (abs(cor_Y)[intersect.ind.X.Y])])
  ind_Y_new = setdiff(ind_Y,
                      intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                          (abs(cor_Y)[intersect.ind.X.Y])])

  CDcML_XtoY_S = mr_cML_DP(b_exp = cor_X[ind_X_new],
                           b_out = cor_Y[ind_X_new],
                           se_exp = se_corX[ind_X_new],
                           se_out = se_corY[ind_X_new],
                           n = min(n_X,n_Y),random_start = random_start,
                           random_start_pert = random_start,
                           num_pert = num_pert)

  CDcML_YtoX_S = mr_cML_DP(b_exp = cor_Y[ind_Y_new],
                           b_out = cor_X[ind_Y_new],
                           se_exp = se_corY[ind_Y_new],
                           se_out = se_corX[ind_Y_new],
                           n = min(n_X,n_Y),random_start = random_start,
                           random_start_pert = random_start,
                           num_pert = num_pert)
  return(list(XtoY.est.NoS = CDcML_XtoY_NoS$MA_BIC_theta,
              XtoY.se.NoS = CDcML_XtoY_NoS$MA_BIC_se,
              YtoX.est.NoS = CDcML_YtoX_NoS$MA_BIC_theta,
              YtoX.se.NoS = CDcML_YtoX_NoS$MA_BIC_se,
              XtoY.est.S = CDcML_XtoY_S$MA_BIC_theta,
              XtoY.se.S = CDcML_XtoY_S$MA_BIC_se,
              YtoX.est.S = CDcML_YtoX_S$MA_BIC_theta,
              YtoX.se.S = CDcML_YtoX_S$MA_BIC_se,
              XtoY.est.NoS.DP = CDcML_XtoY_NoS$MA_BIC_DP_theta,
              XtoY.se.NoS.DP = CDcML_XtoY_NoS$MA_BIC_DP_se,
              YtoX.est.NoS.DP = CDcML_YtoX_NoS$MA_BIC_DP_theta,
              YtoX.se.NoS.DP = CDcML_YtoX_NoS$MA_BIC_DP_se,
              XtoY.est.S.DP = CDcML_XtoY_S$MA_BIC_DP_theta,
              XtoY.se.S.DP = CDcML_XtoY_S$MA_BIC_DP_se,
              YtoX.est.S.DP = CDcML_YtoX_S$MA_BIC_DP_theta,
              YtoX.se.S.DP = CDcML_YtoX_S$MA_BIC_DP_se))

}
