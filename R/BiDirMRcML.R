BiDirMRcML <- function(b_X,b_Y,se_X,se_Y,n_X,n_Y,
                       sig.cutoff = 5e-8, num_pert = 100,
                       random_start = 0, random.seed = 0)
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
  MRcML_XtoY_NoS = mr_cML_DP(b_exp = b_X[ind_X],
                             b_out = b_Y[ind_X],
                             se_exp = se_X[ind_X],
                             se_out = se_Y[ind_X],
                             n = min(n_X,n_Y),random_start = random_start,
                             random_start_pert = random_start,
                             num_pert = num_pert)
  MRcML_YtoX_NoS = mr_cML_DP(b_exp = b_Y[ind_Y],
                             b_out = b_X[ind_Y],
                             se_exp = se_Y[ind_Y],
                             se_out = se_X[ind_Y],
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
  MRcML_XtoY_S = mr_cML_DP(b_exp = b_X[ind_X_new],
                           b_out = b_Y[ind_X_new],
                           se_exp = se_X[ind_X_new],
                           se_out = se_Y[ind_X_new],
                           n = min(n_X,n_Y),random_start = random_start,
                           random_start_pert = random_start,
                           num_pert = num_pert)
  MRcML_YtoX_S = mr_cML_DP(b_exp = b_Y[ind_Y_new],
                           b_out = b_X[ind_Y_new],
                           se_exp = se_Y[ind_Y_new],
                           se_out = se_X[ind_Y_new],
                           n = min(n_X,n_Y),random_start = random_start,
                           random_start_pert = random_start,
                           num_pert = num_pert)
  return(list(XtoY.est.NoS = MRcML_XtoY_NoS$MA_BIC_theta,
              XtoY.se.NoS = MRcML_XtoY_NoS$MA_BIC_se,
              YtoX.est.NoS = MRcML_YtoX_NoS$MA_BIC_theta,
              YtoX.se.NoS = MRcML_YtoX_NoS$MA_BIC_se,
              XtoY.est.S = MRcML_XtoY_S$MA_BIC_theta,
              XtoY.se.S = MRcML_XtoY_S$MA_BIC_se,
              YtoX.est.S = MRcML_YtoX_S$MA_BIC_theta,
              YtoX.se.S = MRcML_YtoX_S$MA_BIC_se,
              XtoY.est.NoS.DP = MRcML_XtoY_NoS$MA_BIC_DP_theta,
              XtoY.se.NoS.DP = MRcML_XtoY_NoS$MA_BIC_DP_se,
              YtoX.est.NoS.DP = MRcML_YtoX_NoS$MA_BIC_DP_theta,
              YtoX.se.NoS.DP = MRcML_YtoX_NoS$MA_BIC_DP_se,
              XtoY.est.S.DP = MRcML_XtoY_S$MA_BIC_DP_theta,
              XtoY.se.S.DP = MRcML_XtoY_S$MA_BIC_DP_se,
              YtoX.est.S.DP = MRcML_YtoX_S$MA_BIC_DP_theta,
              YtoX.se.S.DP = MRcML_YtoX_S$MA_BIC_DP_se))
}
