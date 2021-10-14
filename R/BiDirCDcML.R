BiDirCDcML <- function(b_X,b_Y,se_X,se_Y,n_X,n_Y,
                       sig.cutoff = 5e-8, num_pert = 100,
                       random_start = 0, random.seed = 0
                       )
{
  pvalue.X = pnorm(-abs(b_X/se_X))*2
  pvalue.Y = pnorm(-abs(b_Y/se_Y))*2

  ind_X = which(pvalue.X<(5e-8))
  ind_Y = which(pvalue.Y<(5e-8))

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


}
