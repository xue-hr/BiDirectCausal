#' Bi-directional CD Methods
#'
#' @param b_X Vector of estimated effect sizes for trait 1.
#' @param b_Y Vector of estimated effect sizes for trait 2.
#' @param se_X Vector of standard errors of b_X.
#' @param se_Y Vector of standard errors of b_Y.
#' @param n_X Sample size for trait 1.
#' @param n_Y Sample size for trait 2.
#' @param sig.cutoff Significant cutoff to choose IV, defaul is 5e-8.
#' @param random.seed Random seed, default is 0, i.e. no random seed is set.
#'
#' @return List of estimated correlation ratio (est) and
#' corresponding standard error (se) obtained with CD-Ratio and CD-Egger, for
#' both directions. "NoS" indicates no screening, "S" indicates with screening.
#' @export
#'
#' @examples
BiDirCDMethod <- function(b_X,b_Y,se_X,se_Y,n_X,n_Y,
                          sig.cutoff = 5e-8, random.seed = 0)
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
  CD_XtoY_NoS = CDMethods(b_exp = b_X[ind_X],
                          b_out = b_Y[ind_X],
                          se_exp = se_X[ind_X],
                          se_out = se_Y[ind_X],
                          n_exp = n_X,n_out = n_Y)
  CD_YtoX_NoS = CDMethods(b_exp = b_Y[ind_Y],
                          b_out = b_X[ind_Y],
                          se_exp = se_Y[ind_Y],
                          se_out = se_X[ind_Y],
                          n_exp = n_Y,n_out = n_X)
  # With Screening
  intersect.ind.X.Y = intersect(ind_X,ind_Y)
  ind_X_new = setdiff(ind_X,
                      intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                          (abs(cor_Y)[intersect.ind.X.Y])])
  ind_Y_new = setdiff(ind_Y,
                      intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                          (abs(cor_Y)[intersect.ind.X.Y])])
  CD_XtoY_S = CDMethods(b_exp = b_X[ind_X_new],
                        b_out = b_Y[ind_X_new],
                        se_exp = se_X[ind_X_new],
                        se_out = se_Y[ind_X_new],
                        n_exp = n_X,n_out = n_Y)
  CD_YtoX_S = CDMethods(b_exp = b_Y[ind_Y_new],
                        b_out = b_X[ind_Y_new],
                        se_exp = se_Y[ind_Y_new],
                        se_out = se_X[ind_Y_new],
                        n_exp = n_Y,n_out = n_X)

  return(list(CDRatio.XtoY.est.NoS = CD_XtoY_NoS$CD_Ratio_result$T1toT2[1],
              CDRatio.XtoY.se.NoS = CD_XtoY_NoS$CD_Ratio_result$T1toT2[2],
              CDEgger.XtoY.est.NoS = CD_XtoY_NoS$CD_Egger_result$T1toT2[2],
              CDEgger.XtoY.se.NoS = CD_XtoY_NoS$CD_Egger_result$T1toT2[4],
              CDRatio.YtoX.est.NoS = CD_YtoX_NoS$CD_Ratio_result$T1toT2[1],
              CDRatio.YtoX.se.NoS = CD_YtoX_NoS$CD_Ratio_result$T1toT2[2],
              CDEgger.YtoX.est.NoS = CD_YtoX_NoS$CD_Egger_result$T1toT2[2],
              CDEgger.YtoX.se.NoS = CD_YtoX_NoS$CD_Egger_result$T1toT2[4],
              CDRatio.XtoY.est.S = CD_XtoY_S$CD_Ratio_result$T1toT2[1],
              CDRatio.XtoY.se.S = CD_XtoY_S$CD_Ratio_result$T1toT2[2],
              CDEgger.XtoY.est.S = CD_XtoY_S$CD_Egger_result$T1toT2[2],
              CDEgger.XtoY.se.S = CD_XtoY_S$CD_Egger_result$T1toT2[4],
              CDRatio.YtoX.est.S = CD_YtoX_S$CD_Ratio_result$T1toT2[1],
              CDRatio.YtoX.se.S = CD_YtoX_S$CD_Ratio_result$T1toT2[2],
              CDEgger.YtoX.est.S = CD_YtoX_S$CD_Egger_result$T1toT2[2],
              CDEgger.YtoX.se.S = CD_YtoX_S$CD_Egger_result$T1toT2[4]))


}
