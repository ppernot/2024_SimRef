cl <- makeCluster(detectCores())
methods = c('BS','SimN','SimT','SimRasN','SimRasT')

stats   = names(calScoresBS1(1:M,cbind(E,uE)))

muBS    = c(1, 0, NA)

smc     = matrix(NA, nrow = nMC, ncol = length(stats))

scores = bias = zmatBS  = zmatSimN = zmatSimT =
  zmatSim2N = zmatSim2T = muSimN = seSimN = muSimT = seSimT =
  ciScores = ciRasT = ciRasN =
  matrix(NA, nrow = length(setList), ncol = length(stats))

colnames(zmatBS) = colnames(zmatSimN) = colnames(zmatSimT) =
  colnames(zmatSim2N) = colnames(zmatSim2T) =stats

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  uE = D2$uE
  E  = D2$E
  M  = length(uE)
  mu = mean(E)

  # BS scores and CIs
  bs = fPredBS(cbind(E,uE), calScoresBS1, cl = cl)
  scores[i,] = bs$t0
  bias[i,]   = bs$bias
  ciScores[i,] = paste0(
    '[', signif(bs$bca[1,],3),', ',
    signif(bs$bca[2,],3), ']')

  # Known target
  zmatBS[i,] = fZetaBS(bs, muBS)

  # Normal simulation of target
  for(j in 1:nMC)
    smc[j,] = calScoresBS1(1:M, cbind(uE * rnorm(M), uE))
  muSimN[i,] = apply(smc, 2, mean, na.rm = TRUE)
  seSimN[i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
  zmatSimN[i,] = fZetaBS(bs, muSimN[i,],Utarget = 2*seSimN[i,])

  ciSim = apply(smc, 2, quantile, probs = c(0.025,0.975), na.rm = TRUE)
  ciRasN[i,] = paste0(
    '[', signif(ciSim[1,],3),', ',
    signif(ciSim[2,],3), ']')
  bsSim =  list()
  bsSim$t0 = muSimN[i,]
  bsSim$bias = 0
  bsSim$bca = ciSim
  zmatSim2sN[i,] = -fZetaBS(bsSim, scores[i,])

  df = ZdistPars[i,"df"]
  for(j in 1:nMC)
    smc[j,] = calScoresBS1(
      1:M,
      cbind(
        uE * rt_ls(M, df = df, mu = 0, sigma = 1) /
          sqrt(df/(df-2)),
        uE
      )
    )
  muSimT[i,] = apply(smc, 2, mean, na.rm = TRUE)
  seSimT[i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
  zmatSimT[i,] = fZetaBS(bs, muSimT[i,], Utarget = 2*seSimT[i,])

  ciSim = apply(smc, 2, quantile, probs = c(0.025,0.975), na.rm = TRUE)
  ciRasT[i,] = paste0(
    '[', signif(ciSim[1,],3),', ',
    signif(ciSim[2,],3), ']')
  bsSim =  list()
  bsSim$t0 = muSimT[i,]
  bsSim$bias = 0
  bsSim$bca = ciSim
  zmatSim2T[i,] = -fZetaBS(bsSim, scores[i,])

}
stopCluster(cl)

save(
  stats, methods, setList, scores, bias, ciScores, zmatBS,
  muSimN, seSimN, zmatSimN, ci2N, zmatSim2N,
  muSimT, seSimT, zmatSimT, ci2T, zmatSim2T,
  file = 'zetaScores.Rda'
)

