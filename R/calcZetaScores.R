# Sensitivity of scores to D

cl <- makeCluster(detectCores())
methods = c('BS','SimN','SimT','Sim2N','Sim2T')

nBin  = 20
intrv = ErrViewLib::genIntervals(1:100, nBin)
stats   = names(calScoresBS2(1:M, cbind(E,uE), intrv = intrv))

muBS    = c(1, NA, NA, NA)

smc     = matrix(NA, nrow = nMC, ncol = length(stats))

scores = bias = zmatBS  = zmatSimN = zmatSimT =
  zmatSim2N = zmatSim2T = muSimN = seSimN = muSimT = seSimT =
  ciScores = ci2T = ci2N =
  matrix(NA, nrow = length(setList), ncol = length(stats))

colnames(zmatBS) = colnames(zmatSimN) = colnames(zmatSimT) =
  colnames(zmatSim2N) = colnames(zmatSim2T) =stats


for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  uE = D2$uE
  E  = D2$E
  M  = length(uE)
  mu = mean(E)

  intrv = ErrViewLib::genIntervals(1:M, nBin)
  intrvJack = ErrViewLib::genIntervals(1:(M-1), nBin)

  # BS scores and CIs
  bs = fPredBS(cbind(E,uE), calScoresBS2, cl = cl,
               intrv = intrv, intrvJack = intrvJack)
  scores[i,] = bs$t0
  bias[i,]   = bs$bias
  ciScores[i,] = paste0(
    '[', signif(bs$bca[1,],3),', ',
    signif(bs$bca[2,],3), ']')

  # Known target
  zmatBS[i,] = fZetaBS(bs, muBS)

  # Normal simulation of target
  for(j in 1:nMC)
    smc[j,] = calScoresBS2(1:M, cbind(uE * rnorm(M), uE),
                           intrv = intrv, intrvJack = intrvJack)
  muSimN[i,] = apply(smc, 2, mean, na.rm = TRUE)
  seSimN[i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
  zmatSimN[i,] = fZetaBS(bs, muSimN[i,],Utarget = 2*seSimN[i,])

  ciSim = apply(smc, 2, quantile, probs = c(0.025,0.975),
                na.rm = TRUE)
  ci2N[i,] = paste0(
    '[', signif(ciSim[1,],3),', ',
    signif(ciSim[2,],3), ']')
  bsSim =  list()
  bsSim$t0 = muSimN[i,]
  bsSim$bias = 0
  bsSim$bca = ciSim
  zmatSim2N[i,] = -fZetaBS(bsSim, scores[i,])

  df = 6 #ZdistPars[i,"df"]
  for(j in 1:nMC)
    smc[j,] = calScoresBS2(
      1:M,
      cbind(
        uE * rt_ls(M, df = df, mu = 0, sigma = 1) /
          sqrt(df/(df-2)),
        uE
      ),
      intrv = intrv, intrvJack = intrvJack
    )
  muSimT[i,] = apply(smc, 2, mean, na.rm = TRUE)
  seSimT[i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
  zmatSimT[i,] = fZetaBS(bs, muSimT[i,], Utarget = 2*seSimT[i,])

  ciSim = apply(smc, 2, quantile, probs = c(0.025,0.975), na.rm = TRUE)
  ci2T[i,] = paste0(
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
  file = file.path(tmpDir,"zetaScores.Rda"))

