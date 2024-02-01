
setListCal = setList

nMC  = 10000
nBin = 20

D2 = dataList[[paste0(setListCal[1],'_cal')]]
uE = D2$uE
E  = D2$E
M = length(uE)
intrv = ErrViewLib::genIntervals(M, 10)
stats = names(calScoresBS2(1:M,cbind(E,uE),intrv))

# Sensitivity to uE ####
smc     = matrix(NA, nrow = nMC, ncol = length(stats))
scores1  = uscores1 = matrix(NA, nrow = length(setListCal), ncol = length(stats))
colnames(scores1) = colnames(uscores1) = stats

for(i in seq_along(setListCal)) {
  D2 = dataList[[paste0(setListCal[i],'_cal')]]
  print(setListCal[i])
  uE = D2$uE
  M  = length(uE)
  intrv = ErrViewLib::genIntervals(M, nBin)

  for(k in 1:nMC) {
    Ep = uE * rnorm(M)
    smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
  }
  scores1[i,] = apply(smc, 2, mean, na.rm = TRUE)
  uscores1[i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
}

png(
  file = file.path(figDir, paste0('fig_sensitivity1.png')),
  width  = 2.75*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,3),
  mar = c(4,3,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
X = 1:nrow(scores1)
for(stat in stats) {
  ylim = range(c(scores1[,stat] - 2*uscores1[,stat],scores1[,stat] + 2*uscores1[,stat]))
  plot(X, scores1[,stat],
       type = 'b', col = gPars$cols[5], pch = 1,
       lwd = gPars$lwd, lty = 1,
       xlab = 'Set #',
       ylab = stat,
       ylim = ylim,
       main =''
  )
  grid()
  segments(X, scores1[,stat] - 2*uscores1[,stat],
           X, scores1[,stat] + 2*uscores1[,stat],
           col = gPars$cols[5],
           lwd = gPars$lwd
  )

}
dev.off()

stop()

# Sensitivity to D ####
nuSeq = c(3,4,5,10,15,20)
smc     = matrix(NA, nrow = nMC, ncol = length(stats))
scores =
  array(
    NA,
    dim = c(length(nuSeq), length(setListCal), length(stats)),
    dimnames = list(nuSeq,paste0('Set',1:length(setListCal)),stats)
  )

nBin = 20

for(i in seq_along(setListCal)) {
  D2 = dataList[[paste0(setListCal[i],'_cal')]]
  print(setListCal[i])
  uE = D2$uE
  # E  = uE * rnorm(uE) #D2$E
  M  = length(uE)
  intrv = ErrViewLib::genIntervals(M, nBin)

  for(j in seq_along(nuSeq)) {
    df = nuSeq[j]; print(df)
    for(k in 1:nMC) {
      Ep = uE * rt_ls(M, df = df, mu=0, sigma=1) / sqrt(df/(df-2))
      smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
    }
    muSim = apply(smc, 2, mean, na.rm = TRUE)
    scores[j,i,] = muSim
  }
}

png(
  file = file.path(figDir, paste0('fig_sensitivity.png')),
  width  = 2.75*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,3),
  mar = c(4,3,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.5*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)

for(stat in stats) {
  matplot( nuSeq, scores[,,stat]/scores[length(nuSeq,,stat)],
           type = 'b', col = gPars$cols[1:7],
           lwd = gPars$lwd, lty = 1,
           xlab = expression(nu),
           ylab = stat,
           main =''
  )
  grid()
}
# if(iset == 1)
#   legend(
#     'topleft', bty = 'n',
#     legend = methods,
#     pch = 17:20, col = gPars$cols[2:5], lty = 1
#   )
dev.off()
