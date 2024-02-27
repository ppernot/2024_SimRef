CCfun  = function(X) {
  E  = X[,1]
  uE = X[,2]
  cor(abs(E), uE, method='pearson')
}
ZMSfun  = function(X) {
  E  = X[,1]
  uE = X[,2]
  Z  = E / uE
  mean(Z^2)
}
png(
  file = file.path(figDir, paste0('fig_CC_outl.png')),
  width  = 3*gPars$reso,
  height = 3*gPars$reso
)
par(
  mfrow = c(3,3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.5*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
cl <- makeCluster(detectCores())
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  uE = D2$uE
  E  = D2$E
  M = length(E)
  Z  = E / uE

  pcVec = seq(0,10,by=0.1)
  io    = order(uE, E, decreasing = TRUE)
  X     = cbind(E[io],uE[io])
  ccCC  = Sconf0(X , CCfun, pcVec)
  ccCC = ccCC-ccCC[1]
  ccZMS = Sconf0(X , ZMSfun, pcVec)
  ccZMS = ccZMS-ccZMS[1]
  bs    = fPredBS(X, calScoresBS1, cl = cl)
  ciZMS = c(bs$bca[1,1]-bs$t0[1],bs$bca[2,1]-bs$t0[1])
  ciCC  = c(bs$bca[1,2]-bs$t0[2],bs$bca[2,2]-bs$t0[2])

  ylim = range(c(ccCC,ccZMS,ciCC,ciZMS))
  plot(
    pcVec, ccZMS, log = '', type = 'b',
    col = gPars$cols[5], lwd = 2*gPars$lwd,
    pch = 16, cex = 0.5,
    xlab = 'k% discarded', xlim = c(-1,10),
    ylab = 'ZMS, RCE', ylim = ylim,
    main = paste0('Set ', i)
  )
  grid(equilogs = FALSE)
  points(pcVec, ccCC, col = gPars$cols[2], type = 'b',
         pch = 16, lwd = 2*gPars$lwd, cex = 0.5)
  abline(h=c(0,1), lty =2, col = gPars$cols[1])

  segments(-1, ciCC[1], -1, ciCC[2],
           lwd=40, lend=2, col=gPars$cols_tr2[2])
  segments(-0.5, ciZMS[1], -0.5, ciZMS[2],
           lwd=40, lend=2, col=gPars$cols_tr2[5])

  if(i==1)
    legend(
      'topright', bty = 'n',
      legend = c('ZMS','CC'),
      pch = 16, lty = 0,
      col = gPars$cols[c(5,2)]
    )
}
stopCluster(cl)
dev.off()

