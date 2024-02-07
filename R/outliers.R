MVfun  = function(X) {
  uE = X[,2]
  mean(uE^2)
}
MSEfun  = function(X) {
  E = X[,1]
  mean(E^2)
}
RCEfun  = function(X) {
  E  = X[,1]
  uE = X[,2]
  (sqrt(mean(uE^2)) - sqrt(mean(E^2))) / sqrt(mean(uE^2))
}
ZMSfun  = function(X) {
  E  = X[,1]
  uE = X[,2]
  Z  = E / uE
  mean(Z^2)
}
calScoresBS3 = function(x, data) {
  E  = data[x,1]
  uE = data[x,2]
  Z  = E / uE
  c(
    RCE  = (sqrt(mean(uE^2)) - sqrt(mean(E^2))) / sqrt(mean(uE^2)),
    ZMS  = mean(Z^2)
  )
}
Sconf0 = function( X, stat, pcVec = 0:99) {
  M = NROW(X)
  S0 = stat(X)
  vstat = rep(0.0, length(pcVec))
  vstat[1] = S0
  for (i in 2:length(pcVec)) {
    k = pcVec[i]
    sel = 1:floor(k * M / 100)
    if (length(sel) == 0) {
      vstat[i] = NA
    } else {
      vstat[i] = stat(X[-sel,])
    }
  }
  return(vstat)
}

png(
  file = file.path(figDir, paste0('fig_stability0.png')),
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
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE
  E  = D2$E
  M = length(E)
  Z  = E / uE

  pcVec = seq(0,10,by=0.1)
  io = order(uE, E, decreasing = TRUE)
  X = cbind(E[io],uE[io])
  ccRCE = Sconf0(X , RCEfun, pcVec)
  ccRCE = ccRCE-ccRCE[1]
  ccZMS = Sconf0(X , ZMSfun, pcVec)
  ccZMS = ccZMS-ccZMS[1]

  bs = fPredBS(X, calScoresBS3, cl = cl)
  ciRCE = c(bs$bca[1,1]-bs$t0[1],bs$bca[2,1]-bs$t0[1])
  ciZMS = c(bs$bca[1,2]-bs$t0[2],bs$bca[2,2]-bs$t0[2])

  ylim = range(c(ccRCE,ccZMS,ciRCE,ciZMS))
  plot(
    pcVec, ccZMS, log = '', type = 'b',
    col = gPars$cols[5], lwd = 2*gPars$lwd,
    pch = 16, cex = 0.5,
    xlab = 'k% discarded', xlim = c(-1,10),
    ylab = 'ZMS, RCE', ylim = ylim,
    main = paste0('Set ', i)
  )
  grid(equilogs = FALSE)
  points(pcVec, ccRCE, col = gPars$cols[2], type = 'b',
         pch = 16, lwd = 2*gPars$lwd, cex = 0.5)
  abline(h=c(0,1), lty =2, col = gPars$cols[1])

  segments(-1, ciRCE[1], -1, ciRCE[2],
           lwd=40, lend=2, col=gPars$cols_tr2[2])
  segments(-0.5, ciZMS[1], -0.5, ciZMS[2],
           lwd=40, lend=2, col=gPars$cols_tr2[5])

  if(i==1)
    legend(
      'topright', bty = 'n',
      legend = c('ZMS','RCE'),
      pch = 16, lty = 0,
      col = gPars$cols[c(5,2)]
    )
}
stopCluster(cl)
dev.off()


png(
  file = file.path(figDir, paste0('fig_stabilityAux.png')),
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

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE
  E  = D2$E

  pcVec = seq(0, 10, by = 0.1)
  io = order(uE, decreasing = TRUE)
  X = cbind(E[io],uE[io])
  ccRCE = Sconf0(X , MSEfun, pcVec)
  # ccRCE = ccRCE-ccRCE[1]
  ccZMS = Sconf0(X , MVfun, pcVec)
  # ccZMS = ccZMS-ccZMS[1]

  ylim = range(c(ccRCE,ccZMS))
  plot(
    pcVec, ccZMS, log = '', type = 'b',
    col = gPars$cols[5], lwd = 2*gPars$lwd,
    pch = 16, cex = 0.5,
    xlab = 'k% discarded',
    ylab = 'MSE, MV', ylim = ylim,
    main = paste0('Set ', i)
  )
  grid(equilogs = FALSE)
  points(pcVec, ccRCE, col = gPars$cols[2], type = 'b',
         pch = 16, lwd = 2*gPars$lwd, cex = 0.5)
  abline(h=c(0,1), lty =2, col = gPars$cols[1])

  if(i==1)
    legend(
      'topright', bty = 'n',
      legend = c('MV','MSE'),
      pch = 16, lty = 0,
      col = gPars$cols[c(5,2)]
    )
}
dev.off()




M  = 5000
png(
  file = file.path(figDir, paste0('fig_stability.png')),
  width  = 2*gPars$reso,
  height = 1*gPars$reso
)
par(
  mfrow = c(1, 2),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)

ylim = c(-0.5, 1.5)

for(df in c(2,20)) {
  for(iMC in 1:100) {
    uE = rinvgamma(M, df, df)
    E  = uE * rnorm(M)
    Z = E / uE

    X2 = uE^2
    Y2 = E^2
    Z2 = Z^2
    io = order(X2)
    X2 = X2[io]
    Y2 = Y2[io]
    Z2 = Z2[io]
    cX2 = cumsum(X2) / seq_along(X2)
    cY2 = cumsum(Y2) / seq_along(Y2)
    cXY = (sqrt(cX2)-sqrt(cY2))/ sqrt(cX2)
    cZ2 = cumsum(Z2) / seq_along(Z2)

    if(iMC == 1) {
      plot(
        X2, cZ2, log = 'x',
        col = gPars$cols_tr[5], lwd = 1, pch = 16, cex = 0.5,
        xlab = 'uE^2', ylab = '<Z^2>, <E^2>/<uE^2>', ylim = ylim,
        main = paste0('nu = ', df)
      )
      grid()
    } else {
      points(X2, cZ2, col = gPars$cols_tr[5], pch = 16, cex = 0.5)
    }
    points(X2, cXY, col = gPars$cols_tr[2], pch = 16, cex = 0.5)
  }
  abline(h=c(0,1), lty =2, col = gPars$cols[1])
}
dev.off()
