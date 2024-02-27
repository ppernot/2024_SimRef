ENCEfun = function(x, data, intrv) {
  # !!!! REQUIRES ORDERED uE !!!!!
  E  = data[x,1]
  uE = data[x,2]

  RMV = RMSE = c()
  for(i in 1:intrv$nbr) {
    sel     = intrv$lwindx[i]:intrv$upindx[i]
    RMV[i]  = sqrt(mean(uE[sel]^2))
    RMSE[i] = sqrt(mean(E[sel]^2))
  }
  ENCE  = mean( abs(RMV - RMSE) / RMV )
  # ENCEs = ENCE * sqrt( length(uE) / intrv$nbr )
  return(ENCE)
}
ZMSEfun = function(x, data, intrv) {
  # !!!! REQUIRES ORDERED uE !!!!!
  E  = data[x,1]
  uE = data[x,2]
  Z  = E / uE

  AE = c()
  for(i in 1:intrv$nbr) {
    sel   = intrv$lwindx[i]:intrv$upindx[i]
    AE[i] = abs(log(mean(Z[sel]^2)))
  }
  ZMSE = mean(AE)
  # ZMSEs = ZMSE * sqrt( length(uE) / intrv$nbr )
  return(ZMSE)
}

nMC     = 5000
nuSeq   = c(3:6,12,24)
mSeq    = c(2000,4000,8000,12000,16000)
nBinSeq = 10*(1:5)

# D = N(0,1) ####
smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rnorm(Mt)
        smct[kt,1:2] = ENCEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ENCEfun,
     file = 'testENCE.Rda')

smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rnorm(Mt)
        smct[kt,1:2] = ZMSEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ZMSEfun,
     file = 'testZMSE.Rda')


png(
  file = file.path(figDir, paste0('fig_ENCE_ZMSE.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,2),
  mar = c(4,4,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.0*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
lwd = 1.5 *gPars$lwd
pch = c(0,1,2,3,6)

load(file = 'testENCE.Rda')
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ENCE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ENCE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

load(file = 'testZMSE.Rda')
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ZMSE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ZMSE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

dev.off()


# D = ts(6) ####

smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rT4(Mt, df = 6)
        smct[kt,1:2] = ENCEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ENCEfun,
     file = 'testENCE_T6.Rda')

smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rT4(Mt, df = 6)
        smct[kt,1:2] = ZMSEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ZMSEfun,
     file = 'testZMSE_T6.Rda')


png(
  file = file.path(figDir, paste0('fig_ENCE_ZMSE_T6.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,2),
  mar = c(4,4,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.0*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
lwd = 1.5 *gPars$lwd
pch = c(0,1,2,3,6)

load(file = 'testENCE_T6.Rda')
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ENCE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ENCE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

load(file = 'testZMSE_T6.Rda')
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ZMSE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ZMSE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

dev.off()

# D = ts(12) ####

smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rT4(Mt, df = 12)
        smct[kt,1:2] = ENCEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ENCEfun,
     file = 'testENCE_T12.Rda')

smct =  matrix(NA, nrow = nMC, ncol = 2)
scorest = uscorest =
  array(
    NA,
    dim = c(length(mSeq), length(nuSeq), length(nBinSeq), 2),
    dimnames = list(mSeq, nuSeq, nBinSeq, 1:2)
  )
for(i in seq_along(mSeq)) {
  Mt = mSeq[i]
  for(j in seq_along(nuSeq)) {
    nut = nuSeq[j]
    for(k in seq_along(nBinSeq)) {
      nBint = nBinSeq[k]
      print(c(Mt,nut,nBint))
      intrvt = ErrViewLib::genIntervals(1:Mt, nBint)
      for(kt in 1:nMC) {
        uEt = sort(sqrt(MCMCpack::rinvgamma(Mt,nut,nut)))
        Et  = uEt * rT4(Mt, df = 12)
        smct[kt,1:2] = ZMSEfun(1:Mt, cbind(Et,uEt), intrvt)
      }
      scorest[i,j,k,]  = apply(smct,2,mean)
      uscorest[i,j,k,] = apply(smct,2,sd) / sqrt(nMC)
    }
  }
}
save(nMC, nuSeq, mSeq, nBinSeq, scorest, uscorest, ZMSEfun,
     file = 'testZMSE_T12.Rda')

# Checks ####

load(file = 'testENCE_T4.Rda')

lsize = log(mSeq)
lbin = log(nBinSeq)

sc = scorest[,'24',,1]
for(i in 1:ncol(sc)) {
  lscore = log(sc[,i])
  reg = lm(lscore~lsize)
  print(coefficients(reg))
}

lbin = log(nBinSeq)
for(i in 1:nrow(sc)) {
  lscore = log(sc[i,])
  reg = lm(lscore~lbin)
  print(coefficients(reg))
}


par(mfrow=c(1,2))
load(file = 'testENCE.Rda')
sc1 = scorest
load(file = 'testENCE_T4.Rda')
sc2 = scorest
load(file = 'testENCE_T12.Rda')
sc3 = scorest

x = y1 = y2 = y3 = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x  = c(x, sqrt(nBint/mSeq))
    y1 = c(y1, sc1[,k,j,1])
    y2 = c(y2, sc2[,k,j,1])
    y3 = c(y3, sc3[,k,j,1])
  }
}

lx  = log(x)
ly  = log(y1)
reg = lm(ly~lx)
co1 = coefficients(reg)
ly  = log(y2)
reg = lm(ly~lx)
co2 = coefficients(reg)
ly  = log(y3)
reg = lm(ly~lx)
co3 = coefficients(reg)

ylim = c(0,1.1*max(c(y1,y2,y3)))
plot(x,y1, pch=16, col=gPars$cols[1], log = '',
     xlim = c(0,1.1*max(x)), xaxs = 'i',
     xlab = '(M/N)^1/2',
     ylab = expression(tilde(theta)[list(nu,ref)]),
     ylim = ylim,
     main = 'ENCE')
grid(equilogs = FALSE)
points(x,y2,pch=17,col=gPars$cols[2])
points(x,y3,pch=18,col=gPars$cols[3])
curve(exp(co1[1])*x^co1[2],
      from = 0, to = max(x)*2,
      add = TRUE, lty = 2,col=gPars$cols[1])
curve(exp(co2[1])*x^co2[2],
      from = 0, to = max(x)*2,
      add = TRUE, lty = 2,col=gPars$cols[2])
curve(exp(co3[1])*x^co3[2],
      from = 0, to = max(x)*2,
      add = TRUE, lty = 2,col=gPars$cols[3])
legend(
  'topleft', bty = 'n', cex = 0.85,
  title = expression(beta),
  legend = round(c(co1[2],co2[2],co3[2]),2),
  lty = 2,
  col = gPars$cols[1:3],
  pch = NULL
)
box()

load(file = 'testZMSE.Rda')
sc1 = scorest
load(file = 'testZMSE_T4.Rda')
sc2 = scorest
load(file = 'testZMSE_T12.Rda')
sc3 = scorest

x = y1 = y2 = y3 = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x  = c(x, sqrt(nBint/mSeq))
    y1 = c(y1, sc1[,k,j,1])
    y2 = c(y2, sc2[,k,j,1])
    y3 = c(y3, sc3[,k,j,1])
  }
}

lx  = log(x)
ly  = log(y1)
reg = lm(ly~lx)
co1 = coefficients(reg)
ly  = log(y2)
reg = lm(ly~lx)
co2 = coefficients(reg)
ly  = log(y3)
reg = lm(ly~lx)
co3 = coefficients(reg)

ylim = c(0,1.1*max(c(y1,y2,y3)))
plot(x,y1, pch=16, col=gPars$cols[1], log = '',
     xlim = c(0,1.1*max(x)), xaxs = 'i',
     xlab = '(M/N)^1/2',
     ylab = expression(tilde(theta)[list(nu,ref)]),
     ylim = ylim,
     main = 'ZMSE')
grid(equilogs = FALSE)
points(x,y2,pch=17,col=gPars$cols[2])
points(x,y3,pch=18,col=gPars$cols[3])
curve(exp(co1[1])*x^co1[2],
      from = 0, to = max(x)*2,
      add = TRUE, lty = 2,col=gPars$cols[1])
curve(exp(co2[1])*x^co2[2],
      from = 0, to = max(x)*2,
      add = TRUE, lty = 2,col=gPars$cols[2])
curve(exp(co3[1])*x^co3[2],
      from = 0, to = max(x)*2,
      add = TRUE, lty = 2,col=gPars$cols[3])
legend(
  'topleft', bty = 'n', cex = 0.85,
  title = expression(beta),
  legend = round(c(co1[2],co2[2],co3[2]),2),
  lty = 2,
  col = gPars$cols[1:3],
  pch = NULL
)
box()

## Lit. Datasets ####

png(
  file = file.path(figDir, paste0('fig_ZMSEval.png')),
  width  = 2.5*gPars$reso,
  height = 2.5*gPars$reso
)
par(
  mfrow = c(3,3),
  mar = c(4,4,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = 2*gPars$lwd
)
lwd = 1.5 *gPars$lwd
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  io = order(D2$uE)
  uEt = D2$uE[io]
  Et  = D2$E[io]
  res = varBinSize(Et, uEt)
  Mt = length(uEt)

  x = sqrt(res$nBins / Mt)
  y = res$ZMSE
  plot(x, y,
       col = gPars$cols[5], pch = 16, cex = 1,
       xlab = '(N/M)^1/2', xlim = c(0,1.1*max(x)), xaxs = 'i',
       ylab = 'ZMSE', ylim = c(0,1.1*max(y)), yaxs = 'i',
       main = paste0('Set ',i))
  grid()

  sel = res$nBins > 20
  x   = x[sel]
  y   = y[sel]
  reg = lm(y~x)
  int = summary(reg)$coefficients[1,1]
  Uint = 2*summary(reg)$coefficients[1,2]
  abline(reg = reg, lty = 1, lwd = lwd, col = gPars$cols[2])
  segments(
    0.0005,int-Uint,
    0.0005,int+Uint,
    lwd = 4*lwd,
    lty = 1,
    lend = 2,
    col = gPars$cols[2])
  abline(a=0, b=1.14, lty = 2, lwd = lwd,
         col = gPars$cols[1])
  abline(a=0.006, b=1.577, lty = 3, lwd = lwd,
         col = gPars$cols[1])
  box()
  if(i==1)
    legend(
      'topleft', bty = 'n', cex = 0.65,
      legend = c('Data','Linear fit','NIG ref.','T6IG ref.'),
      lty = c(0,1,2,3),
      pch = c(16,NULL,NULL,NULL),
      col = gPars$cols[c(5,2,1,1)]
    )
}
dev.off()







for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  io = order(D2$uE)
  uEt = D2$uE[io]
  Et  = D2$E[io]
  res = varBinSize(Et, uEt)
  Mt = length(uEt)

  x = sqrt(res$nBins / Mt)
  y = res$ENCE
  plot(x, y, log = '',
       xlim = c(0,max(x)),
       ylim = c(0,max(y)),
       main = setList[i])
  grid()
  abline(reg=lm(y~x), lty = 2, col = 1)

  Et  = uEt * rnorm(Mt)
  res = varBinSize(Et, uEt)
  y = res$ENCE
  points(x, y, pch = 2, col = 2)
  grid()
  abline(reg=lm(y~x), lty =2 , col = 2)

  Et  = uEt * rT4(Mt, df = 6)
  res = varBinSize(Et, uEt)
  y = res$ENCE
  points(x, y, pch = 2, col = 2)
  grid()
  abline(reg=lm(y~x), lty =2 , col = 2)

}


# Mt  = 10000
# uEt = sort(sqrt(MCMCpack::rinvgamma(Mt, nut, nut)))
# Et  = uEt * rT4(Mt, df = 12)
# Mt  = length(uEt)
# res = varBinSize(Et, uEt)
# x = 1/sqrt(Mt/res$nBins)
# y = res$ENCE
# plot(x, y, log = '',
#      xlim = c(0,max(x)),
#      ylim = c(0,max(y)),
#      main = '')
# grid()
# co = coefficients(lm(y~x)); print(co)
# abline(reg=lm(y~x), lty =2, col = 2)




load(file = 'testENCE.Rda')
sc1 = scorest
load(file = 'testENCE_T4.Rda')
sc3 = scorest
load(file = 'testENCE_T12.Rda')
sc2 = scorest

x = y1 = y2 = y3 = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x  = c(x, mSeq/nBint)
    y1 = c(y1, sc1[,k,j,1])
    y2 = c(y2, sc2[,k,j,1])
    y3 = c(y3, sc3[,k,j,1])
  }
}
x = 1/sqrt(x)
plot(x, y1, log = '',
     xlim = c(0,max(x)),
     ylim = c(0,max(c(y1,y2,y3))),
     main = 'ENCE')
points(x,y2, pch = 2, col = 2)
points(x,y3, pch = 3, col = 3)
grid()
sel = x <= 0.05
x = x[sel]; y1 = y1[sel]; y2 = y2[sel]; y3 = y3[sel]
abline(reg=lm(y1~x), lty =2, col = 1)
abline(reg=lm(y2~x), lty =2, col = 2)
abline(reg=lm(y3~x), lty =2, col = 3)
