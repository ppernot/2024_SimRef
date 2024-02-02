
setListCal = setList

nMC  = 10000
nBin = 20

# Get stats names
D2 = dataList[[paste0(setListCal[1],'_cal')]]
uE = D2$uE
E  = D2$E
M = length(uE)
intrv = ErrViewLib::genIntervals(M, nBin)
stats = names(calScoresBS2(1:M,cbind(E,uE),intrv))


# Synthetic datasets ####

## Sensitivity to uE ####
nuSeq = c(3,4,5,10,15,20)
M = 10000
nBin = 100
intrv = ErrViewLib::genIntervals(M, nBin)

smc     = matrix(NA, nrow = nMC, ncol = length(stats))
scores2  = uscores2 = matrix(NA, nrow = length(nuSeq), ncol = length(stats))
colnames(scores2) = colnames(uscores2) = stats
fvZMSscores2 = fvRCEscores2 = matrix(NA, nrow = length(nuSeq), ncol = 3)

for(i in seq_along(nuSeq)) {
  nu = nuSeq[i]; print(nu)
  uE = sqrt(MCMCpack::rinvgamma(M, nu/2, nu/2 ))
  for(k in 1:nMC) {
    Ep = uE * rnorm(M)
    smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
  }
  scores2[i,] = apply(smc, 2, mean, na.rm = TRUE)
  uscores2[i,] = apply(smc, 2, sd, na.rm = TRUE) / sqrt(nMC)

  res = ErrViewLib::plotLZMS(uE, Ep/uE, nBin = nBin, score = TRUE,
                             plot = FALSE, parallel = TRUE)
  fvZMSscores2[i,]=c(res$fVal, res$lofVal, res$upfVal)
  res = ErrViewLib::plotLRCE(uE, uE, Ep, nBin = nBin,
                             plot = FALSE, parallel = TRUE)
  fvRCEscores2[i,]=c(res$fVal, res$lofVal, res$upfVal)
}



X = nuSeq
sc = scores2
usc = uscores2
sclw = sc - 2*usc
scup = sc + 2*usc
png(
  file = file.path(figDir, paste0('fig_sens_Synth_uE.png')),
  width  = 2.75*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,3),
  mar = c(4,3,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
for(stat in stats) {
  ylim = range(c(sclw[,stat],scup[,stat]))
  plot(X, sc[,stat], log = 'x',
       type = 'b', col = gPars$cols[5], pch = 1,
       lwd = 2*gPars$lwd, lty = 1,
       xlab = expression(nu),
       ylab = stat,
       ylim = ylim,
       main =''
  )
  grid()
  abline(h = c(0,1), lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
  if(stat %in% c('ENCE','ZMSE'))
    abline(h = mean(sc[,stat]), lwd = gPars$lwd, lty = 2,
           col = gPars$cols[1] )
  segments(X, sclw[,stat],
           X, scup[,stat],
           col = gPars$cols[5],
           lwd = 2*gPars$lwd
  )
}
plot(X, fvZMSscores2[,1], log = 'x',
     type = 'b', col = gPars$cols[5], pch = 1,
     lwd = 2*gPars$lwd, lty = 1,
     xlab = expression(nu),
     ylab = expression(f[v]),
     ylim = c(0.8,1),
     main =''
)
grid()
points(X, fvRCEscores2[,1],
     type = 'b', col = gPars$cols[2], pch = 1,
     lwd = 2*gPars$lwd, lty = 1
)
abline(h = 0.95, lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
segments(X, fvZMSscores2[,2],
         X, fvZMSscores2[,3],
         col = gPars$cols[5],
         lwd = 2*gPars$lwd
)
segments(X, fvRCEscores2[,2],
         X, fvRCEscores2[,3],
         col = gPars$cols[2],
         lwd = 2*gPars$lwd
)
box()

legend(
  'bottomright', bty = 'n',
  legend = c("RCE","ZMS"),
  col = gPars$cols[c(2,5)],
  pch = 1, lty = 1, lwd = 2*gPars$lwd
)
dev.off()


## Sensitivity to D ####
set.seed(123)
nBin = 50
smc     = matrix(NA, nrow = nMC, ncol = length(stats))
scores3  = uscores3 = matrix(NA, nrow = length(nuSeq), ncol = length(stats))
colnames(scores3) = colnames(uscores3) = stats
fvZMSscores3 = fvRCEscores3 = matrix(NA, nrow = length(nuSeq), ncol = 3)

uE = sqrt(MCMCpack::rinvgamma(M, 2, 2 ))

for(i in seq_along(nuSeq)) {
  df = nuSeq[i]; print(df)
  for(k in 1:nMC) {
    Ep = uE * rt_ls(M, df = df, mu=0, sigma=1) / sqrt(df/(df-2))
    smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
  }
  scores3[i,] = apply(smc, 2, mean, na.rm = TRUE)
  uscores3[i,] = apply(smc, 2, sd, na.rm = TRUE) / sqrt(nMC)

  res = ErrViewLib::plotLZMS(uE, Ep/uE, nBin = nBin, score = TRUE,
                             plot = FALSE, parallel = TRUE)
  fvZMSscores3[i,]=c(res$fVal, res$lofVal, res$upfVal)
  res = ErrViewLib::plotLRCE(uE, uE, Ep, nBin = nBin,
                             plot = FALSE, parallel = TRUE)
  fvRCEscores3[i,]=c(res$fVal, res$lofVal, res$upfVal)
}

X = nuSeq
sc = scores3
usc = uscores3
sclw = sc - 2*usc
scup = sc + 2*usc
png(
  file = file.path(figDir, paste0('fig_sens_Synth_D.png')),
  width  = 2.75*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,3),
  mar = c(4,3,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
for(stat in stats) {
  ylim = range(c(sclw[,stat],scup[,stat]))
  plot(X, sc[,stat], log = 'x',
       type = 'b', col = gPars$cols[5], pch = 1,
       lwd = 2*gPars$lwd, lty = 1,
       xlab = expression(nu),
       ylab = stat,
       ylim = ylim,
       main =''
  )
  grid()
  abline(h = c(0,1), lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
  segments(X, sclw[,stat],
           X, scup[,stat],
           col = gPars$cols[5],
           lwd = 2*gPars$lwd
  )
}
plot(X, fvZMSscores3[,1], log = 'x',
     type = 'b', col = gPars$cols[5], pch = 1,
     lwd = 2*gPars$lwd, lty = 1,
     xlab = expression(nu),
     ylab = expression(f[v]),
     ylim = c(0.8,1),
     main =''
)
grid()
points(X, fvRCEscores3[,1],
       type = 'b', col = gPars$cols[2], pch = 1,
       lwd = 2*gPars$lwd, lty = 1
)
abline(h = 0.95, lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
segments(X, fvZMSscores3[,2],
         X, fvZMSscores3[,3],
         col = gPars$cols[5],
         lwd = 2*gPars$lwd
)
segments(X, fvRCEscores3[,2],
         X, fvRCEscores3[,3],
         col = gPars$cols[2],
         lwd = 2*gPars$lwd
)
box()

legend(
  'bottomright', bty = 'n',
  legend = c("RCE","ZMS"),
  col = gPars$cols[c(2,5)],
  pch = 1, lty = 1, lwd = 2*gPars$lwd
)
dev.off()



# Real datsets ####

# ## Sensitivity to uE ####
# smc     = matrix(NA, nrow = nMC, ncol = length(stats))
# scores1  = uscores1 = matrix(NA, nrow = length(setListCal), ncol = length(stats))
# colnames(scores1) = colnames(uscores1) = stats
#
# for(i in seq_along(setListCal)) {
#   D2 = dataList[[paste0(setListCal[i],'_cal')]]
#   print(setListCal[i])
#   uE = D2$uE
#   M  = length(uE)
#   intrv = ErrViewLib::genIntervals(M, nBin)
#
#   for(k in 1:nMC) {
#     Ep = uE * rnorm(M)
#     smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
#   }
#   scores1[i,] = apply(smc, 2, mean, na.rm = TRUE)
#   uscores1[i,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)
# }
#
# X = 1:length(setListCal)
# sc = scores1
# usc = uscores1
# sclw = sc - 2*usc
# scup = sc + 2*usc
# png(
#   file = file.path(figDir, paste0('fig_sens_Real_uE.png')),
#   width  = 2.75*gPars$reso,
#   height = 2*gPars$reso
# )
# par(
#   mfrow = c(2,3),
#   mar = c(4,3,1,1),
#   mgp = gPars$mgp,
#   pty = 's',
#   tcl = gPars$tcl,
#   cex = 1.25*gPars$cex,
#   cex.main = 1,
#   lwd = gPars$lwd
# )
# for(stat in stats) {
#   ylim = range(c(sclw[,stat],scup[,stat]))
#   plot(X, sc[,stat], log = '',
#        type = 'b', col = gPars$cols[5], pch = 1,
#        lwd = 2*gPars$lwd, lty = 1,
#        xlab = 'Set #',
#        ylab = stat,
#        ylim = ylim,
#        main =''
#   )
#   grid()
#   abline(h = c(0,1), lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
#   segments(X, sclw[,stat],
#            X, scup[,stat],
#            col = gPars$cols[5],
#            lwd = 2*gPars$lwd
#   )
# }
# dev.off()


## Sensitivity to D ####
nuSeq = c(3,4,5,10,15,20)
smc     = matrix(NA, nrow = nMC, ncol = length(stats))
scores = uscores =
  array(
    NA,
    dim = c(length(nuSeq), length(setListCal), length(stats)),
    dimnames = list(nuSeq,paste0('Set',1:length(setListCal)),stats)
  )
fvZMSscores = fvRCEscores =
  array(
    NA,
    dim = c(length(nuSeq), length(setListCal), 3),
    dimnames = list(nuSeq,paste0('Set',1:length(setListCal)),1:3)
  )

for(i in seq_along(setListCal)) {
  D2 = dataList[[paste0(setListCal[i],'_cal')]]
  print(setListCal[i])
  uE = D2$uE
  # E  = uE * rnorm(uE) #D2$E
  M  = length(uE)
  nBin  = 50
  intrv = ErrViewLib::genIntervals(M, nBin)

  for(j in seq_along(nuSeq)) {
    df = nuSeq[j]; print(df)
    for(k in 1:nMC) {
      Ep = uE * rt_ls(M, df = df, mu=0, sigma=1) / sqrt(df/(df-2))
      smc[k,] = calScoresBS2(1:M, cbind(Ep,uE), intrv)
    }
    scores[j,i,] = apply(smc, 2, mean, na.rm = TRUE)
    uscores[j,i,] = apply(smc, 2, sd, na.rm = TRUE)

    res = ErrViewLib::plotLZMS(uE, Ep/uE, nBin = nBin, popMin = 30,
                               score = TRUE,plot = FALSE, parallel = TRUE)
    fvZMSscores[j,i,]=c(res$fVal, res$lofVal, res$upfVal)
    # res = ErrViewLib::plotLRCE(uE, uE, Ep, nBin = nBin,
    #                            plot = FALSE, parallel = TRUE)
    # fvRCEscores[j,i,]=c(res$fVal, res$lofVal, res$upfVal)
  }
}
save(nuSeq, stats, setListCal, scores, uscores,
     fvZMSscores,fvRCEscores, file ='sensitivity_new.Rda')

# load(file ='sensitivity.Rda')
png(
  file = file.path(figDir, paste0('fig_sens_Real_D.png')),
  width  = 2.75*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,3),
  mar = c(4,3,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)

for(stat in stats) {
  if(stat == 'RCE') # Do not normalize by 0
    sc = t(t(scores[,,stat])-scores[length(nuSeq),,stat])
  else
    sc = t(t(scores[,,stat])/scores[length(nuSeq),,stat])
  matplot( nuSeq, sc, log = 'x',
           type = 'b', col = gPars$cols[1:7],
           lwd = 2* gPars$lwd, lty = 1,
           xlab = expression(nu),
           ylab = stat,
           main =''
  )
  grid()
  abline(h = c(0,1), lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
  box()
}
matplot(X, t(t(fvZMSscores[,,1])-fvZMSscores[length(nuSeq),,1]),
      log = 'x',type = 'b', col = gPars$cols[1:7],
     lwd = 2*gPars$lwd, lty = 1,
     xlab = expression(nu),
     ylab = expression(f[v]),
     # ylim = c(0.7,1),
     main =''
)
grid()
abline(h = 0, lwd = gPars$lwd, lty = 2, col = gPars$cols[1] )
# for(i in seq_along(setListCal)) {
#   segments(X, fvZMSscores[,i,2],
#            X, fvZMSscores[,i,3],
#            col = gPars$cols[1:7],
#            lwd = 2*gPars$lwd
#   )
# }
box()
dev.off()
