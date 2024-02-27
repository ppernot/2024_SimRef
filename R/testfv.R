# nBin = 10
#
# df = 20
# Ep = uE * rT4(M, df = df)
# ErrViewLib::plotEvsPU(uE,Ep/uE,type='horiz',runQuant = TRUE)
# res = ErrViewLib::plotLZMS(
#   uE, Ep/uE, nBin = nBin, popMin = 30,
#   score = TRUE,plot = TRUE, parallel = TRUE)
# print(res)
#
# df = 3
# Ep = uE * rT4(M, df = df)
# ErrViewLib::plotEvsPU(uE,Ep/uE,type='horiz',runQuant = TRUE)
# res = ErrViewLib::plotLZMS(
#   uE, Ep/uE, nBin = nBin, popMin = 30,
#   score = TRUE,plot = TRUE, parallel = TRUE)
# print(res)


nBinSeq = seq(20, 100, by = 20)
fvZMSval =
  array(
    NA,
    dim = c(length(nBinSeq), length(setList), 3),
    dimnames = list(nBinSeq,paste0('Set',1:length(setList)),1:3)
  )
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  print(setList[i])
  E  = D2$E
  uE = D2$uE
  M  = length(uE)

  for(j in seq_along(nBinSeq)) {
    nBin = nBinSeq[j]; print(nBin)
    res = ErrViewLib::plotLZMS(
      uE, E/uE, nBin = nBin, popMin = 20,
      score = TRUE, plot = FALSE, parallel = TRUE)
    fvZMSval[j,i,]=c(res$fVal, res$lofVal, res$upfVal)
  }

}

png(
  file = file.path(figDir, paste0('fig_valid_fv.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(3, 3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1 * gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
for(i in seq_along(setList)) {
  ylim = c(min(fvZMSval[,i,]),1)
  plot(nBinSeq, fvZMSval[,i,1], type = 'b',
       pch = 16, col = gPars$cols[5],
       xlab = expression(N[b]),
       ylim = ylim, ylab = expression(f[list(v,ZMS)]),
       main = paste0('Set ',i))

  grid()
  abline(h=0.95, lty = 2, col = gPars$cols[1])
  abline(v=sqrt(length(dataList[[paste0(setList[i],'_cal')]]$E)), lty = 2,
         col = gPars$cols[2])
  for(j in seq_along(nBinSeq))
  arrows(nBinSeq[j],fvZMSval[j,i,2],nBinSeq[j],fvZMSval[j,i,3],
         angle = 90, lwd = 2* gPars$lwd, col = gPars$cols[5],
         code = 3)
  box()
}
dev.off()
