nMC = 10000
nuSeq = 3:10

M1 = 5000
M2 = 10000

tabSig2_1 = tabSig2_2 = matrix(NA, ncol = length(nuSeq), nrow = nMC)
tabExact = c()
for(j in seq_along(nuSeq)) {
  nu = nuSeq[j]
  for(i in 1:nMC) {
    sig2 = MCMCpack::rinvgamma(M1, nu/2, nu/2 )
    tabSig2_1[i,j] = mean(sig2)
    sig2 = MCMCpack::rinvgamma(M2, nu/2, nu/2 )
    tabSig2_2[i,j] = mean(sig2)
  }
  tabExact[j] = nu/(nu-2)
}

data1 = t(tabSig2_1) - tabExact
data1 = as.data.frame(t(data1))
names(data1) = nuSeq
data2 = t(tabSig2_2) - tabExact
data2 = as.data.frame(t(data2))
names(data2) = nuSeq

png(
  file = file.path(figDir, 'fig_InvGammaMean.png'),
  width  = 1.5 * gPars$reso,
  height = 1* gPars$reso
)
par(
  mfrow = c(1,2),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  lwd = gPars$lwd
)

boxplot(data1,
        col = gPars$cols[5],
        cex = 0.75, pch = '+',
        xlab = expression(nu),
        ylab = '<x^2> - exact',
        ylim = c(-0.6,2),
        main = paste0('M = ',M1)
)
grid()
abline(h=0, lty = 2)
box()

boxplot(data2,
        col = gPars$cols[5],
        cex = 0.75, pch = '+',
        xlab = expression(nu),
        ylab = '<x^2> - exact',
        ylim = c(-0.6,2),
        main = paste0('M = ',M2)
)
grid()
abline(h=0, lty = 2)
box()

dev.off()


