# parstab = matrix(NA, nrow = length(setList), ncol = 2)
# colnames(parstab) = c("shape","scale")
#
# gi = bgm = c()

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE

  fit.t<-fitdistrplus::fitdist(
    1/uE^2, "gamma", method = "mle",
    start = list(shape = 10, scale = 3),
    keepdata = FALSE)
  pars = summary(fit.t)$estimate
  shape = pars["shape"]
  scale = pars["scale"]

  X = log(uE^2)
  h1 = hist(X, nclass = 55, plot = FALSE)
  plot(h1, freq = FALSE,
       main = paste0('Set ',i), xlab = 'log(uE^2)'
  )
  grid()
  curve(
    MCMCpack::dinvgamma(exp(x),shape = shape, scale = 1/scale)*exp(x),
    from = min(X), to = max(X), lwd = 1, #3*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)
  box()
  # parstab[i,]=c(shape, 1/scale)
  #
  # gi[i] = ErrViewLib::gimc(uE)
  # bgm[i] = ErrViewLib::skewgm(uE)
}
