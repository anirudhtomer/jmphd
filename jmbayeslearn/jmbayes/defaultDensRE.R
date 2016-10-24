defaultDensRE = function (b, mu, D = NULL, invD = NULL, log = FALSE, prop = TRUE) {
  dmvnorm(b, mu = mu, Sigma = D, invSigma = invD, log = log,
          prop = prop)
}