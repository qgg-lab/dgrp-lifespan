# make qq plot for a vector of p values
# ============================================================

# assume the null distribution is a uniform distribution
# thus the expected ith order p value is i/n (n = total number of p values)
# the distribution of the ith order p value is beta(i, n + 1 - i)
# see wiki http://en.wikipedia.org/wiki/Order_statistic

# created 07/23/2014 by Wen Huang

# ============================================================

pQQ <- function(pval, down = TRUE, down.rate = 0.1, n.cut = 1000, p.cut = 1e-4, plot.conf = FALSE, conf = c(0.025, 0.975), ...) {
  # if there are more than n.cut p values <= p.cut
  # only plot the first n.cut
  # extra arguments passed to base R plot
  n.pval <- length(pval)
  logp <- -log10(sort(pval))
  loge <- -log10((1:n.pval)/n.pval )
  p.cut <- -log10(p.cut)
  # make an effort to down-sample the p values to speed up plotting
  if (down & n.pval > n.cut) {
    if (logp[n.cut] < p.cut) {
      n.cut <- sum(logp >= p.cut)
    }
    down.sample <- c(1:n.cut, sort(sample((n.cut + 1):n.pval, size = floor((n.pval - n.cut)*down.rate), replace = FALSE)))
  } else {
    down.sample <- 1:n.pval
  }
  # plot diagnol line and confidence interval
  plot.max <- max(c(logp[1], loge[1]))
  if (plot.conf) {
    plot(c(0, plot.max), c(0, plot.max), type = "n", xlab = expression(Expected~~-log[10](italic(P))),
         ylab = expression(Observed~~-log[10](italic(P))), ...)
    lines(rev(loge[down.sample]), -log10(rev(qbeta(conf[1], down.sample, n.pval - down.sample + 1))), lty = 2)
    lines(rev(loge[down.sample]), -log10(rev(qbeta(conf[2], down.sample, n.pval - down.sample + 1))), lty = 2)
  } else {
    plot(c(0, plot.max), c(0, plot.max), type = "n", xlab = expression(Expected~~-log[10](italic(P))),
         ylab = expression(Observed~~-log[10](italic(P))), ...)
  }  
  abline(a = 0, b = 1, col = "red")
  points(loge[down.sample], logp[down.sample], pch = 19, cex = 0.6)
}
