# ============================================
# = make figure to show overall distribution =
# ============================================

args <- commandArgs(TRUE) # args <- c("reportData/xqtl.eff.txt", "report/figure_xqtl.eff.size.pdf")
library("RColorBrewer")

# read data
# ============================================================

xqtl.gsi.eff <- read.csv(args[2], header = TRUE, as.is = TRUE)
xqtl.gei.eff <- read.csv(args[3], header = TRUE, as.is = TRUE)

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[4], width = file.width/25.4, height = file.width*0.5/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(2.4, 2.5, 1.1, 0.5), ps = 7, lwd = 0.5, mfrow = c(1, 2))


# plot data
# ============================================================

plot(xqtl.gsi.eff[, 1], xqtl.gsi.eff[, 2], xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE, pch = 21, bg = "grey", main = "", xlab = "", ylab = "", cex = 0.5)

axis(side = 1, at = seq(-1, 1, 0.5), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

abline(a = 0, b = 1, lty = 2, lwd = 0.5)

title(xlab = "Allele frequency difference\n(Old - Young) in males", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.1, 0, 0))
title(ylab = "Allele frequency difference\n(Old - Young) in females", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

text(grconvertX(0.05 , from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)




plot(xqtl.gei.eff[, 1], xqtl.gei.eff[, 2], xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE, pch = 21, bg = "grey", main = "", xlab = "", ylab = "", cex = 0.5)

axis(side = 1, at = seq(-1, 1, 0.5), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

abline(a = 0, b = 1, lty = 2, lwd = 0.5)

title(xlab = "Allele frequency difference\n(Old - Young) at Temp 1", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.1, 0, 0))
title(ylab = "Allele frequency difference\n(Old - Young) at Temp 2", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
