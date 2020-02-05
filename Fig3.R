# ============================================
# = make figure to show overall distribution =
# ============================================

args <- commandArgs(TRUE) # args <- c("reportData/life.gwas.txt", "reportData/lifespan.eff.ave.txt", "reportData/lifespan.lns.eff.ave.txt", "report/figure_eff.size.pdf")
library("RColorBrewer")

# read data
# ============================================================

life.gwas.res <- read.table(args[2], header = FALSE, as.is = TRUE)
life.gwas.res <- life.gwas.res[life.gwas.res[, 1] %in% c("18c.mean", "25c.mean", "28c.mean"), ]

eff.size.sig <- c(life.gwas.res[life.gwas.res[, 10] < 1e-5, 9],
                  life.gwas.res[life.gwas.res[, 12] < 1e-5, 11])

# prepare file to plot
# ============================================================

file.width = 89*0.6
cairo_pdf(file = args[5], width = file.width/25.4, height = file.width/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(2, 3, 1, 0.5), ps = 7, lwd = 0.5)

# plot data
# ============================================================

hist(eff.size.sig, breaks = seq(-20, 20, 1), axes = FALSE, col = "grey", main = "", xlab = "", ylab = "")

axis(side = 1, at = seq(-20, 20, 5), cex.axis = 6/par("ps")/par("cex"), mgp = c(2.1, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(xlab = "Allelic effect for lifespan", cex.lab = 7/par("ps")/par("cex"), mgp = c(0.7, 0, 0))
title(ylab = "Number of significant variants", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))


dev.off()
