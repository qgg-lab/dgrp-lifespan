# ============================================
# = make figure to show overall distribution =
# ============================================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "reportData/gwas.pheno/female.18c.lns.pheno", "reportData/gwas.pheno/female.25c.lns.pheno", "reportData/gwas.pheno/female.28c.lns.pheno", "reportData/gwas.pheno/male.18c.lns.pheno", "reportData/gwas.pheno/male.25c.lns.pheno", "reportData/gwas.pheno/male.28c.lns.pheno")
library("RColorBrewer")

# read data
# ============================================================

female.ln.18c <- read.table(args[2], header = FALSE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))
female.ln.25c <- read.table(args[3], header = FALSE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))
female.ln.28c <- read.table(args[4], header = FALSE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))

male.ln.18c <- read.table(args[5], header = FALSE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))
male.ln.25c <- read.table(args[6], header = FALSE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))
male.ln.28c <- read.table(args[7], header = FALSE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))

rownames(female.ln.18c) <- paste("line_", female.ln.18c[, 1])
rownames(female.ln.25c) <- paste("line_", female.ln.25c[, 1])
rownames(female.ln.28c) <- paste("line_", female.ln.28c[, 1])

rownames(male.ln.18c) <- paste("line_", male.ln.18c[, 1])
rownames(male.ln.25c) <- paste("line_", male.ln.25c[, 1])
rownames(male.ln.28c) <- paste("line_", male.ln.28c[, 1])

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[8], width = file.width/25.4, height = file.width*1.3/3*2/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(2, 2, 1, 0.5), ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 2, 3, 3), ncol = 2, nrow = 2, byrow = TRUE), width = c(0.54, 0.46))

# reaction norms for different temperatures
# ============================================================

unique.lines <- unique(c(rownames(female.ln.18c), rownames(female.ln.25c), rownames(female.ln.28c), rownames(male.ln.18c), rownames(male.ln.25c), rownames(male.ln.28c)))


par(mai = c(file.width/25.4*0.04, file.width/25.4*0.1, file.width/25.4*0.05, file.width/25.4*0.02))
plot(c(0.8, 3.2), c(0, 4), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in unique.lines) {
  
  points(c(1, 2, 3), c(female.ln.18c[i, 2], female.ln.25c[i, 2], female.ln.28c[i, 2]), type = "b", cex = 0.1, pch = 21, col = "#FF000040")
  points(c(1, 2, 3), c(female.ln.18c[i, 2], female.ln.25c[i, 2], female.ln.28c[i, 2]), type = "p", cex = 0.5, pch = 21, bg = brewer.pal(9, "Set1")[c(2, 9, 1)], col = NA)
  
}

axis(side = 1, at = c(1, 2, 3), labels = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), expression(paste("28 ", {}*degree, "C"))), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(ylab = "Micro-environmental", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))
title(ylab = expression(paste("variance (ln", italic(sigma), " lifespan)")), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.6, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

par(family = "Arial")
text(1.5, 0.5, "\u2640", col = "red", cex = 4)
par(family = args[1])


par(mai = c(file.width/25.4*0.04, file.width/25.4*0.02, file.width/25.4*0.05, file.width/25.4*0.02))
plot(c(0.8, 3.2), c(0, 4), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in unique.lines) {
  
  points(c(1, 2, 3), c(male.ln.18c[i, 2], male.ln.25c[i, 2], male.ln.28c[i, 2]), type = "b", cex = 0.1, pch = 21, col = "#0000FF40")
  points(c(1, 2, 3), c(male.ln.18c[i, 2], male.ln.25c[i, 2], male.ln.28c[i, 2]), type = "p", cex = 0.5, pch = 21, bg = brewer.pal(9, "Set1")[c(2, 9, 1)], col = NA)
  
}

axis(side = 1, at = c(1, 2, 3), labels = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), expression(paste("28 ", {}*degree, "C"))), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), at = seq(0, 4, 1), label = rep("", 5), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

par(family = "Arial")
text(1.5, 0.5, "\u2642", col = "blue", cex = 4)
par(family = args[1])



# within each environment
# ============================================================

par(mai = c(file.width/25.4*0.04, file.width/25.4*0.1, file.width/25.4*0.05, file.width/25.4*0.02))
plot(c(0.8, 3.2), c(0, 4), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in unique.lines) {

  points(c(0.8, 1.2), c(female.ln.18c[i, 2], male.ln.18c[i, 2]), type = "b", cex = 0.1, pch = 21, col = brewer.pal(9, "Set1")[2])
  points(c(0.8, 1.2), c(female.ln.18c[i, 2], male.ln.18c[i, 2]), type = "p", cex = 0.5, pch = 21, bg = c("red", "blue"), col = NA)

}


for (i in unique.lines) {

  points(c(1.8, 2.2), c(female.ln.25c[i, 2], male.ln.25c[i, 2]), type = "b", cex = 0.1, pch = 21, col = brewer.pal(9, "Set1")[9])
  points(c(1.8, 2.2), c(female.ln.25c[i, 2], male.ln.25c[i, 2]), type = "p", cex = 0.5, pch = 21, bg = c("red", "blue"), col = NA)

}

for (i in unique.lines) {

  points(c(2.8, 3.2), c(female.ln.28c[i, 2], male.ln.28c[i, 2]), type = "b", cex = 0.1, pch = 21, col = brewer.pal(9, "Set1")[1])
  points(c(2.8, 3.2), c(female.ln.28c[i, 2], male.ln.28c[i, 2]), type = "p", cex = 0.5, pch = 21, bg = c("red", "blue"), col = NA)

}


axis(side = 1, at = c(1, 2, 3), labels = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), expression(paste("28 ", {}*degree, "C"))), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(ylab = "Micro-environmental", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))
title(ylab = expression(paste("variance (ln", italic(sigma), " lifespan)")), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.6, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

par(family = "Arial")
legend("bottomleft", bty = "n", pch = 21, col = c("red", "blue"), pt.cex = 0.6, pt.bg = c("red", "blue"), legend = c("\u2640", "\u2642"), text.col = c("red", "blue"), x.intersp = 0.5, y.intersp = 0.6, cex = 1.5)
par(family = args[1])


dev.off()
