# ============================================
# = make figure to show overall distribution =
# ============================================

args <- commandArgs(TRUE) # args <- c("reportData/Formatted_18Degree_RNAi.csv", "reportData/Formatted_25Degree_RNAi.csv", "reportData/Formatted_28Degree_RNAi.csv")
library("RColorBrewer")

# read data
# ============================================================
gene.name <- c("bru-3", "CG11828", "CG13921", "CG42750", "CG9265", "Cht-2", "Cka", "Doa", "E23", "Eip75B", "olf413", "PGRP-LA", "PGRP-LC", "Rdl", "tou")

rnai.18c <- read.csv(args[2], header = TRUE, as.is = TRUE)
rnai.18c[rnai.18c[, 2] == "cka", 2] <- "Cka"
rnai.18c[rnai.18c[, 2] == "Tou", 2] <- "tou"

rnai.25c <- read.csv(args[3], header = TRUE, as.is = TRUE)
rnai.25c[rnai.25c[, 2] == "Cht2", 2] <- "Cht-2"
rnai.25c[rnai.25c[, 2] == "KK control", 2] <- "Control"

rnai.28c <- read.csv(args[4], header = TRUE, as.is = TRUE)
rnai.28c[rnai.28c[, 2] == "cka", 2] <- "Cka"
rnai.28c[rnai.28c[, 2] == "Tou", 2] <- "tou"


# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[5], width = file.width/25.4, height = file.width/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(2.4, 2.5, 1.5, 0.5), ps = 7, lwd = 0.5, mfrow = c(3, 1), xpd = TRUE)

gene.order <- c("Control", setdiff(names(sort(with(rnai.18c[rnai.18c[, 4] == "Female",], sapply(split(Age, Gene), mean)))), "Control"))



plot(c(1, 16), c(-30, 30), type = "n", axes = FALSE, xlab = "", ylab = "")
base.line <- sort(with(rnai.18c[rnai.18c[, 4] == "Female",], sapply(split(Age, Gene), mean)))["Control"]


for (i in 1:length(gene.order)) {
	
	if (i %% 2 == 0) {
		
		rect(i - 0.5, par('usr')[3], i + 0.5, par('usr')[4], col = "grey90", border = NA)
		
	}
	
	this.gene.female <- rnai.18c[rnai.18c[, 4] == "Female" & rnai.18c[, 2] == gene.order[i], 5] - base.line
	#points(i - 0.2 + runif(length(this.gene.female), -0.1, 0.1), this.gene.female, pch = 1, col = brewer.pal(9, "Reds")[6], cex = 0.1, lwd = 0.2)
	
	segments(i - 0.2 - 0.1, mean(this.gene.female), i - 0.2 + 0.1, mean(this.gene.female), col = brewer.pal(9, "Reds")[9])
	segments(i - 0.2, mean(this.gene.female) - sd(this.gene.female)/sqrt(length(this.gene.female)), i - 0.2, mean(this.gene.female) + sd(this.gene.female)/sqrt(length(this.gene.female)), col = brewer.pal(9, "Reds")[9])
	
	this.gene.male <- rnai.18c[rnai.18c[, 4] == "Male" & rnai.18c[, 2] == gene.order[i], 5] - base.line
	#points(i + 0.2 + runif(length(this.gene.male), -0.1, 0.1), this.gene.male, pch = 1, col = brewer.pal(9, "Blues")[6], cex = 0.1, lwd = 0.2)
	
	segments(i + 0.2 - 0.1, mean(this.gene.male), i + 0.2 + 0.1, mean(this.gene.male), col = brewer.pal(9, "Blues")[9])
	segments(i + 0.2, mean(this.gene.male) - sd(this.gene.male)/sqrt(length(this.gene.male)), i + 0.2, mean(this.gene.male) + sd(this.gene.male)/sqrt(length(this.gene.male)), col = brewer.pal(9, "Blues")[9])
	
	
}

legend("top", lwd = 1, col = c(brewer.pal(9, "Reds")[9], brewer.pal(9, "Blues")[9]), legend = c("Females", "Males"), bty = "n", x.intersp = 0.5, y.intersp = 0.5)

axis(side = 1, at = 1:16, label = NA, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
text(1:16 + 0.5, -38, pos = 2, gene.order, srt = 45, cex = 5/par("ps")/par("cex"))
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")
title(ylab = "Deviation from female controls", mgp = c(1.4, 0, 0), cex = 7/par("ps")/par("cex"))
text(0.2, 30, expression(paste("18 ", {}*degree, "C")), pos = 4, cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05 , from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)




plot(c(1, 16), c(-15, 15), type = "n", axes = FALSE, xlab = "", ylab = "")
base.line <- sort(with(rnai.25c[rnai.25c[, 4] == "Female",], sapply(split(Age, Gene), mean)))["Control"]

for (i in 1:length(gene.order)) {
	
	if (i %% 2 == 0) {
		
		rect(i - 0.5, par('usr')[3], i + 0.5, par('usr')[4], col = "grey90", border = NA)
		
	}
	
	this.gene.female <- rnai.25c[rnai.25c[, 4] == "Female" & rnai.25c[, 2] == gene.order[i], 5] - base.line
	#points(i - 0.2 + runif(length(this.gene.female), -0.1, 0.1), this.gene.female, pch = 1, col = brewer.pal(9, "Reds")[6], cex = 0.1, lwd = 0.2)
	
	segments(i - 0.2 - 0.1, mean(this.gene.female), i - 0.2 + 0.1, mean(this.gene.female), col = brewer.pal(9, "Reds")[9])
	segments(i - 0.2, mean(this.gene.female) - sd(this.gene.female)/sqrt(length(this.gene.female)), i - 0.2, mean(this.gene.female) + sd(this.gene.female)/sqrt(length(this.gene.female)), col = brewer.pal(9, "Reds")[9])
	
	this.gene.male <- rnai.25c[rnai.25c[, 4] == "Male" & rnai.25c[, 2] == gene.order[i], 5] - base.line
	#points(i + 0.2 + runif(length(this.gene.male), -0.1, 0.1), this.gene.male, pch = 1, col = brewer.pal(9, "Blues")[6], cex = 0.1, lwd = 0.2)
	
	segments(i + 0.2 - 0.1, mean(this.gene.male), i + 0.2 + 0.1, mean(this.gene.male), col = brewer.pal(9, "Blues")[9])
	segments(i + 0.2, mean(this.gene.male) - sd(this.gene.male)/sqrt(length(this.gene.male)), i + 0.2, mean(this.gene.male) + sd(this.gene.male)/sqrt(length(this.gene.male)), col = brewer.pal(9, "Blues")[9])
	
	
}

axis(side = 1, at = 1:16, label = NA, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
text(1:16 + 0.5, -19, pos = 2, gene.order, srt = 45, cex = 5/par("ps")/par("cex"))
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")
title(ylab = "Deviation from female controls", mgp = c(1.4, 0, 0), cex = 7/par("ps")/par("cex"))
text(0.2, 15, expression(paste("25 ", {}*degree, "C")), pos = 4, cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05 , from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)






plot(c(1, 16), c(-15, 15), type = "n", axes = FALSE, xlab = "", ylab = "")
base.line <- sort(with(rnai.28c[rnai.28c[, 4] == "Female",], sapply(split(Age, Gene), mean)))["Control"]

for (i in 1:length(gene.order)) {
	if (i %% 2 == 0) {
		
		rect(i - 0.5, par('usr')[3], i + 0.5, par('usr')[4], col = "grey90", border = NA)
		
	}
	
	this.gene.female <- rnai.28c[rnai.28c[, 4] == "Female" & rnai.28c[, 2] == gene.order[i], 5] - base.line
	#points(i - 0.2 + runif(length(this.gene.female), -0.1, 0.1), this.gene.female, pch = 1, col = brewer.pal(9, "Reds")[6], cex = 0.1, lwd = 0.2)
	
	segments(i - 0.2 - 0.1, mean(this.gene.female), i - 0.2 + 0.1, mean(this.gene.female), col = brewer.pal(9, "Reds")[9])
	segments(i - 0.2, mean(this.gene.female) - sd(this.gene.female)/sqrt(length(this.gene.female)), i - 0.2, mean(this.gene.female) + sd(this.gene.female)/sqrt(length(this.gene.female)), col = brewer.pal(9, "Reds")[9])
	
	this.gene.male <- rnai.28c[rnai.28c[, 4] == "Male" & rnai.28c[, 2] == gene.order[i], 5] - base.line
	#points(i + 0.2 + runif(length(this.gene.male), -0.1, 0.1), this.gene.male, pch = 1, col = brewer.pal(9, "Blues")[6], cex = 0.1, lwd = 0.2)
	
	segments(i + 0.2 - 0.1, mean(this.gene.male), i + 0.2 + 0.1, mean(this.gene.male), col = brewer.pal(9, "Blues")[9])
	segments(i + 0.2, mean(this.gene.male) - sd(this.gene.male)/sqrt(length(this.gene.male)), i + 0.2, mean(this.gene.male) + sd(this.gene.male)/sqrt(length(this.gene.male)), col = brewer.pal(9, "Blues")[9])
	
	
}

axis(side = 1, at = 1:16, label = NA, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
text(1:16 + 0.5, -19, pos = 2, gene.order, srt = 45, cex = 5/par("ps")/par("cex"))
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")
title(ylab = "Deviation from female controls", mgp = c(1.4, 0, 0), cex = 7/par("ps")/par("cex"))
text(0.2, 15, expression(paste("28 ", {}*degree, "C")), pos = 4, cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05 , from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
