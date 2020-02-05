args <- commandArgs(TRUE) # args <- c("Myriad Pro", "reportData/all.snp.eff.txt", "mean.eff.corr.pdf")

# file
# ============================================================
plot.con <- data.frame(temp = c("18", "18", "25", "25", "28", "28"), sex = rep(c("Female", "Male"), 3), stringsAsFactors = F)

all.eff <- read.table(args[2], header = F, as.is = T)

#text(3:3, eval(parse(text = paste("expression(paste(", a, ", {}*degree))"))))


file.width = 120
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width/25.4, family = args[1])
par(las = 1, tcl = -0.2, ps = 7, lwd = 0.5, xpd = TRUE)
layout(mat = matrix(1:36, ncol = 6, nrow = 6, byrow = TRUE), widths = c(0.25, rep(0.15, 5)), heights = c(rep(0.15, 5), 0.25))

k = 1

for (i in 1:6) {
	
	mai3 <- 0;
	mai4 <- 0;
	
	for (j in 1:6) {
		
		if (j == 1) { mai2 = file.width/25.4*0.10 } else { mai2 = 0 }
		if (i == 6) { mai1 = file.width/25.4*0.10 } else { mai1 = 0 }
		
		par(mai = c(mai1, mai2, mai3, mai4))		
		
		if (j == i) { 
			
			plot(1:5, 1:5, type = "n", xlab = "", ylab = "", axes = FALSE)
			text(3, 3, eval(parse(text = paste("expression(paste(", plot.con[i, 1], ", {}*degree, \"C \",", plot.con[i, 2], "))"))), cex = 8/par("ps")/par("cex"))
			box()
			
		} else if (j > i) {
			
			
			plot(1:5, 1:5, type = "n", xlab = "", ylab = "", axes = FALSE)
			box()
			
		} else {
			
			smoothScatter(all.eff[, c(j + 1, i + 1)], axes = FALSE, xlab = "", ylab = "", cex = 0.002, pch = 16, xlim = c(-15, 15), ylim = c(-15, 15))
			if (j == 1) { axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(1.5, 0.3, 0), at = seq(-10, 10, 10), lwd = 0.5) }
				if (i == 6) { axis(side = 1, cex.axis = 6/par("ps")/par("cex"), mgp = c(1.5, -0.1, 0), at = seq(-10, 10, 10), lwd = 0.5) }
			text(-8, 12, paste("r = ", formatC(cor(all.eff[, c(i + 1, j + 1)], method = "spearman")[1, 2], format = "f", digit = 2)))
			
		}
				
		}

		
}

dev.off()

