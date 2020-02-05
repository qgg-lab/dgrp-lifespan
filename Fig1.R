# ============================================
# = make figure to show overall distribution =
# ============================================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "reportData/DGRP18CLifespanRawData.csv", "reportData/DGRP25CLifespanRawData.csv", "reportData/DGRP28CLifespanRawData.csv")
library("RColorBrewer")

# read data
# ============================================================

life.18c <- read.csv(args[2], header = TRUE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))
life.25c <- read.csv(args[3], header = TRUE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))
life.28c <- read.csv(args[4], header = TRUE, as.is = TRUE, na.strings = c(".", "", " ", "NA"))

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[5], width = file.width/25.4, height = file.width*1.3/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(2, 2, 1, 0.5), ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 1, 2, 3, 4, 4), ncol = 2, nrow = 3, byrow = TRUE), width = c(0.54, 0.46))

# filter data, require at least 2 females and 2 males
# ============================================================

life.18c <- within(life.18c, vial <- paste("RAL", Line, "-", Rep, sep = ""))
life.25c <- within(life.25c, vial <- paste("RAL", Line, "-", Rep, sep = ""))
life.28c <- within(life.28c, vial <- paste("RAL", Line, "-", Rep, sep = ""))

life.18c <- subset(life.18c, vial %in% intersect(names(which(with(na.omit(life.18c), table(vial[Sex == "Female"])) >= 2)),
                                                 names(which(with(na.omit(life.18c), table(vial[Sex == "Male"])) >= 2))))

life.25c <- subset(life.25c, vial %in% intersect(names(which(with(na.omit(life.25c), table(vial[Sex == "Female"])) >= 2)),
                                                 names(which(with(na.omit(life.25c), table(vial[Sex == "Male"])) >= 2))))

life.28c <- subset(life.28c, vial %in% intersect(names(which(with(na.omit(life.28c), table(vial[Sex == "Female"])) >= 2)),
                                                 names(which(with(na.omit(life.28c), table(vial[Sex == "Male"])) >= 2))))

cat("total number of flies: ", sum(c(nrow(na.omit(life.18c)), nrow(na.omit(life.25c)), nrow(na.omit(life.28c)))), "\n", sep = "")
cat("summary of life span: \n")
print(summary(c(life.18c$Age, life.25c$Age, life.28c$Age)))

# plot data
# ============================================================

life.18c.hist <- hist(life.18c$Age, breaks = seq(0, 200, 10), plot = FALSE)$counts
life.25c.hist <- hist(life.25c$Age, breaks = seq(0, 200, 10), plot = FALSE)$counts
life.28c.hist <- hist(life.28c$Age, breaks = seq(0, 200, 10), plot = FALSE)$counts

par(mai = c(file.width/25.4*0.08, file.width/25.4*0.1, file.width/25.4*0.05, file.width/25.4*0.01))

barplot(rbind(life.18c.hist, life.25c.hist, life.28c.hist), xlim = c(2, 80.5), ylim = c(-100, 8000), beside = TRUE, col = c(brewer.pal(9, "Set1")[c(2, 9, 1)]), border = NA, axes = FALSE, xlab = "", ylab = "")
axis(side = 1, at = c(0, 20, 40, 60, 80) + 0.5, labels = c(0, 50, 100, 150, 200), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(xlab = "Lifespan (days)", cex.lab = 7/par("ps")/par("cex"), mgp = c(0.7, 0, 0))
title(ylab = "Number of flies", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.6, 0, 0))

# add text
# ============================================================

text(16, 7000, expression(paste("28 ", {}*degree, "C")), pos = 4, col = brewer.pal(9, "Set1")[1], cex = 7/par("ps")/par("cex"))
text(20, 5000, expression(paste("25 ", {}*degree, "C")), pos = 4, col = brewer.pal(9, "Set1")[9], cex = 7/par("ps")/par("cex"))
text(40, 2500, expression(paste("18 ", {}*degree, "C")), pos = 4, col = brewer.pal(9, "Set1")[2], cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# reaction norms for different temperatures
# ============================================================

life.18c.female <- with(subset(life.18c, Sex == "Female"), sapply(split(Age, paste("line_", Line, sep = "")), mean, na.rm = TRUE))
life.18c.male <- with(subset(life.18c, Sex == "Male"), sapply(split(Age, paste("line_", Line, sep = "")), mean, na.rm = TRUE))

life.25c.female <- with(subset(life.25c, Sex == "Female"), sapply(split(Age, paste("line_", Line, sep = "")), mean, na.rm = TRUE))
life.25c.male <- with(subset(life.25c, Sex == "Male"), sapply(split(Age, paste("line_", Line, sep = "")), mean, na.rm = TRUE))

life.28c.female <- with(subset(life.28c, Sex == "Female"), sapply(split(Age, paste("line_", Line, sep = "")), mean, na.rm = TRUE))
life.28c.male <- with(subset(life.28c, Sex == "Male"), sapply(split(Age, paste("line_", Line, sep = "")), mean, na.rm = TRUE))

unique.lines <- unique(paste("line_", c(life.18c$Line, life.25c$Line, life.28c$Line), sep = ""))


par(mai = c(file.width/25.4*0.04, file.width/25.4*0.1, file.width/25.4*0.05, file.width/25.4*0.02))
plot(c(0.8, 3.2), c(0, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in unique.lines) {
  
  points(c(1, 2, 3), c(life.18c.female[i], life.25c.female[i], life.28c.female[i]), type = "b", cex = 0.1, pch = 21, col = "#FF000040")
  points(c(1, 2, 3), c(life.18c.female[i], life.25c.female[i], life.28c.female[i]), type = "p", cex = 0.5, pch = 21, bg = brewer.pal(9, "Set1")[c(2, 9, 1)], col = NA)
  
}

axis(side = 1, at = c(1, 2, 3), labels = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), expression(paste("28 ", {}*degree, "C"))), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(ylab = "Lifespan (days)", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.6, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

par(family = "Arial")
text(3, 120, "\u2640", col = "red", cex = 5)
par(family = args[1])


par(mai = c(file.width/25.4*0.04, file.width/25.4*0.02, file.width/25.4*0.05, file.width/25.4*0.02))
plot(c(0.8, 3.2), c(0, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in unique.lines) {
  
  points(c(1, 2, 3), c(life.18c.male[i], life.25c.male[i], life.28c.male[i]), type = "b", cex = 0.1, pch = 21, col = "#0000FF40")
  points(c(1, 2, 3), c(life.18c.male[i], life.25c.male[i], life.28c.male[i]), type = "p", cex = 0.5, pch = 21, bg = brewer.pal(9, "Set1")[c(2, 9, 1)], col = NA)
  
}

axis(side = 1, at = c(1, 2, 3), labels = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), expression(paste("28 ", {}*degree, "C"))), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), at = seq(0, 150, 50), label = rep("", 4), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(ylab = "Lifespan (days)", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.6, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

par(family = "Arial")
text(3, 120, "\u2642", col = "blue", cex = 5)
par(family = args[1])

# within each environment
# ============================================================

par(mai = c(file.width/25.4*0.04, file.width/25.4*0.1, file.width/25.4*0.05, file.width/25.4*0.02))
plot(c(0.8, 3.2), c(0, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in unique.lines) {
  
  points(c(0.8, 1.2), c(life.18c.female[i], life.18c.male[i]), type = "b", cex = 0.1, pch = 21, col = brewer.pal(9, "Set1")[2])
  points(c(0.8, 1.2), c(life.18c.female[i], life.18c.male[i]), type = "p", cex = 0.5, pch = 21, bg = c("red", "blue"), col = NA)
  
}


for (i in unique.lines) {
  
  points(c(1.8, 2.2), c(life.25c.female[i], life.25c.male[i]), type = "b", cex = 0.1, pch = 21, col = brewer.pal(9, "Set1")[9])
  points(c(1.8, 2.2), c(life.25c.female[i], life.25c.male[i]), type = "p", cex = 0.5, pch = 21, bg = c("red", "blue"), col = NA)
  
}

for (i in unique.lines) {
  
  points(c(2.8, 3.2), c(life.28c.female[i], life.28c.male[i]), type = "b", cex = 0.1, pch = 21, col = brewer.pal(9, "Set1")[1])
  points(c(2.8, 3.2), c(life.28c.female[i], life.28c.male[i]), type = "p", cex = 0.5, pch = 21, bg = c("red", "blue"), col = NA)
  
}


axis(side = 1, at = c(1, 2, 3), labels = c(expression(paste("18 ", {}*degree, "C")), expression(paste("25 ", {}*degree, "C")), expression(paste("28 ", {}*degree, "C"))), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

title(ylab = "Lifespan (days)", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.6, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 8/par("ps")/par("cex"), pos = 1, xpd = TRUE)

par(family = "Arial")
legend("topright", bty = "n", pch = 21, col = c("red", "blue"), pt.cex = 0.6, pt.bg = c("red", "blue"), legend = c("\u2640", "\u2642"), text.col = c("red", "blue"), x.intersp = 0.5, y.intersp = 0.6, cex = 1.5)
par(family = args[1])


dev.off()
