args <- commandArgs(TRUE)

geno <- read.table(args[1], as.is = TRUE, header = FALSE, na.strings = "0")
rownames(geno) <- geno[, 1]
map <- read.table(args[2], as.is = TRUE, header = FALSE)

heat.col <- colorRampPalette(c("white", "red"))(20)

png(file = args[3], width = 8, height = 8, unit = "in", res = 300)

# Reorder things such that chromosomes are arranged as X, 2L, 3L, 3L, 3R, 4
map[map[, 1] == 5, 1] <- 0
geno <- geno[, order(map[, 1], map[, 4])]
map <- map[order(map[, 1], map[, 4]), ]

# Find gene boundaries
chrX.start <- suppressWarnings(min(which(map[, 1] == 0)))
chrX.end <- suppressWarnings(max(which(map[, 1] == 0)))

chr2L.start <- suppressWarnings(min(which(map[, 1] == 1)))
chr2L.end <- suppressWarnings(max(which(map[, 1] == 1)))

chr2R.start <- suppressWarnings(min(which(map[, 1] == 2)))
chr2R.end <- suppressWarnings(max(which(map[, 1] == 2)))

chr3L.start <- suppressWarnings(min(which(map[, 1] == 3)))
chr3L.end <- suppressWarnings(max(which(map[, 1] == 3)))

chr3R.start <- suppressWarnings(min(which(map[, 1] == 4)))
chr3R.end <- suppressWarnings(max(which(map[, 1] == 4)))

chr4.start <- suppressWarnings(min(which(map[, 1] == 6)))
chr4.end <- suppressWarnings(max(which(map[, 1] == 6)))

geno.cor <- cor(matrix(as.numeric(unlist(geno)), nrow = nrow(geno), byrow = T), use = "pairwise")^2
geno.cor[is.na(geno.cor)] <- 0

plot(c(0, nrow(geno.cor)*1.2), c(0, nrow(geno.cor)*1.2), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
image(1:nrow(geno.cor), 1:ncol(geno.cor), geno.cor, col = heat.col, breaks = seq(0, 1, 0.05), add = TRUE)
rect(0.5, 0.5, nrow(geno.cor) + 0.5, nrow(geno.cor) + 0.5, lwd = 1)

# Plot boundaries
if (length(chrX.start) > 0) {
  segments(chrX.end + 0.5, 0.5, chrX.end + 0.5, nrow(geno.cor) + 0.5, lwd = 1)
  segments(0.5, chrX.end + 0.5, nrow(geno.cor) + 0.5, chrX.end + 0.5, lwd = 1)
  text((chrX.start + chrX.end)/2, nrow(geno.cor), labels = "X", pos = 3)
}

if (length(chr2L.start) > 0) {
  segments(chr2L.end + 0.5, 0.5, chr2L.end + 0.5, nrow(geno.cor) + 0.5, lwd = 1)
  segments(0.5, chr2L.end + 0.5, nrow(geno.cor) + 0.5, chr2L.end + 0.5, lwd = 1)
  text((chr2L.start + chr2L.end)/2, nrow(geno.cor), labels = "2L", pos = 3)
}

if (length(chr2R.start) > 0) {
  segments(chr2R.end + 0.5, 0.5, chr2R.end + 0.5, nrow(geno.cor) + 0.5, lwd = 1)
  segments(0.5, chr2R.end + 0.5, nrow(geno.cor) + 0.5, chr2R.end + 0.5, lwd = 1)
  text((chr2R.start + chr2R.end)/2, nrow(geno.cor), labels = "2R", pos = 3)
}

if (length(chr3L.start) > 0) {
  segments(chr3L.end + 0.5, 0.5, chr3L.end + 0.5, nrow(geno.cor) + 0.5, lwd = 1)
  segments(0.5, chr3L.end + 0.5, nrow(geno.cor) + 0.5, chr3L.end + 0.5, lwd = 1)
  text((chr3L.start + chr3L.end)/2, nrow(geno.cor), labels = "3L", pos = 3)
}

if (length(chr3R.start) > 0) {
  segments(chr3R.end + 0.5, 0.5, chr3R.end + 0.5, nrow(geno.cor) + 0.5, lwd = 1)
  segments(0.5, chr3R.end + 0.5, nrow(geno.cor) + 0.5, chr3R.end + 0.5, lwd = 1)
  text((chr3R.start + chr3R.end)/2, nrow(geno.cor), labels = "3R", pos = 3)
}

if (length(chr4.start) > 0) {
  text((chr4.start + chr4.end)/2, nrow(geno.cor), labels = "4", pos = 3)
}

# plot chromosome and connecting lines
# ============================================================

chr.len <- c(22422827, 23011544, 21146708, 24543557, 27905053, 1351857)
names(chr.len) <- c("X", "2L", "2R", "3L", "3R", "4")
chr.gap <- c(0, 2000000, 2000000, 4000000, 4000000, 6000000)
chr.cum <- cumsum(c(0, chr.len[1:5])) + chr.gap
names(chr.cum) <- names(chr.len)

# the coordinates need to be converted to the axis of the heatmap
# ============================================================

n.max <- ncol(geno.cor)

for (chr in names(chr.len)) {
  axis(side = 1, at = c(chr.cum[chr], chr.len[chr] + chr.cum[chr])/(chr.cum["4"] + chr.len["4"])*n.max, mgp = c(2, 2, 1), label = c("", ""))
  text((chr.cum[chr] + chr.len[chr]/2)/(chr.cum["4"] + chr.len["4"])*n.max, -150, chr, xpd = T, cex = 0.8)
}

for (i in 1:nrow(map)) {
  
  info = unlist(strsplit(map[i, 2], split = "_"))
  segments(i, 0, (chr.cum[info[1]] + as.numeric(info[2]))/(chr.cum["4"] + chr.len["4"])*n.max, -100, xpd = T, lwd = 0.2)
  
}


# Plot scale bar
r20 <- seq(0.01, 1, by = 0.01)
n.col <- length(r20)
seg.width <- nrow(geno.cor)/n.col
seg.pos <- 0.5 + 0.5*seg.width + seg.width*(0:(n.col - 1))

image(c(nrow(geno.cor)*1.02, nrow(geno.cor)*1.04), seg.pos, matrix(rep(r20, 2), nrow = 2, byrow = T), col = heat.col, add = T, useRaster = T)
rect(nrow(geno.cor)*1.01, 0.5, nrow(geno.cor)*1.05, nrow(geno.cor) + 0.5)

axis(side = 2, pos = nrow(geno.cor)*1.01, at = seq(0, 1, 0.2)*nrow(geno.cor) + 0.5, tcl = 0.2, labels = rep("", 6))

axis(side = 4, pos = nrow(geno.cor)*1.05, at = seq(0, 1, 0.2)*nrow(geno.cor) + 0.5, tcl = 0.2, labels = seq(0, 1, 0.2), mgp = c(2, 0.2, 0), las = 2, cex.axis = 0.8)

# plot snp position


dev.off()
