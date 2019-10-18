# Compute p value with combined evidence from multiple
# replicates within each group
# ============================================================

args <- commandArgs(TRUE) # args <- c("pileup/pileup.18c.female.filtered.out", "192,192", "192,192", "female", 10, 1000, "pileup/18c.female.bsa.out")
freq.file <- args[1]
nchr1 <- args[2]
nchr2 <- args[3]
sex <- args[4]
tc.min <- as.numeric(args[5])
tc.max <- as.numeric(args[6])
output.file <- args[7]

# read data
# ============================================================
freq.out <- read.table(freq.file, header = FALSE, as.is = TRUE)
nchr1 <- as.numeric(unlist(strsplit(nchr1, split = ",")))
nchr2 <- as.numeric(unlist(strsplit(nchr2, split = ",")))

# loop through each site
# ============================================================

result <- matrix(NA, ncol = 9+length(nchr1)+length(nchr2), nrow = nrow(freq.out))
for (i in 1:5) { result[, i] <- freq.out[, i] }

for (i in 1:nrow(result)) {
  
  alleles1 <- as.numeric(unlist(strsplit(unlist(strsplit(freq.out[i, 6], split = ":")), split = ",")))
  alleles2 <- as.numeric(unlist(strsplit(unlist(strsplit(freq.out[i, 7], split = ":")), split = ",")))
  
  cov1 <- alleles1[seq(1, length(alleles1), 2)] + alleles1[seq(2, length(alleles1), 2)]
  cov2 <- alleles2[seq(1, length(alleles2), 2)] + alleles2[seq(2, length(alleles2), 2)]
  
  if (sum(cov1 >= tc.min & cov1 <= tc.max) == length(cov1) & sum(cov2 >= tc.min & cov2 <= tc.max) == length(cov2)) {
  
    p1 <- alleles1[seq(2, length(alleles1), 2)]/cov1
    p2 <- alleles2[seq(2, length(alleles2), 2)]/cov2
  
    w1 <- cov1/sum(cov1)
    w2 <- cov2/sum(cov2)
  
    p1.ave <- sum(w1*p1)
    p2.ave <- sum(w2*p2)
  
    diff <- p1.ave - p2.ave
  
    if (freq.out[i, 1] == "X" & sex == "male") {
    
      var <- sum(w1^2 * p1*(1-p1) * (1/nchr1*2 + 1/cov1)) + sum(w2^2 * p2*(1-p2) * (1/nchr2*2 + 1/cov2))

    } else {
    
      var <- sum(w1^2 * p1*(1-p1) * (1/nchr1 + 1/cov1)) + sum(w2^2 * p2*(1-p2) * (1/nchr2 + 1/cov2))
  
    }
  
    if (var > 0) {
      test.stat <- diff^2/var
      p <- pchisq(test.stat, df = 1, lower.tail = FALSE)
    } else {
      p <- NA
    } 
  
    result[i, 6:(6+length(nchr1))] <- c(p1, p1.ave)
    result[i, (7+length(nchr1)):(7+length(nchr1) + length(nchr2))] <- c(p2, p2.ave)
    result[i, (8+length(nchr1) + length(nchr2)):(9+length(nchr1)+length(nchr2))] <- c(diff, p)
  
  }  
  
  cat(i, "\n")
  
}

colnames(result) <- c("chr", "pos", "ref", "major", "minor", paste("A", 1:length(nchr1), sep = ""), "A.ave", paste("B", 1:length(nchr2), sep = ""), "B.ave", "A-B", "p")

write.table(result, output.file, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
