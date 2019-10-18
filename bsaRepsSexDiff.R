# Compute p value with combined evidence from multiple
# replicates within each group
# test for difference of difference
# for example
# one xQTL experiemnt tests for p11 - p12
# the other also tests for p21 - p22
# the combined test tests whether there is a difference between
# these two tests
# wp11-wp12 - wp21 + wp22, where w = weight based on coverage 
# ============================================================

args <- commandArgs(TRUE) # args <- c("pileup/pileup.25c.female.filtered.out", "192,192", "192,192", "pileup/pileup.18c.female.filtered.out", "192,184", "196,196", "female", 10, 1000, "pileup/18c.25c.diff.female.bsa.out")
freq.file1 <- args[1]
nchr11 <- args[2]
nchr12 <- args[3]

freq.file2 <- args[4]
nchr21 <- args[5]
nchr22 <- args[6]

tc.min <- as.numeric(args[7])
tc.max <- as.numeric(args[8])
output.file <- args[9]

# read data
# ============================================================
freq.out1 <- read.table(freq.file1, header = FALSE, as.is = TRUE)
nchr11 <- as.numeric(unlist(strsplit(nchr11, split = ",")))
nchr12 <- as.numeric(unlist(strsplit(nchr12, split = ",")))

freq.out2 <- read.table(freq.file2, header = FALSE, as.is = TRUE)
nchr21 <- as.numeric(unlist(strsplit(nchr21, split = ",")))
nchr22 <- as.numeric(unlist(strsplit(nchr22, split = ",")))

# loop through each site
# ============================================================

result <- matrix(NA, ncol = 7, nrow = nrow(freq.out1))
for (i in 1:5) { result[, i] <- freq.out1[, i] }

for (i in 1:nrow(result)) {
  
  alleles11 <- as.numeric(unlist(strsplit(unlist(strsplit(freq.out1[i, 6], split = ":")), split = ",")))
  alleles12 <- as.numeric(unlist(strsplit(unlist(strsplit(freq.out1[i, 7], split = ":")), split = ",")))
  
  alleles21 <- as.numeric(unlist(strsplit(unlist(strsplit(freq.out2[i, 6], split = ":")), split = ",")))
  alleles22 <- as.numeric(unlist(strsplit(unlist(strsplit(freq.out2[i, 7], split = ":")), split = ",")))
  
  
  cov11 <- alleles11[seq(1, length(alleles11), 2)] + alleles11[seq(2, length(alleles11), 2)]
  cov12 <- alleles12[seq(1, length(alleles12), 2)] + alleles12[seq(2, length(alleles12), 2)]
  
  cov21 <- alleles21[seq(1, length(alleles21), 2)] + alleles21[seq(2, length(alleles21), 2)]
  cov22 <- alleles22[seq(1, length(alleles22), 2)] + alleles22[seq(2, length(alleles22), 2)]
  
  
  if (sum(cov11 >= tc.min & cov11 <= tc.max) == length(cov11) & sum(cov12 >= tc.min & cov12 <= tc.max) == length(cov12) &
      sum(cov21 >= tc.min & cov21 <= tc.max) == length(cov21) & sum(cov22 >= tc.min & cov22 <= tc.max) == length(cov22)) {
  
    p11 <- alleles11[seq(2, length(alleles11), 2)]/cov11
    p12 <- alleles12[seq(2, length(alleles12), 2)]/cov12
    
    p21 <- alleles21[seq(2, length(alleles21), 2)]/cov21
    p22 <- alleles22[seq(2, length(alleles22), 2)]/cov22
      
    w11 <- cov11/sum(cov11)
    w12 <- cov12/sum(cov12)
  
    p11.ave <- sum(w11*p11)
    p12.ave <- sum(w12*p12)
  
    diff1 <- p11.ave - p12.ave
    
    w21 <- cov21/sum(cov21)
    w22 <- cov22/sum(cov22)
  
    p21.ave <- sum(w21*p21)
    p22.ave <- sum(w22*p22)
  
    diff2 <- p21.ave - p22.ave
    
    diff <- diff1 - diff2
    
  
    if (freq.out1[i, 1] == "X") {
    
      var <- sum(w11^2 * p11*(1-p11) * (1/nchr11 + 1/cov11)) + sum(w12^2 * p12*(1-p12) * (1/nchr12 + 1/cov12)) +
             sum(w21^2 * p21*(1-p21) * (1/nchr21*2 + 1/cov21)) + sum(w22^2 * p22*(1-p22) * (1/nchr22*2 + 1/cov22))
    
    } else {
    
      var <- sum(w11^2 * p11*(1-p11) * (1/nchr11 + 1/cov11)) + sum(w12^2 * p12*(1-p12) * (1/nchr12 + 1/cov12)) +
             sum(w21^2 * p21*(1-p21) * (1/nchr21 + 1/cov21)) + sum(w22^2 * p22*(1-p22) * (1/nchr22 + 1/cov22))
      
    }
  
    if (var > 0) {
      test.stat <- diff^2/var
      p <- pchisq(test.stat, df = 1, lower.tail = FALSE)
    } else {
      p <- NA
    } 
  
    result[i, 6:7] <- c(diff, p)
  
  }  
  
  cat(i, "\n")
  
}

colnames(result) <- c("chr", "pos", "ref", "major", "minor", "diff", "p")

write.table(result, output.file, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
