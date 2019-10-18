#! /usr/bin/perl -w

# filter pileup coverage from extractReps.pl

use Getopt::Long;
use List::Util qw(sum);
use Text::NSP::Measures::2D::Fisher::twotailed;

$minTC = 20; # minimum number of bases
$minGB = 20; # minimum number of good bases
$minPercentGB = 0.8; # minimum number of percent good bases
$maxPercentDel = 0.05; # maximum percent indel
$maxTC = 250; # maximum number of coverage
$minPCher = 0.00001; # for p value of the chernoff larger than this number, will be considered sequencing error

$minBA = 0.95; # minimum biallelic proportion
$minAF = 0.05; # minimum allele frequency
$maxFisherP = 0.00001; # maximum fisher p for strand bias

GetOptions("minGB=i" => \$minGB,
	      "minTC=i" => \$minTC,
	      "minPG=f" => \$minPercentGB,
	      "maxPD=f" => \$maxPercentDel,
	      "maxTC=i" => \$maxTC,
	      "minBA=f" => \$minBA,
	      "minAF=f" => \$minAF,
	          "minPC=f" => \$minPCher,
	   "maxFP=f" => \$maxFisherP);

# get lines
while (<>) {
    chomp $_;
    @counts = split /\t/, $_;
    @snp_info = @counts[0..4];
    @sample_counts = split /,/, $counts[5];
    # check sample 1 and sample 2
    print join("\t", @snp_info);
    for (my $i = 0; $i <= $#sample_counts; $i++) {
	@sample_check = countRepFilter(($sample_counts[$i], $minTC, $maxTC, $minGB, $minPercentGB, $maxPercentDel, $minBA, $minAF, $minPCher, $maxFisherP));
	print "\t", $sample_check[0], "\t", $sample_check[1], ",", $sample_check[2]; 
    }
    print "\n";
    
}


sub countRepFilter {
    my @sample_counts = split /:/, shift(@_);
    my ($min_TC, $max_TC, $min_GB, $min_PG, $max_DP, $minBA, $minAF, $minPCher, $maxFisherP) = @_;
    my $tc = $sample_counts[2];
    my $mu = $sample_counts[4];
    my @major_counts = split /\//, $sample_counts[5];
    my @minor_counts = split /\//, $sample_counts[6];
    my $minor_allele_count = sum(@minor_counts);
    my $major_allele_count = sum(@major_counts);

    if ($tc < $min_TC || $tc > $max_TC) {
	return ("TC", $major_allele_count, $minor_allele_count);
    }
    
    my $del_percent = $sample_counts[1]/$tc;

    if ($del_percent >= $max_DP) {
	return ("DP", $major_allele_count, $minor_allele_count);
    }

    my $gb = $sample_counts[3];

    if ($gb < $min_GB) {
	return ("GB", $major_allele_count, $minor_allele_count);
    }

    my $gb_percent = $sample_counts[3]/$sample_counts[2];
    
    if ($gb_percent < $min_PG) {
	return ("PG", $major_allele_count, $minor_allele_count);
    }

    

    
    
    my $ba = ($major_allele_count + $minor_allele_count)/$gb;
    
    if ($ba < $minBA) {
	return("BA", $major_allele_count, $minor_allele_count);
    }

    my $logcherno = log($minPCher/10);

    if ($minor_allele_count > 0) {

	$logcherno = $minor_allele_count - $minor_allele_count*log($minor_allele_count/$mu) - $mu;

    }

    if ($logcherno > log($minPCher)) {
	return("ER", $major_allele_count, $minor_allele_count);
    }

    my $fisherp = calculateStatistic(n11=>$minor_counts[0],n1p=>$minor_allele_count,
				      np1=>$minor_counts[0] + $major_counts[0],
				     npp=>$major_allele_count + $minor_allele_count);
    if ($fisherp < $maxFisherP) {
	return("SB", $major_allele_count, $minor_allele_count);
    }

    my $af = $minor_allele_count/($major_allele_count + $minor_allele_count);
    
    if ($af < $minAF) {
	return("AF", $major_allele_count, $minor_allele_count);
    }

    

    return ("PASS", $major_allele_count, $minor_allele_count);
}
