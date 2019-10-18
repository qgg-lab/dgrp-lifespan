#! /usr/bin/perl -w

# summarize pileups
# also get average of -

# get options
use Getopt::Long;
use List::Util qw(sum);

#use Statistics::Basic qw(mean median stddev);
#use List::MoreUtils qw(pairwise);
#our ($a, $b);

$minQ = 13; # default minimum base quality
$track = "trackfile.out";
$ref = "qq";
$rmvmono = 1;

GetOptions("minQ=i" => \$minQ,
	      "track=s" => \$track,
	      "mono=i"=> \$rmvmono,
	   "ref=s" => \$ref);

if (!(-e $ref)) {
    die "cannot find ref\n";
}

# read reference
open DMEL, "$ref";
%genome = (); 
while (<DMEL>) {
    chomp($_); 
  if ($_ =~ /^>(.*?) /) {
    $chr = $1;
    $genome{$chr} = "";
    next;
  }
    $genome{$chr} .= uc($_);
}
close DMEL;

if ($rmvmono == 1) {
    open TRACK, ">$track";
}

# main program
while (<>) {
    chomp $_;
    # get sample bases and metrics
    @pileup_line = split /\t/, $_;
    $chr = $pileup_line[0];
    $pos = $pileup_line[1];
    @original_bases = ();
    @original_quals = ();
    @good_bases = ();
    @good_bases_count = ();
    @good_bases_meanq = ();
    @total_length = ();
    @ins_count = ();
    @del_count = ();
    $col = 4;
    while ($col <= $#pileup_line - 1) {
	$current_sample_bases = $pileup_line[$col];
	$current_sample_quals = $pileup_line[$col + 1];

	# 0. clean up get total length of the pileup string;
	$current_sample_bases =~ s/\^.|\$//g; # remove the $ and ^. characters that indicate ends of sequences, this may be in the future a filter.
	# this has a very minor problem for the beginning of a chromosome where one of the samples is represented by * and * for base and qual
	# while the cov should be zero, it is OK to treat it as 1 for computational ease.
	push(@total_length, length($current_sample_quals));
	
	# 1. remove -\d+[ATCGNatcgn]+ in the bases, this does not affect the quality scores
	# deletion will be dealt with at sites where they are deleted.
	@delstring = $current_sample_bases =~ m/-(\d+)[ATCGNatcgn]+/g;
	@uniq_delstring = uniq(@delstring);
	foreach my $delcount (@uniq_delstring) {
	    eval("\$current_sample_bases =~ s/-$delcount\[ATCGNatcgn\]\{$delcount\}//g");
	}

	# 2. remove \+\d+[ATCGNatcgn]+
	@instring = $current_sample_bases =~ m/\+(\d+)[ATCGNatcgn]+/g;
	@uniq_instring = uniq(@instring);
	push(@ins_count, $#instring + 1);
	foreach my $inscount (@uniq_instring) {
	    eval("\$current_sample_bases =~ s/\\+$inscount\[ATCGNatcgn\]\{$inscount\}//g");
	}
	
	# 3. remove "*"s and corresponding quality values
	$current_del_count = 0;
	while ($current_sample_bases =~ m/\*+/) {
	    $current_del_count += $+[0] - $-[0];
	    substr($current_sample_bases, $-[0], $+[0] - $-[0], "");
	    substr($current_sample_quals, $-[0], $+[0] - $-[0], "");
	}
	push(@del_count, $current_del_count);
	push(@original_bases, $current_sample_bases);
	push(@original_quals, $current_sample_quals);
	$col += 3;
    }
    
    # Filter low quality reads
    for (my $sample = 0; $sample <= $#original_bases; $sample++) {
	$current_sample_bases = $original_bases[$sample];
	$current_sample_quals = $original_quals[$sample];
	($current_sample_goodbases, $mean_good_q) = filterQ($current_sample_bases, $current_sample_quals, $minQ);
	push(@good_bases_count, length($current_sample_goodbases));
	push(@good_bases, $current_sample_goodbases);
	push(@good_bases_meanq, $mean_good_q);
    }
    
    # identify major and minor alleles among good bases
    @major_forward_count = ();
    @major_reverse_count = ();
    @minor_forward_count = ();
    @minor_reverse_count = ();
    ($major_allele, $minor_allele) = alleles(join("", @good_bases));
    $ref_allele = substr($genome{$chr}, $pos - 1, 1);

    if ($major_allele eq $minor_allele) {
	if ($rmvmono == 1) {
	    print TRACK $chr, "\t", $pos, "\t", $ref_allele, "\t", $major_allele, "\n";
	    next;
	} else {
	    $minor_allele = "U";
	}
    }

    for (my $sample = 0; $sample <= $#good_bases; $sample++) {
	# get major and minor alleles;
	# count major and minor alleles, also count strands;
	push(@major_forward_count, eval("\$good_bases\[\$sample\] =~ tr/".$major_allele."//"));
	push(@major_reverse_count, eval("\$good_bases\[\$sample\] =~ tr/".lc($major_allele)."//"));
	push(@minor_forward_count, eval("\$good_bases\[\$sample\] =~ tr/".$minor_allele."//"));
	push(@minor_reverse_count, eval("\$good_bases\[\$sample\] =~ tr/".lc($minor_allele)."//"));
    }
    
    # print per sample statistics
    print $chr, "\t", $pos, "\t", $ref_allele, "\t", $major_allele, "\t", $minor_allele, "\t";
    
    @per_sample_output = ();
    
    for (my $sample = 0; $sample <= $#ins_count; $sample++) {
	  push(@per_sample_output, $ins_count[$sample].":".$del_count[$sample].":".$total_length[$sample].":".$good_bases_count[$sample].":".
	              $good_bases_meanq[$sample].":".
	       $major_forward_count[$sample]."/".$major_reverse_count[$sample].":".$minor_forward_count[$sample]."/".$minor_reverse_count[$sample]);
    }

    print join(",", @per_sample_output), "\n";

}

if ($rmvmono == 1) {
    close (TRACK);
}

# subroutines
# unique elements
sub uniq {
    my %seen;
    my @unique = grep { ! $seen{$_}++ } @_;
    return (@unique);
}

# filter bases by quality
sub filterQ {
    my ($base_string, $qual_string, $q) = @_;
    $base_string =~ s/(\$|\^.)//g;
    my @base_split = split //, $base_string;
    my @qual_split = split //, $qual_string;
    my @good_qual = (0);
    for (my $i = 0; $i <= $#qual_split; $i++) {
	my $current_qual = ord($qual_split[$i]) - 33;
	if ($current_qual < $q) {
	    $base_split[$i] = "";
	}
	push(@good_qual, 1/(10**($current_qual/10)));
    }
    my $mean_q = sum(@good_qual);
    return ((join("", @base_split), $mean_q));
}

sub alleles {
    my $allele_str = uc($_[0]);
    my @allele_count = (($allele_str =~ tr/A//), ($allele_str =~ tr/C//),
			($allele_str =~ tr/G//), ($allele_str =~ tr/T//));
    my @allele_code = ("A", "C", "G", "T");
    my $major = -1;
    my $minor = -1;
    my $major_c = 0;
    my $minor_c = 0;

    for (my $i = 0; $i <= 3; $i++) {
	if ($allele_count[$i] == 0) {
	    next;
	}
	    
	if ($allele_count[$i] >= $major_c) {
	    $minor = $major;
	    $minor_c = $major_c;
	    $major = $i;
	    $major_c = $allele_count[$i];
	    next;
	}
	    
	if ($allele_count[$i] < $major_c && $allele_count[$i] >= $minor_c) {
	    $minor = $i;
	    $minor_c = $allele_count[$i];
	}
    }

    if ($major == -1) {
	return(("N", "N"));
    }
        
    if ($minor_c == 0) {
	$minor = $major;
    }

    return (($allele_code[$major], $allele_code[$minor]));
}


