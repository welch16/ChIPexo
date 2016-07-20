#!/usr/bin/env perl;
use warnings;
use strict;

my ( $infile, $outfile, $div_num, $precision, $seed ) = @ARGV;

# extract number of lines

my $wc_out = `wc -l $infile`;
my @wc_split = split /\s/, $wc_out;
my $nline = $wc_split[0];

print "Number of reads: $nline\n";

# calculate random number

my $rand_num = int( $precision * $nline / $div_num );

# generate index

open OUT, ">$outfile" or die "Cannot open $outfile\n";

srand( $seed );

for ( my $i=0 ; $i < $nline ; $i++ ) {
    
	# random number generation
	
	my $read_index = int( rand($rand_num) / $precision ) + 1;
	
    # write the results
    
    print OUT "$read_index\n";
}
    
close OUT;
