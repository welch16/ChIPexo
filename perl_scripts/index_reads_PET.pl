#!/usr/bin/env perl;
use warnings;
use strict;

my ( $infile, $outfile, $div_num, $precision, $seed ) = @ARGV;

# extract number of lines

my $wc_out = `wc -l $infile`;
my @wc_split = split /\s/, $wc_out;
my $nline = $wc_split[0];
my $nfrag = $nline / 2;

print "Number of reads: $nline\n";
print "Number of fragments: $nfrag\n";

# calculate random number (PET)

my $rand_num = int( ( $precision * $nline ) / ( 2 * $div_num ) );

# generate index (PET)

open OUT, ">$outfile" or die "Cannot open $outfile\n";

srand( $seed );

for ( my $i=0 ; $i < $nfrag ; $i++ ) {
    
	# random number generation
	
	my $read_index = int( rand($rand_num) / $precision ) + 1;
	
    # write the results
    
    print OUT "$read_index\n";
    print OUT "$read_index\n";
}
    
close OUT;
