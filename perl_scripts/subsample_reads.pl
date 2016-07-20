#!/usr/bin/env perl;
use warnings;
use strict;

my ( $infile, $index_file, @option ) = @ARGV;

# extract number of lines

my $wc_out = `wc -l $infile`;
my @wc_split = split /\s/, $wc_out;
my $nline = $wc_split[0];

print "Number of lines: $nline\n";

# extract max index

my $index_max = 0;

if ( scalar(@option) == 0 ) {
    print "Max index is NOT specified.\n";
    
    open INDEX, "$index_file" or die "Cannot open $index_file\n";

    while ( my $value = <INDEX> ) {	
	
        chomp $value;
	
        if ( $value > $index_max ) {
            $index_max = $value;
        }
    }
		
    close INDEX;
} else {    
    $index_max = $option[0];
    
    print "Max index is specified as $index_max.\n"
}

print "Maximum index value: $index_max\n";
print "----------------------------\n";

# generate index

for ( my $ind=1 ; $ind <= $index_max ; $ind++ ) {
	my $outfile = $infile."_".$index_file."_subsample_".$ind.".txt";
	
	open IN, "$infile" or die "Cannot open $infile\n";
	open INDEX, "$index_file" or die "Cannot open $index_file\n";
	open OUT, ">$outfile" or die "Cannot open $outfile\n";
	
	my $n_processed = 0;
	
	for ( my $i=1 ; $i <= $nline ; $i++ ) {
	    
		my $line = <IN>;
		my $line_index = <INDEX>;
		
		chomp $line;
		chomp $line_index;
		
	    # write the results
	    
	    if ( $line_index <= $ind ) {
	    	print OUT "$line\n";
	    	$n_processed++;
    	}
	}
	
	close IN;
	close INDEX;    
	close OUT;
	
	print "Index: $ind,\t";
	print "# lines: $n_processed\n";
}
