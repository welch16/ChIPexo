#!/usr/bin/env perl;
use warnings;
use strict;

my ($infile, $outfile, $format, $L) = @ARGV;

# load read file and process it
# (chromosome information is extracted from read file)

open IN, "$infile" or die "Cannot open $infile\n";
open OUT, ">$outfile" or die "Cannot open $outfile\n";

while( my $line = <IN> ){
    chomp $line;
    
    # procee read file, based on "format" option
    
    my @parsed;
    
    if ( $format eq "eland_result" ) {
        @parsed = eland_result( $line, $L );
    } elsif ( $format eq "eland_extended" ) {
        @parsed = eland_extended( $line, $L );
    } elsif ( $format eq "eland_export" ) {
        @parsed = eland_export( $line, $L );
    } elsif ( $format eq "bowtie" ) {
        @parsed = bowtie( $line, $L );
    } elsif ( $format eq "sam" ) {
        @parsed = sam( $line, $L );
    } elsif ( $format eq "bed" ) {
        @parsed = bed( $line, $L );
    } elsif ( $format eq "csem" ) {
        @parsed = csem( $line, $L );
    } else {
        # Unsupported file format -> exit and return 1 to R environment
        
        exit 1;
    }
    
    # skip if invalid line
    
    next if ( $parsed[0] == 0 );
    
    # write the results
    
    print OUT "$line\n";
}
    
close IN;
close OUT;

################################################################################
#                                                                              #
#                              subroutines                                     #
#                                                                              #
################################################################################

# -------------------------- subroutine: eland result -------------------------- 

sub eland_result {
    my $line = shift;   
    my $L = shift;
    
    # parse line    
    
    my ($t1, $seq, $map, $t3, $t4, $t5, $chrt, $pos, $str, @rest) = split /\t/, $line;
    
    # exclude invalid lines
    
    if ( $map eq "QC" || $map eq "NM" ) {
        my @status = 0;
        return @status;
    }
    
    if ( $chrt eq "QC" || $chrt eq "NM" ) {
        my @status = 0;
        return @status;
    }
    
    # exclude multi-reads
    
    if ( $map eq "R0" || $map eq "R1" || $map eq "R2" ) {
        my @status = 0;
        return @status;
    }
    
    if ( $chrt eq "R0" || $chrt eq "R1" || $chrt eq "R2" ) {
        my @status = 0;
        return @status;
    }
    
    my @pos_sets = split( /,/, $pos );
    if ( scalar(@pos_sets) > 1 ) {
        my @status = 0;
        return @status;
    }
    
    # fragment length adjustment
    
    my $read_length = length $seq;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: eland extended ------------------------ 

sub eland_extended {
    my $line = shift;   
    my $L = shift;
    
    # parse line    
    
    my ($t1, $seq, $map, $map_result) = split /\t/, $line;
    
    # exclude invalid lines
    
    if ( $map eq "QC" || $map eq "NM" || $map eq "RM" ) {
        my @status = 0;
        return @status;
    }
    
    # exclude multi-reads
    
    if ( $map ne "1:0:0" && $map ne "0:1:0" && $map ne "0:0:1" ) {
        my @status = 0;
        return @status;
    }
    
    # process chromosome, position, and strand
    
    $map_result =~ m/(chr\w+.fa):(\d+)([FR])/;
        # use "chr\w+.fa" in order to cover chrX & chrY as well
    my $chrt = $1;
    my $pos = $2;
    my $str = $3;
    
    # fragment length adjustment
    
    my $read_length = length $seq;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: eland export -------------------------- 

sub eland_export {
    my $line = shift;   
    my $L = shift;
    
    # parse line    
    
    my ($t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8, $seq, $t9, $chrt, $pos, $str, @rest) = split /\t/, $line;
    
    # exclude invalid lines
    
    if ( $chrt eq "QC" || $chrt eq "NM" ) {
        my @status = 0;
        return @status;
    }
    
    # exclude multi-reads
    
    my @chrt_split = split( /\:/, $chrt );
    if ( scalar(@chrt_split) > 1 ) {
        my @status = 0;
        return @status;
    }
    
    # fragment length adjustment
    
    my $read_length = length $seq;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: bowtie default ------------------------ 

sub bowtie {
    my $line = shift;   
    my $L = shift;
    
    # parse line    
    
    my ( $read_idx, $str_org, $chrt, $pos, $seq, @rest ) = split /\t/, $line;
    
    # exclude invalid lines?
    # exclude multi-reads?
    
    # pos: 0-based offset -> adjust
    
    $pos++;
    
    # str: "+" and "-" -> convert to "F" and "R", respectively
    
    my $str;    
    if ( $str_org eq "+" ) {
        $str = "F";
    } elsif ( $str_org eq "-" ) {
        $str = "R";
    }
    
    # fragment length adjustment
    
    my $read_length = length $seq;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: SAM -----------------------------------

sub sam {
    my $line = shift;   
    my $L = shift;
    
    # parse line
    
    if ( $line =~ /^[^@].+/ ) { 
        # exclude lines starting with "@" (comments)
        
        my ($t1, $bwflag, $chrt, $pos, $t2, $t3, $t4, $t5, $t6, $seq, $t7 ) = split /\t/, $line;
        $pos = int($pos);
        my $str;
        
        if ( $bwflag & 4 or $bwflag & 512 or $bwflag & 1024 ) {
            # exclude invalid lines
                
            my @status = 0;
            return @status;
        } elsif ( $bwflag & 16 ) {
            # reverse strand
            
            $str = "R";
        } else {
            $str = "F";
        }   
    
        # exclude multi-reads?
        
        # fragment length adjustment
        
        my $read_length = length $seq;
        my $L_tmp = $L;  
        #if ( $L < $read_length ) {
        #    $L_tmp = $read_length;
        #}
        
        # return
        
        my $status = 1;
    	my $prob = 1;
        my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
        return @parsed;         
    } else {
        my @status = 0;
        return @status; 
    }
}

# -------------------------- subroutine: BED -------------------------- 

sub bed {
    my $line = shift;   
    my $L = shift;
    
    # skip track line
    
    my ($element1, @element2) = split /\t/, $line;
    
    if ( $element1 eq "track" ) {
        my @status = 0;
        return @status;
    }
    
    # parse line    
    
    my ($chrt, $start, $end, $t1, $t2, $str_org, @rest) = split /\t/, $line;
    
    # pos: 0-based offset -> adjust
    
    my $pos = $start;
    $pos++;
    
    # str: "+" and "-" -> convert to "F" and "R", respectively
    
    my $str;
    if ( $str_org eq "+" ) {
        $str = "F";
    } elsif ( $str_org eq "-" ) {
        $str = "R";
    }
    
    # fragment length adjustment
    
    my $read_length = $end - $start + 1;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: CSEM -------------------------- 

sub csem {
    my $line = shift;   
    my $L = shift;
    
    # skip track line
    
    my ($element1, @element2) = split /\s+/, $line;
    
    if ( $element1 eq "track" ) {
        my @status = 0;
        return @status;
    }
    
    # parse line    
    
    my ($chrt, $start, $end, $t1, $score, $str_org, @rest) = split /\t/, $line;
    
    # pos: 0-based offset -> adjust
    
    my $pos = $start;
    $pos++;
    
    # str: "+" and "-" -> convert to "F" and "R", respectively
    
    my $str;
    if ( $str_org eq "+" ) {
        $str = "F";
    } elsif ( $str_org eq "-" ) {
        $str = "R";
    }
    
    # fragment length adjustment
    
    my $read_length = $end - $start + 1;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = $score / 1000;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}
