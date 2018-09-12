#! /usr/bin/perl 
#use warnings;
use strict;
use List::UtilsBy qw(max_by);

my $counter=1;
my %SEQS;

pop(@ARGV);

if (looks_like_number($ARGV[0])) {

} else {

}

my $FILE=$ARGV[0];
my $LENGTH=$ARGV[1];

open (my $FILE,'<', shift) or die "Could not open FILE: $ARGV[0]";

do {
	
	my $seq=<$FILE>;	

	if(($counter+2)%4==0) {
			$SEQS{substr($seq,0,($LENGTH-1))}++
	}

} until eof $FILE;

my $adapter = max_by { $SEQS{$_} } keys %SEQS;

print "$adapter";