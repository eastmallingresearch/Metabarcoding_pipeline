#! /usr/bin/perl 
#use warnings;
use strict;
use Scalar::Util qw(looks_like_number);
use List::Util qw(sum);;

#################################################################################
#
# Demultiplex paired fastq with a set of primers
# Unidentified sequences are written ambig output files
#
# Usage: demulti.pl forward_read reverse_read mismatch i1 i2 i..n i..n+1
#
################################################################################

pop(@ARGV);

my $r1=$ARGV[0];
my $r2=$ARGV[1];

open (my $R1,'<', shift) or die "Could not open R1: $ARGV[0]";
open (my $R2,'<', shift) or die "Could not open R2: $ARGV[1]";

my $allowmiss;
if (looks_like_number($ARGV[0])) {$allowmiss = shift;}	
else {$allowmiss = 0;}

my @ind_f;
my @ind_r;
while(scalar(@ARGV)>0) {
	push @ind_f,shift;
	push @ind_r,shift;		
}
chomp(@ind_f);
chomp(@ind_r);

my @lf = map length, @ind_f;
my @lr = map length, @ind_r;

my %fileNames;
my $x=0;
for (my $i=1;$i<=scalar(@ind_f);$i++) {
	$fileNames{$ind_f[$i-1]}="$r1.ps$i.fastq";
	$fileNames{$ind_r[$i-1]}="$r2.ps$i.fastq";
}
$fileNames{"ambig_f"}="$r1.ambig.fastq";
$fileNames{"ambig_r"}="$r2.ambig.fastq";


my %fileHandles;
for my $filename (keys %fileNames) {
    open $fileHandles{$filename}, '>', $fileNames{$filename} or die "Can't open $filename: $!";
}

my $R1_l = "";
my $R2_l = "";
my $counter=1;
my %output;
my $destringyfied;

do {
	my $l1 = <$R1>;
	my $l2 = <$R2>;

	$R1_l.=$l1;
	$R2_l.=$l2;
	
	if(($counter+2)%4==0) {
	
		my $countf;
		my $countr;
		
		for(my $i=0;$i<scalar(@ind_f);$i++){
			$countf=count_match(substr($l1,0,$lf[$i]),$ind_f[$i]);
			$countr=count_match(substr($l2,0,$lr[$i]),$ind_r[$i]);
			 
			if($countf <= $allowmiss && $countr <= $allowmiss){
				$output{$i}=1;
			} else {
				$output{$i}=0;
			}
		}
	}
	
	#my $sum = sum(%output);

	if($counter%4==0) {
		my $x=sum values %output;
		if ($x==1) {
			for my $key (keys %output) {
				if($output{$key}==1){
					$destringyfied = $fileHandles{$ind_f[$key]};
					print $destringyfied "$R1_l";
					$destringyfied = $fileHandles{$ind_r[$key]};
					print $destringyfied "$R2_l";				
				}
				
			}
		} else {
			$destringyfied = $fileHandles{"ambig_f"};
			print $destringyfied "$R1_l";
			$destringyfied = $fileHandles{"ambig_r"};
			print $destringyfied "$R2_l";
		}
		$destringyfied="";
		$R1_l="";
		$R2_l="";
		%output=();
	}

	$counter++;
} until eof $R1;


sub count_match {
	my $str1 = uc(shift);
	my $str2 = uc(shift);
	my $count = length($str1);
	my $l = length($str1);
	for (my $i=0;$i<$l;$i++) {
		my $x = to_bin(substr($str1,$i,1));
		my $y = to_bin(substr($str2,$i,1));
		$count -- if $x&$y;
 	}
	return($count);
	
}

sub to_bin{
	my $t = @_[0];
	my $score = 0;
	$score ++ if $t=~/[ARWMDHVN]/;
	$score =$score+2 if $t=~/[CYSMBHVN]/;
	$score =$score+4 if $t=~/[GRSKBDVN]/;
	$score =$score+8 if $t=~/[TYWKBDHN]/;
	return($score);	
}

