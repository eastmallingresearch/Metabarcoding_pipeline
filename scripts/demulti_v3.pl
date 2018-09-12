#! /usr/bin/perl 
#use warnings;
use strict;
use Scalar::Util qw(looks_like_number);
use List::Util qw(sum);
use List::UtilsBy qw(max_by);


#############################################################################################
#
# Demultiplex paired fastq with a set of (IPAC) primers
# Unidentified sequences are written ambig output files
#
# Usage: demulti_v3.pl forward_read reverse_read mismatch* i1 i2 i..n i..n+1
#
# * Mismatch has max of 9. If negative forward reads only. Final digit is no. allowed mismatches
#  
#
###################################################################################################

# globals
my $allowmiss;
my @ind_f;
my @ind_r;
my %fileNames;
my %fileHandles;
my $forward_only=0;
my $adapter_length=0;
my $R1_l = "";
my $R2_l = "";
my $counter=1;
my %output;
my $destringyfied;
my @f4000;
my @r4000;

# constants
use constant N => 4000;

# final command line argument when run from PIPELINE.sh must be removed
pop(@ARGV);

# First two arguments are the sequence files
my $r1=$ARGV[0];
my $r2=$ARGV[1];

# open the seqeunce filesa and assign to filehandlers
open (my $R1,'<', shift) or die "Could not open R1: $ARGV[0]";
open (my $R2,'<', shift) or die "Could not open R2: $ARGV[1]";


# check if arg3 (mismatches) is set 
if (looks_like_number($ARGV[0])) {
	$allowmiss = shift;
	if($allowmiss<0) {
		$allowmiss = $allowmiss * -1;
		$forward_only=1;
	}
	if($allowmiss>9) {
		$adapter_length = ($allowmiss - ( $allowmiss % 10 ))/10;
		$allowmiss= $allowmiss % 10;
	}
}	
else {$allowmiss = 0;}

# add the forward and reverse primers to arrays
while(scalar(@ARGV)>0) {
	push @ind_f,shift;
	push @ind_r,shift;		
}
chomp(@ind_f);
chomp(@ind_r);

while (<$R1>) {
	push @f4000, $_;
	last if @f4000 == N || eof;
}
while (<$R2>) {
	push @r4000, $_;
	last if @r4000 == N || eof;
}

# if the mismatch argument is also set to look for adapters - find them and append to the primer arrays
grab_adapters() if $adapter_length; # (\@ind_f, \@ind_r ) no point passing by ref, set globaly

#my $x=0;

for (my $i=1;$i<=scalar(@ind_f);$i++) {
	$fileNames{$ind_f[$i-1]}="$r1.ps$i.fastq";
	$fileNames{$ind_r[$i-1]}="$r2.ps$i.fastq";
}
$fileNames{"ambig_f"}="$r1.ambig.fastq";
$fileNames{"ambig_r"}="$r2.ambig.fastq";

my @lf = map length, @ind_f;
my @lr = map length, @ind_r;

for my $filename (keys %fileNames) {
    open $fileHandles{$filename}, '>', $fileNames{$filename} or die "Can't open $filename: $!";
}

for($counter;$counter<=scalar(@r4000);$counter++){
	do_things($f4000[($counter-1)],$r4000[($counter-1)])
}

do {
	my $in1 = <$R1>;
	my $in2 = <$R2>;
	do_things($in1,$in2);
	$counter++;
} until eof $R1;

sub do_things {
	my $l1 = shift;
	my $l2 = shift;
	
	$R1_l.=$l1;
	$R2_l.=$l2;
	
	if(($counter+2)%4==0) {
	
		my $countf;
		my $countr;
		
		for(my $i=0;$i<scalar(@ind_f);$i++){
			$countf=count_match(substr($l1,0,$lf[$i]),$ind_f[$i]);
			$countr=count_match(substr($l2,0,$lr[$i]),$ind_r[$i]);
			if(($forward_only&&$countf <= $allowmiss)||($countf <= $allowmiss&&$countr <= $allowmiss) ){ 
				$output{$i}=1;
			} else {
				$output{$i}=0;
			}
		}
	}

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
}



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

sub grab_adapters {
# finds (probable) adapter sequence by searching first 4000 lines of fastq files
# NOTE: is not primer specific - i.e. different multiplexed primer sets can't have unique adapters
# implementation of this is not too hard, convert adapter_match_x to a hash of arrays keyed by primer would work  - number of lines to search should probably be increased as well
	my %F;
	my %R;
	for(my $c=1;$c<scalar(@f4000);$c+=4) {
		my $l1 = $f4000[$c];
		my $l2 = $r4000[$c];
		my $x = length($l1);my $y=length($l2);
		if((length($l1)>30)&&(length($l2)>30)){
			my @adapter_match_f;
			my @adapter_match_r;

			for(my $i=0;$i<20;$i++) {
				foreach(@ind_f) {
					$adapter_match_f[$i]=count_match(substr($l1,$i,5),substr($_,0,5));
				}
				foreach(@ind_r) {
					$adapter_match_r[$i]=count_match(substr($l2,$i,5),substr($_,0,5));
				}			
			}

			my $idxMax_f = 0;
			my $idxMax_r = 0;

			$adapter_match_f[$idxMax_f] < $adapter_match_f[$_] or $idxMax_f = $_ for 1 .. $#adapter_match_f;
			$adapter_match_r[$idxMax_r] < $adapter_match_r[$_] or $idxMax_r = $_ for 1 .. $#adapter_match_r;

			$F{substr($l1,0,($idxMax_f))}++;
			$R{substr($l2,0,($idxMax_r))}++;
		}
	} 
	
	my $adapter_f = max_by { $F{$_} } keys %F;
	my $adapter_r = max_by { $R{$_} } keys %R;

	# add adapters to start of primer sequence
	@ind_f = map $adapter_f.$_, @ind_f;
	@ind_r = map $adapter_r.$_, @ind_r;

	# reset filehandlers to begining of file - good idea won't work in metabarcoding pipeline due to use of fifo input files
	# seek $R1, 0, 0;
	# seek $R2, 0, 0;

	print"Forward adapter $adapter_f\n";
	foreach(@ind_f) {
		print"Full forward    $_\n";
	}
	print"Reverse adapter $adapter_r\n";
	foreach(@ind_r) {
		print"Full reverse    $_\n";
	}
}
