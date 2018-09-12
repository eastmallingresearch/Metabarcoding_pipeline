#! /usr/bin/perl -s

use Switch;

my $d ="k__Fungi";
my $p ="p__unidentified";
my $c ="c__unidentified";
my $o ="o__unidentified";
my $f ="f__unidentified";
my $g ="g__unidentified";
my $s ="s__unidentified";

while(<>) {

	my @A = split(/\t/,$_);
	my @l = split(/;/,$A[1]);
	test_it($l[0]);
	test_it($l[1]);
	test_it($l[2]);
	test_it($l[3]);
	test_it($l[4]);
	test_it($l[5]);
	test_it($l[6]);
	chomp $d;
	chomp $p;
	chomp $c;
	chomp $o;
	chomp $f;
	chomp $g;
	chomp $s;

	if(q_i($s)) {
		if(q_i($g)) {
			if(q_i($f)) {
				if(q_i($o)) {
					if(q_i($c)) {
						$s=$d;
						$s=~s/.__/s__/;
					} else {
						$s=$c;
						$s=~s/.__/s__/;		
					}
				} else {
					$s=$o;
					$s=~s/.__/s__/;		
				}
			} else {
				$s=$f;
				$s=~s/.__/s__/;		
			}
		} else {
			$s=$g;
			$s=~s/.__/s__/;		
		}
	}	
	
	print"$A[0],$d,$p,$c,$o,$f,$g,$s\n";
	$d ="k__Fungi";
	$p ="p__unidentified";
	$c ="c__unidentified";
	$o ="o__unidentified";
	$f ="f__unidentified";
	$g ="g__unidentified";
	$s ="s__unidentified";

}


sub q_i {
	return 1 if ($_[0]=~/unidentified/);
	return 0;
}

sub test_it {
	
	my @s = split(/_/,$_[0]);
#print "$s[0]\n";
	switch	($s[0]) {
		case "k" {$d="k__$s[2]";}
		case "p" {$p="p__$s[2]";}
		case "c" {$c="c__$s[2]";}
		case "o" {$o="o__$s[2]";}
		case "f" {$f="f__$s[2]";}
		case "g" {$g="g__$s[2]";}
		case "s" {$s="s__$s[2]";}
	}
}

