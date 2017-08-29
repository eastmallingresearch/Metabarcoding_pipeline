#! /usr/bin/perl -s

use Switch;

my $d ="unknown,0";
my $p ="unknown,0";
my $c ="unknown,0";
my $o ="unknown,0";
my $f ="unknown,0";
my $g ="unknown,0";
my $s ="unknown,0";

while(<>) {

	my @line = split(/\t/,$_);
	$line[1]=~s/\(|\)|:/\,/g;
	my @l = split(/\,/,$line[1]);
	test_it(($l[0],$l[1],$l[2]));
	test_it(($l[4],$l[5],$l[6]));
	test_it(($l[8],$l[9],$l[10]));
	test_it(($l[12],$l[13],$l[14]));
	test_it(($l[16],$l[17],$l[18]));
	test_it(($l[20],$l[21],$l[22]));
	test_it(($l[24],$l[25],$l[26]));
	chomp $d;
	chomp $p;
	chomp $c;
	chomp $o;
	chomp $f;
	chomp $g;
	chomp $s;
	print"$line[0],$d,$p,$c,$o,$f,$g,$s\n";
	$d ="unknown,0";
	$p ="unknown,0";
	$c ="unknown,0";
	$o ="unknown,0";
	$f ="unknown,0";
	$g ="unknown,0";
	$s ="unknown,0";

}

sub test_it {

	switch	($_[0]) {
		case "d" {$d="$_[1],$_[2]";}
		case "k" {$d="$_[1],$_[2]";}
		case "p" {$p="$_[1],$_[2]";}
		case "c" {$c="$_[1],$_[2]";}
		case "o" {$o="$_[1],$_[2]";}
		case "f" {$f="$_[1],$_[2]";}
		case "g" {$g="$_[1],$_[2]";}
		case "s" {$s="$_[1],$_[2]";}
	}
}
