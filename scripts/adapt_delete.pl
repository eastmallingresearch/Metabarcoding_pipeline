#!/usr/bin/perl -w -s

my $headers = do { local $/; <STDIN> };

#$headers =~ s/\R//g;
my ($inpFile) = $ARGV[0] || die "Please provide fastq file" ;
unless(open(INFILE, $inpFile) ) {
	print("Cannot open input file \"$inpFile\"\n\n");
     	exit;
}
my @fastq = <INFILE>;
chomp(@fastq);
close INFILE;

my $skipper = 0;
my $counter = 4;

foreach(@fastq) {
	if ($counter%4==0) {
		my @new = split(" ",$_,);
		substr($new[0], 0, 1) = "" ;
		if(index($headers,$new[0]) !=-1) {
			$skipper = 4;
		}
	}
	if(!$skipper) {
		print"$_\n";	
	}else {
		$skipper--;
	}
	$counter++;
}


