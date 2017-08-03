#!/usr/bin/perl -s -w
#######################################
# merges biom table with taxonomy data
#
# ./biom_maker.pl taxafile biom_table
#	
######################################

my ($taxaFile) = $ARGV[0] || die "Please provide taxa file" ; 
unless(open(INFILE, $taxaFile) ) {
	print("Cannot open taxa file \"$taxaFile\"\n\n");
     	exit;
}
my @taxa = <INFILE>;
chomp(@taxa);
close INFILE;
my %seq=();
foreach(@taxa) {
 $_ =~ s/"//g;
 @t = split(/,/);
 $seq{$t[0]}="{\"taxonomy\": [\"$t[1]\", \"$t[3]\", \"$t[5]\", \"$t[7]\", \"$t[9]\", \"$t[11]\", \"$t[13]\", $t[2],$t[4],$t[6],$t[8],$t[10],$t[12],$t[14]]}";

#$seq{$t[0]}="{\"taxonomy\": [\"k__$t[1]\", \"p__$t[3]\", \"c__$t[5]\", \"o__$t[7]\", \"f__$t[9]\", \"g__$t[11]\", \"s__$t[13]\"]}";
}  


my ($inpFile) = $ARGV[1] || die "Please provide biom table" ;
unless(open(INFILE, $inpFile) ) {
	print("Cannot open biom table \"$inpFile\"\n\n");
     	exit;
}
my @biom = <INFILE>;
chomp(@biom);
close INFILE;

foreach(@biom) {
  if(($OTU)=($_=~/.*(OTU\d*)",/)) {
    #$OTU=($_=~/.*(OTU\d*)",/)
    $_=~s/null/$seq{$OTU}/;
  }
  print"$_\n";
}