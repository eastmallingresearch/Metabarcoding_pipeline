#!/usr/bin/perl -s -w 
######################################################
# makes a biom table from sparse data
#
# ./otu_to_biom.pl row_file col_file data_file 
#	
#######################################################

open my $file, "<", $ARGV[0];
chomp(my @rows=<$file>);
close $file;

open $file, "<", $ARGV[1];
chomp(my @cols=<$file>);
close $file;

open $file, "<", $ARGV[2];
chomp(my @data=<$file>);
close $file;

my $shape=scalar @rows;
$shape.=",". scalar @cols;

print"{\n";
print"\t\"id\":\"otu_table.biom\",\n";
print"\t\"format\": \"Biological Observation Matrix 1.0\",\n";
print"\t\"format_url\": \"http:\/\/biom-format.org\",\n";
print"\t\"generated_by\": \"otu_to_biom\.pl_v1\",\n";
print"\t\"type\": \"OTU table\",\n";
my $dt= scalar localtime;
print"\t\"date\": \"$dt\",\n";
print"\t\"matrix_type\": \"sparse\",\n";
print"\t\"matrix_element_type\": \"float\",\n";
print"\t\"shape\": [$shape],\n";


my $last = pop @rows;
print"\t\"rows\":[\n";
foreach(@rows) {
	print"\t\t{\"id\":\"$_\",\"metadata\":null},\n";
}
print"\t\t{\"id\":\"$last\",\"metadata\":null}\n";
print"\t],\n";


$last = pop @cols;
print"\t\"columns\":[\n";	
foreach(@cols) {
	print"\t\t{\"id\":\"$_\",\"metadata\":null},\n";
}
print"\t\t{\"id\":\"$last\",\"metadata\":null}\n";
print"\t],\n";

$last = pop @data;
print"\t\"data\":[\n";	
foreach(@data) {
	print"\t\t[$_],\n";
}
print"\t\t[$last]\n";


print"\t]\n}";