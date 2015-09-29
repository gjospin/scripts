#!/usr/bin/perl -w
use strict;
use warnings;

# Author: Guillaume Jospin
# Created on Sept 25th 2015


my $usage = "m_to_z_bibtex_converter.pl only takes 1 argument, the name of the file to be converted\n\$perl m_to_z_bibtex_converter.pl <filename> > <filename2>\n\n";


die "Arguments not valid\n\n$usage\n" if scalar(@ARGV) != 1;

my $file = shift;


open(IN,$file);
while(<IN>){

    if($_ =~ /^keywords/){
	# do not print
    }elsif($_ =~ m/^mendeley-tags/){
	# change the tag and print the line 
	$_ =~ s/^mendeley-tags/keywords/;
	print $_;
    }else{
	# print normally
	print $_;
    }

}
close(IN);
