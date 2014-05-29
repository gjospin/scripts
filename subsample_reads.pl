#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

#########################################
#
# Author:  Guillaume Jospin
# Date : 5/29/2014
#
# Purpose : Subsample PAIRED-END reads randomly down to a particular #
#
# Usage : subsample_reads.pl [options] <File1> (<File2>) <# of PAIRS of reads to keep> <output file>
#
# Output is written in INTERLEAVED format
#
# Option : --trim ##     Cuts the reads and quality scores to a specified length
#          --fastain     Forces the format of the input file(s) as FastA
#          --fastaout    Forces the format of the output file as FastA
#          --fastqin     Forces the format of the input file(s) as FastQ
#          --fastqout    Forces the format of the output file as FastQ
#          --interleaved Only 1 input file is submitted in interleaved format 
#          -h --help     Prints this help message
# Default : No trimming
#           Fastq IN
#           Fastq OUT
#
# Example : >subsample_reads.pl --trim 150 --interleaved --fastaout myfile.fastq 50 myfile_subsampled.fasta
# Example : >subsample_reads.pl myfile_1.fastq myfile_2.fastq 250 myfile_subsampled.fastq
#
############################################
my $usage = qq{
Usage : subsample_reads.pl [options] <File1> (<File2>) <# of PAIRS of reads to keep> <output file>

 Output is written in INTERLEAVED format

 Option : --trim ##     Cuts the reads and quality scores to a specified length
          --fastain     Force the format of the input file(s) as FastA
          --fastaout    Forces the format of the output file as FastA
          --fastqin     Forces the format of the input file(s) as FastQ
          --fastqout    Forces the format of the output file as FastQ
          --interleaved Only 1 input file is submitted in interleaved format
          --h --help    Prints this message
 Default : No trimming
           FastQ IN
           FastQ OUT
           Non interleaved input file

 Example : >subsample_reads.pl --trim 150 --interleaved --fastaout myfile.fastq 50 myfile_subsampled.fasta

 Example : >subsample_reads.pl myfile_1.fastq myfile_2.fastq 250 myfile_subsampled.fastq
 
 Example : >subsample_reads.pl --interleaved --fastain --fastaout myfile.fasta 50 myfile_subsampled.fasta

};

my $trim = 0;
my $interleaved = 0;
my $fastqin = 1;
my $fastqout = 1;
my $fastain = 0;
my $fastaout = 0;
my $help = 0;
my $h = 0;
GetOptions( "trim=i" => \$trim, # Length of resulting reads
	    "interleaved" => \$interleaved,
	    "fastain" => \$fastain,
	    "fastqin" => \$fastqin,
	    "fastaout"=> \$fastaout,
	    "fastqout"=> \$fastqout,
	    "help" => \$help,
	    "h" => \$h,
    );

$fastqin = 0 if $fastain;
$fastqout = 0 if $fastaout;

die "\n\n$usage" if @ARGV < 3 || $h || $help;

print STDERR "Detected options as ";
print STDERR "Non-" if $interleaved ==0;
my $char = $fastqin ? "Q" : "A";
print STDERR "interleaved file(s) in fast$char format\n";

die "\n\nERROR - Can't make up quality scores.\nFastA in and FastQ out are not compatible.\nPlease change the input and output formats\n$usage" if $fastain && $fastqout;
die "\n\nERROR - The --interleaved flag was used but the wrong number of arguments were submitted\n\n\n$usage" if @ARGV != 3 && $interleaved;
die "\n\nERROR - Wrong number of arguments when not using --interleaved\n\n$usage" if @ARGV!=4 && $interleaved==0;


my $file1 = shift;
my $file2 = shift unless $interleaved;
my $readpairs = shift;
my $outfile = shift;
$file1 =~ m/\/?([^\/]+)\.[^\.]$/;
my $core =$1;
#print "Core: $core\n$file\n";
#exit;
open(IN1,"< $file1") or die "Could not open $file1 for reading\n";
my @total_reads=();

# Reading the input file to collect all the reads in an array.
while(<IN1>){
    my $read_header = $_;
    my $read_core;
    if($fastqin){
	$read_header =~ m/^@([^\s]+)(\/\d)?/;
	$read_core = $1;
    }elsif($fastain){
	$read_header =~ m/^>([^\s]+)(\/\d)?/;
	$read_core = $1;
    }
    push(@total_reads, $read_core);
    <IN1>;
    <IN1> unless $fastain;
    <IN1> unless $fastain;
    if($interleaved){
	#need to skip the 2nd read in each pair
	<IN1>;
	<IN1>;
	<IN1> unless $fastain;
	<IN1> unless $fastain;
    }
}
close(IN1);
print STDERR "Finished reading the headers\n";
print STDERR "Shuffling the header array\n";
fisher_yates_shuffle(\@total_reads);
print STDERR "Finished shuffling the header array\n";
my $count1k=0;
my %keep1k=();
# saving the # of headers that we want to keep
for(my $i=0;$i< scalar(@total_reads);$i++){
    $count1k++;
    $keep1k{$total_reads[$i]}=1;
    if($count1k ==$readpairs){
	$i = scalar(@total_reads);
    }
}

my @read1=();
my @read2=();
my ($IN1,$IN2);
open($IN1,"< $file1");
unless($interleaved){
    open($IN2,"< $file2") or die "Could not open $file2 for reading\n";
}
open(OUT1k,">$outfile") or die "Could not open $outfile for printing\n";
print STDERR "Printing the subsampled reads to the output file\n";
while(<$IN1>){
    $read1[0] = $_;
    $read1[1] = <$IN1>;
    $read1[2] = <$IN1> if $fastqin;
    $read1[3] = <$IN1> if $fastqin;
    
    if($interleaved){
	$read2[0] = <$IN1>;
        $read2[1] = <$IN1>;
        $read2[2] = <$IN1> if $fastqin;
        $read2[3] = <$IN1> if $fastqin;
    }else{
	$read2[0] = <$IN2>;
	$read2[1] = <$IN2>;
	$read2[2] = <$IN2> if $fastqin;
	$read2[3] = <$IN2> if $fastqin;
    }
    
    if($trim > 0){
	#need to chomp in case the reads are shorter than the $trim value
	chomp($read1[1]);
	chomp($read1[3]) if $fastqin;
	chomp($read2[1]);
	chomp($read2[3]) if $fastqin;
	$read1[1] = substr($read1[1],0,$trim)."\n";
	$read1[3] = substr($read1[3],0,$trim)."\n" if $fastqin; 
	$read2[1] = substr($read2[1],0,$trim)."\n"; 
	$read2[3] = substr($read2[3],0,$trim)."\n" if $fastqin;
    }

    my $head;
    
    if($fastqin){
	$read1[0] =~ m/^@([^\s]+)(\/\d)?/;
	$head = $1;
    }elsif($fastain){
	$read1[0] =~ m/^>([^\s]+)(\/\d)?/;
    }
    if(exists $keep1k{$head}){
	if($fastqin && $fastaout){
	    # need to change the headers a bit if fastaout and fastqin
	    $read1[0] =~ s/^@/>/;
	    $read2[0] =~ s/^@/>/;
	}
	if($fastqout){
	    print OUT1k @read1;
	    print OUT1k @read2;
	}elsif($fastaout){
	    print OUT1k $read1[0].$read1[1];
	    print OUT1k $read2[0].$read2[1];
	}
    }
}
print STDERR "All done.\n";
1;

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

