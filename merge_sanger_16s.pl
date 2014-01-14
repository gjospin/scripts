#!/usr/bin/perl -w
use strict;
use warnings;

###################################
#
# Author : Guillaume Jospin
# Date : Jan 13th 2014
#
# Usage : >perl merge_sanger_16s.pl file1.fastq file2.fastq
#
# Pre-reqs : Bio::SeqIO (Bio perl) http://search.cpan.org/~cjfields/BioPerl-1.6.901/Bio/SeqIO.pm
#            Muscle - http://www.drive5.com/muscle/
# Description : Reads in 2 ab1 files, forward and reverse reads in ab1 format.
#               Outputs fastq files for each read
#               Qtrim the reads according to a default score of 20
#               Aligns the reads using Muscle
#               Merge the reads from the Muscle output into 1 sequence.
#
#     Assumes that forward is the first file submitted and reverse is the second.
#     Assumes that the reads overlap after the quality trimming.
#     Assumes that muscle is accessible from the user's path.
#     Assumes quality scores of in the PHRED33 format/scale
#
#
# WARNING: This script was designed for use with primers 1391 R and 27 F
#
###################################

die "Number of argument is incorrect. Need excatly 2 arguments. 1 forward read and 1 reverse read in that order.\n" if @ARGV != 2 ;

my $file1 = shift;
my $file2 = shift;
my $qtrim_threshold = 20;
my $fragment_length = 1350;

$file1 =~ /([^\/]*).fa?s?t?q$/;
my $id = $1;
my @read1 = ();
my @read2 = ();
open(IN1,$file1);
@read1 = <IN1>;
close(IN1);
open(IN2,$file2);
@read2 = <IN2>;
close(IN2);

# change the IDs to make them simple and consisten over the 2 reads
$read1[0] = $id."\n";
$read2[0] = $id."\n";
# qtrim the reads from the end
qtrim_read(read=>\@read1,quality => $qtrim_threshold, readtype => "phred33");
qtrim_read(read=>\@read2,quality => $qtrim_threshold, readtype => "phred33");
# make sure there are no trailing new line characters
for(my $i = 0 ;$i < scalar(@read1);$i++){
    chomp($read1[$i]);
    chomp($read2[$i]);
}

# Qtrim the reads from the front
# Reversing to use the same algorithm.
$read1[1] = reverse($read1[1]);
$read1[3] = reverse($read1[3]);
$read2[1] = reverse($read2[1]);
$read2[3] = reverse($read2[3]);
qtrim_read(read=>\@read1,quality => $qtrim_threshold, readtype => "phred33");
qtrim_read(read=>\@read2,quality => $qtrim_threshold, readtype => "phred33");
for(my $i = 0 ;$i < scalar(@read1);$i++){
    chomp($read1[$i]);
    chomp($read2[$i]);
}
# reversing to get back to the original arrangment.
$read1[1] = reverse($read1[1]);
$read1[3] = reverse($read1[3]);
$read2[1] = reverse($read2[1]);
$read2[3] = reverse($read2[3]);

# Reverse complementing the reverse read.
my $rev = reverse($read2[1]);
$rev =~ tr/ACGTacgt/TGCAtgca/;
$read2[1] = "$rev";
$read2[3] = reverse($read2[3]);

if(length($read1[1]) + length($read2[1]) < 1350){
    die "Not enough data to overlap confidently.\nExiting.\n";
}



# Print temp file for alignment with Muscle.
open(ALNOUT,">temp_aln.fa");
print ALNOUT ">forward\n$read1[1]\n";
print ALNOUT ">reverse\n$read2[1]\n";
close(ALNOUT);

# aligning the two sequences.
`muscle -in temp_aln.fa -out temp_aln.aln > /dev/null 2> /dev/null`;

# read in the alignment.
open(INALN, "temp_aln.aln");
my $cur_seq = "";
my $cur_id = "";
my $seq1 ="";
while(<INALN>){
    chomp($_);
    if($_ =~ m/^>(\S+)/){
	if($cur_id eq "forward"){
	    $seq1 = $cur_seq;
	    $cur_seq = "";
	}
	$cur_id = $1;
    }else{
	$cur_seq .= $_;
    }
}
my $seq2 = $cur_seq;
close(INALN);
`rm temp_aln.fa temp_aln.aln`;
my $merged_seq = merge_alignment(read1=>$seq1, read2=>$seq2, quality1 => $read1[3], quality2=> $read2[3]);


print STDERR "Final merged sequence can be found in $id"."_merged.fasta\n";
open(OUT,">$id"."_merged.fasta");
print OUT ">$id"."_merged\n$merged_seq\n";
close(OUT);


=head2 merge_alignment

Merges 2 aligned reads using Quality scores to decide tie breakers.

Input:  Read1, Read2, Quality1, Quality2
Output: Merged sequence.

=cut

sub merge_alignment{
    my %args = @_;
    my $forward = $args{read1};
    my $reverse = $args{read2};
    my $qseq_1 = $args{quality1};
    my $qseq_2 = $args{quality2};
    chomp($forward);
    chomp($reverse);
    chomp($qseq_1);
    chomp($qseq_2);
    my $q_1=0;
    my $q_2=0;
    my @qseq1 = split(//, $qseq_1);
    my @qseq2 = split(//, $qseq_2);
    my @seq1 = split( //, $forward );
    my @seq2 = split( //, $reverse );
    my $result_seq = "";
    my $len = scalar(@seq1);
    my $special=0;
    for ( my $i = 0; $i < $len ; $i++ ) {
	if ( $seq1[$i] eq $seq2[$i] ) {
	    $result_seq .= $seq1[$i];
	} elsif ( $seq1[$i] ne $seq2[$i] && $seq1[$i] =~ /[-\.]/) {
	    $result_seq .= $seq2[$i];
	} elsif ( $seq1[$i] ne $seq2[$i] && $seq2[$i] =~ /[-\.]/ ) {
	    $result_seq .= $seq1[$i];
	} else {
	    #print STDERR  "FOUND A SPECIAL CASE $seq1[$i] $seq2[$i]\n";
	    #print STDERR  "FOUND A SPECIAL CASE $qseq1[$q_1] $qseq2[$q_2]\n";
	    #print STDERR  "FOUND A SPECIAL CASE $q_1 $q_2\n";
	    #print STDERR  "FOUND A SPECIAL CASE ".ord($qseq1[$q_1])." ".ord($qseq2[$q_2])."\n";
	    $special++;
	    if(ord($qseq1[$q_1]) > ord($qseq2[$q_2])){
		# quality is better for seq1, use that residue
		$result_seq .= $seq1[$i];
	    }elsif(ord($qseq1[$q_1]) < ord($qseq2[$q_2])){
		# quality is better for seq2, use that residue
		$result_seq .= $seq2[$i];
	    }else{
		#qual scores are equal, use the one from seq 1
		$result_seq .= $seq1[$i];
	    }
	}
	$q_1++ if $seq1[$i] =~ m/[A-Za-z]/;
	$q_2++ if $seq2[$i] =~ m/[A-Za-z]/;
	#print "$result_seq\n";
    }
    print STDERR "Found $special conflicting case(s) during merging of $len residues\n";
    return $result_seq;
    
}

=head2 qtrim_read

trims a fastq read to a particular quality score using Heng Li's algorithm from bwa.
code based on SGA's implementation.

=cut

sub qtrim_read {
    my %args     = @_;
    my $read     = $args{read};
    my $q        = $args{quality};
    my $readtype = $args{readtype};

    $q += 33 if $readtype eq "phred33";
    $q += 64 if $readtype eq "phred64";

        # Perform a soft-clipping of the sequence by removing low quality bases from the
        # 3' end using Heng Li's algorithm from bwa

    my $seq = @$read[1];
    chomp $seq;
    my $qq = @$read[3];
    chomp $qq;
    my @qual          = split( //, $qq );
    my $endpoint      = 0;                  # not inclusive
    my $max           = 0;
    my $i             = length($seq) - 1;
    my $terminalScore = ord( $qual[$i] );

        # Only perform soft-clipping if the last base has qual less than $q
    return if ( $terminalScore >= $q );

    my $subSum = 0;
    while ( $i >= 0 ) {
	my $ps    = ord( $qual[$i] );
	my $score = $q - $ps;
	$subSum += $score;
	if ( $subSum > $max ) {
	    $max      = $subSum;
	    $endpoint = $i;
	}
	$i--;
    }

        # Clip the read
    @$read[1] = substr( $seq, 0, $endpoint )."\n";
    @$read[3] = substr( @$read[3], 0, $endpoint )."\n";
}
