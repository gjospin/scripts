#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
##########################################
#
# Author: Guillaume Jospin
# Workflow design: David Coil - Guillaume Jospin
#
# Date started: Jan 30th 2017
#
# Purpose: Confirm presence of isolated 16S sequences
#          in an OTU list from qiime 
#          Currently 1.9
#
# Input: List of 16S sequences in fastA format
#        Point to Qiime files
#             - rep_set.fna
#             - raw sequences used as qiime input
#             - OTU mapping file (which reads belong to what OTU
#
#
#
# Output: List of SeqHeader<tab>OTU#<tab>Metric<tab>Metric<tab>Taxonomy
#
#
#
##########################################


# Some parameters to adjust if the experimental 
# design or machine changes

my $threads = 14; # Uses 5 cpus when possible
my $fragsize = 253; # targets filters to a 253bp amplicon size
my $targets = 5; # Only look at a max of 5 OTUs

my $repset_file="";
my $qiime_input="";
my $qiime_otus="";
my $file_16s="";
my $working_dir="OTU_check_wd";
my $project_name="OTU_check_blastdb";
my $qiime_tax="";
my $debug=0;
GetOptions(
    "rep-set=s"        => \$repset_file,
    "qiime-input=s"           => \$qiime_input,
    "16s=s" => \$file_16s,
    "qiime-otus=s" => \$qiime_otus,
    "qiime-tax=s" => \$qiime_tax,
    "working-dir=s" => \$working_dir,
    "project-name=s" => \$project_name,
    "debug" => \$debug,
    );
my $usage = "Usage: perl OTU_check.pl --rep-set <file1> --qiime-input <file2> --16s <file3> --qiime-otus <file4>";
# Command line option checks
if($repset_file =~ m/^$/){
    die "Required argument --rep-set file was not set\n Terminating.\n\n$usage\n\n";
}
if($file_16s =~ m/^$/){
    die "Required argument --16s file was not set\n Terminating.\n\n$usage\n\n";
}
if($qiime_input =~ m/^$/){
    die "Required argument --qiime-input file was not set\n Terminating.\n\n$usage\n\n";
}
if($qiime_otus =~ m/^$/){
    die "Required argument --qiime-otus file was not set\n Terminating.\n\n$usage\n\n";
}if($qiime_tax =~ m/^$/){
    die "Required argument --qiime-tax file was not set\n Terminating.\n\n$usage\n\n";
}

# create working directory if needed
unless (-e "$working_dir"){
    `mkdir -p $working_dir`;
}
print STDERR "Created working directory $working_dir\n" if $debug;


# Parsing qiime taxonomy for later use
open(INTAX,$qiime_tax) or die "Could not read $qiime_tax for reading\n";
my %taxonomy=();
while(<INTAX>){
    chomp($_);
    my @line = split(/\t/,$_);
    $taxonomy{$line[0]}=$line[1];
}
close(INTAX);


# make a blastable DB from the rep_set file if needed
unless(-e "$working_dir/$project_name.blastdb.nhr"){
    `makeblastdb -out $working_dir/$project_name.blastdb -in $repset_file -input_type fasta -dbtype nucl`;
}
print STDERR "Finished creating the blast DB\n" if $debug;
# blast the 16S sequences against the repset file

unless(-e "$working_dir/$project_name.blastn"){
    `blastn -query $file_16s -db $working_dir/$project_name.blastdb -outfmt 6 -out $working_dir/$project_name.blastn -num_threads $threads -max_target_seqs $targets`;
}
print STDERR "Finished blasting the 16S queries against the repset provided\n" if $debug;
# Parse the output to identify the OTUs of interest for each sequence if any.

print STDERR "Parsing Blast output\n" if $debug;
my %blast_results=();
my %candidate_otus=();
open(INBLAST,"$working_dir/$project_name.blastn");
while(<INBLAST>){
    chomp($_);
#    skip comment lines
    next if $_ =~ m/^#/;

    my @line = split(/\t/,$_);
    $blast_results{$line[0]}{$line[1]}=$line[11] if $line[3] >= $fragsize-3;
    $candidate_otus{$line[1]}="" if $line[3] >= $fragsize-3;
}
close(INBLAST);
print STDERR "Finished parsing the blast output\n" if $debug;
print STDERR "Found ".scalar(keys %blast_results)." 16S sequences with results\n" if $debug;
print STDERR "Found ".scalar(keys %candidate_otus)." OTUs with results\n" if $debug;
print STDERR "Reading $qiime_otus file\n" if $debug;
# for each candidate OTU, extract the reads that were used to make the OTU.
# Find which reads make up what OTUs from the 
my %seqs_to_keep=();
open(IN,$qiime_otus) or die "Could not open $qiime_otus for reading.\nTerminating.\n";
while(<IN>){
    chomp($_);
#    next if ($_ =~ m/^/);
    my @line = split(/\t/,$_);
    next unless exists $candidate_otus{$line[0]};
    my $arraylen= scalar(@line);
    for(my $i = 1 ; $i < $arraylen; $i++){
	$seqs_to_keep{$line[$i]}=$line[0];
    }
}
close(IN);
print STDERR "Done reading qiime OTU file\n" if $debug;
# read the qiime input file for sequences and write them to their separate files
# Qiime input is assume to be in FastA format
print STDERR "Reading the raw qiime input sequences\n" if $debug;
open(INQIIMEINPUT,"zcat $qiime_input |") or die "Could not open $qiime_input for reading.\bTerminating.\n";

my %unique_seq=(); #record the uniqueness of sequences per OTUs
while(<INQIIMEINPUT>){
    chomp($_);
    if($_ =~ m/^>(\S+)\s?/){
	my $seq = <INQIIMEINPUT>;

	if(exists $seqs_to_keep{$1}){
	    #if it exists add the sequence to the hash for that OTU
	    unless(exists $unique_seq{$seqs_to_keep{$1}}{$seq}){
		#print STDERR "Seq comp: $seqs_to_keep{$1}\t$seq\n";
		$unique_seq{$seqs_to_keep{$1}}{$seq}= 1;
		$candidate_otus{$seqs_to_keep{$1}}.=$_."\n".$seq;
	    }
	}
    }else{
	# Shouldn't run into this but if we do, do nothing.
    }
}

close(INQIIMEINPUT);
print STDERR "Done reading Qiime raw input sequence file\n" if $debug;
# output 1 file per OTU containing all the files
print STDERR "Printing OTU files\n" if $debug;
foreach my $otu (keys %candidate_otus){
    open(OUTOTU,">$working_dir/$otu.fa") or die "Could not open $working_dir/$otu.fa for writting.\nTerminating.\n";
    print OUTOTU $candidate_otus{$otu};
    close(OUTOTU);
}
print STDERR "Adding query 16S to OTU files\n" if $debug;
# add the 16S queries to each candidate OTU that matched.
# First read in the 16S sequences from the input file
# match them to the OTUs
open(IN16S,$file_16s) or die "Could not open $file_16s for reading.\nTerminating.\n";
my $cur_seq="";
my $cur_head = "";
my %query16S=();
while(<IN16S>){
    chomp($_);
    if($_ =~ /^>(\S+)/){
	my $new_head=$1;
	$query16S{$cur_head}=$cur_seq;
	$cur_head = $new_head;
	$cur_seq = "";
    }else{
	$cur_seq .= $_;
    }
}
$query16S{$cur_head}=$cur_seq;
print STDERR "Found ".scalar(keys %query16S)." sequences in the 16S file\n" if $debug;
# Read blast results to assign the queries to the right OTU
foreach my $query(keys %blast_results){
    foreach my $otu (keys %{$blast_results{$query}}){
	#print STDERR "Processing $query\t$otu\n";
	open(OUTOTU,">>$working_dir/$otu.fa") or die "Could not open $working_dir/$otu.fa for appending.\nTerminating.\n";
	print OUTOTU ">$query.16S\n$query16S{$query}\n";
	close(OUTOTU);
    }
}


print STDERR "Aligning the queries and OTU repsets using cmalign\n" if $debug;
# Need to align the sequences using Muscle
my $realStart=0;
my $realEnd=0;
foreach my $otu(keys %candidate_otus){
#    next if $otu eq "974121";
    print STDERR "Processing ssualign for $otu\n" if $debug;
    `ssu-align -f --dna $working_dir/$otu.fa $working_dir/$otu.cmalign` unless -e "$working_dir/$otu.cmalign/$otu.cmalign.bacteria.stk";
    # Need to change the stockholm format into FastA
    my $alnIN = Bio::AlignIO->new(-file => "$working_dir/$otu.cmalign/$otu.cmalign.bacteria.stk",
	-format => "stockholm" );
    my $alnOUT = Bio::AlignIO->new(-file => ">$working_dir/$otu.cmalign.fa",
			       -format => "fasta" );
    while ( my $aln = $alnIN->next_aln ) {
	$alnOUT->write_aln($aln);
    }

    open(INALN,"$working_dir/$otu.cmalign.fa") or die "Could not open $working_dir/$otu.cmalign.fa for reading\n";
    my $i=0;
    my $starttotal=1;
    my $endtotal=0;
    my $curhead = "";
    my $curseq = "";
    while(<INALN>){
	chomp($_);
        next if ($_ =~ m/^$/);
        if($_ =~ m/^>(\S+)/){
            my $newhead=$1;
            if($curhead !~ m/^$/){
		if($curhead =~ /1-253/){
#		if($curhead =~ /1-903/){
		    $i++;
		    my $testseq = $curseq;
#           print STDERR "|$otu||".$seq->seq."|||\n";
		    if($testseq =~ /^([-\.]*)[^-\.]{10}.*[^-\.]{10}([-\.]*)$/){
			my $start=length($1);
			my $end=length($2);
#			print STDERR "Begin with ".$start." END with ".$end." -\n\n\n" if $otu eq "818181";
			$starttotal+=$start unless $start < 300;
			$endtotal+=$end unless $end < 300;
			last if $i >= 100;
		    }
		}
            }
            $curhead=$newhead;
            $curseq="";
        }else{
            $curseq .= $_;
        }
    }
    $i =1 if $i == 0;
    $realStart = int($starttotal/$i);
    $realEnd=int($endtotal/$i);
    print STDERR "$i LOOP\nStart: $realStart\tEnd: $realEnd\n";
    
#    exit;

    # Identifying a masking strategy
    open(INSSU,"$working_dir/$otu.cmalign.fa") or die "Could not open $working_dir/$otu.cmalign.fa";
    open(OUTMASK,">$working_dir/$otu.mask.fa") or die "Could not open $working_dir/$otu.mask.fa for writing";
    my $curseq ="";
    my $curhead = "";
    while(<INSSU>){
	chomp($_);
	next if ($_ =~ m/^$/);
	if($_ =~ m/^>(\S+)/){
	    my $newhead=$1;
	    if($curhead !~ m/^$/){
#		print STDERR "String leng: ".length($curseq)."\tSTART: $realStart\tEND: $realEnd\n$curseq\n" if $otu eq "818181";
		my $substr = substr($curseq,$realStart,-$realEnd);
#		print STDERR "$substr\nLen:".length($substr)."\n" if $otu eq "818181";
#		exit if $otu eq "818181";
		$substr =~ s/[\.-]//g;
		print OUTMASK ">$curhead\n$substr\n";
	    }
	    $curhead=$newhead;
	    $curseq="";
	}else{
	    $curseq .= $_;
	}
    }
    my $substr = substr($curseq,$realStart,-$realEnd);
    $substr =~ s/[\.-]//g;
    print OUTMASK ">$curhead\n$substr\n";
    close(INSSU);
    $starttotal = 0;
    $endtotal = 0;
}

# Blast all V all for each OTU to get % identitfy


print STDERR "Finished aligning and masking sequences\n" if $debug;
print STDERR "Starting blastn all V all for each OTU\n" if $debug;
#exit;
# Opening file handle for all stats in 1 file

open(ALLSTATS,">$working_dir/allstats.txt");
print ALLSTATS "\tOTUlabel\t99pct\t97pct\t90pct\n";
foreach my $otu(keys %candidate_otus){
    `blastn -query $working_dir/$otu.mask.fa -subject $working_dir/$otu.mask.fa -num_threads $threads -outfmt 6 -out $working_dir/$otu.blastn` unless -e "$working_dir/$otu.blastn";

# Parse all V all blast for each OTU and output
# % identity summary stats for each sequence.
    
    open(INBLASTN, "$working_dir/$otu.blastn") or die "Could not open $working_dir/$otu.blastn for reading\n";
    open(OUTSTATS,">$working_dir/$otu.stats") or die "Could not open $working_dir/$otu.stats";
    my %totalmatches=();
    my %higher90=();
    my %higher97=();
    my %higher99=();
    while(<INBLASTN>){
	# skip comments
	next if $_ =~ m/^#/;
	chomp($_);
	my @line= split(/\t/,$_);
	
	my $query = $line[0];
	my $subject = $line[1];
	# don't need to compare self identity
	next if $query eq $subject;
	# Don't want to compare the input 
	# 16S sequences to each other
	if($query =~ /\.16S$/ && $subject =~m/\.16S$/){
	    next;
	}
	$totalmatches{$line[0]}=0 unless defined $totalmatches{$line[0]};
	$higher90{$line[0]}=0 unless defined $higher90{$line[0]};
	$higher97{$line[0]}=0 unless defined $higher97{$line[0]};
	$higher99{$line[0]}=0 unless defined $higher99{$line[0]};
	if($line[2] > 99){
	    $higher99{$query}++;
	}
	if($line[2] > 97){
	    $higher97{$query}++;
	}
	if($line[2] > 90){
	    $higher90{$query}++;
	}
	$totalmatches{$query}++;
    }
    close(INBLASTN);
# header line formatted as:
# <tabl>OTUlabel<tab>statsheaders<tab>...<tab>statsheader\n

# Stats line fomartted as follow:
# SeqID<tab>matchingOTU<tab>stats<tab>...<tab>stats\n

# What pct of the sequences in the OTU does our queries
# match to at a 99% (97% and 90%) identity or higher
    print OUTSTATS "\tOTUlabel\t99pct\t97pct\t90pct\n";
    foreach my $query(keys %totalmatches){
	next if $query !~ /.16S/;
	print OUTSTATS $query."\t";
	my $pct90 =  sprintf("%.2f", 100*$higher90{$query}/$totalmatches{$query});
	my $pct97 = sprintf("%.2f",100*$higher97{$query}/$totalmatches{$query});
	my $pct99 = sprintf("%.2f",100*$higher99{$query}/$totalmatches{$query});
	print OUTSTATS $otu."\t".$pct99."\t".$pct97."\t".$pct90;
	print OUTSTATS "\t".$taxonomy{$otu} if exists $taxonomy{$otu};
	print OUTSTATS "\n";


	print ALLSTATS $query."\t";
	print ALLSTATS $otu."\t".$pct99."\t".$pct97."\t".$pct90;
        print ALLSTATS "\t".$taxonomy{$otu} if exists $taxonomy{$otu};
        print ALLSTATS "\n";

    }
    
}
close(ALLSTATS);

