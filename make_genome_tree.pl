#!/usr/bin/perl -w
use strict;
use warnings;

################################
#
# make_genome_tree.pl
#
# Author: Guillaume Jospin
# 
# usage: >perl make_genome_tree_vibrio.pl <genomedirectory> <Name of output trees>
# 
# WARNING : This scripts looks for genome files ending in .fasta 
#           in the genome directory used as the first argument
# 
# ****   USE AT YOUR OWN RISK    ****
#    paths for PhyloSift were hard coded - need to change to match your system.
#              PhyloSift config files NOT included
#
#   If intermediate files exists, thing will NOT get regenerated
#   Clean up your genome directory when restarting a run appropriately.
#
# Workflow outline
# 
# 0 - Concatenate the sequences into a single sequence
# 1 - Run Phylosift
# 2 - Extract the aligned concat sequences (nucl AND prot) AND 16S
#   - Rename the sequences according to file name 
# 3 - Make a tree (using FastTree or RaxML)
#
###################################3


my $dir = shift;
my $core_name = shift;
my @files = <$dir/*.fasta>;

foreach my $file (@files){

    $file =~ m/([^\/]+).fasta/;
    my $core = $1;
    print STDERR "concatenating $file\n";
    next if -e "$dir/$core"."_concat.fa";
    open(IN,$file);
    open(OUT,">$dir/$core"."_concat.fa");
    print OUT ">$core\n";
    while(<IN>){
        chomp($_);
        if($_ =~ m/^>/){

        }else{
            print OUT $_;
        }
    }
    close(OUT);
    close(IN);

}

@files = <$dir/*_concat.fa>;
foreach my $file(@files){

    $file =~ m/\/([^\/]+)_concat.fa/;
    my $core = $1;
    print STDERR "Running PhyloSift on $core\n";
    next if -e "$dir/PS_temp/$core/alignDir/concat.codon.updated.1.fasta";
    `mkdir -p $dir/PS_temp/$core`;
    
    `phylosift search --config /Users/gjospin/Documents/Berg54_Christenella/wholegenome_tree/phylosiftrc --isolate --besthit --output $dir/PS_temp/$core $file`;
    `phylosift align -f --config /Users/gjospin/Documents/Berg54_Christenella/wholegenome_tree/phylosiftrc --output $dir/PS_temp/$core $file`;
}

print STDERR "Gathering concat nucleotide sequences\n";
@files = <$dir/PS_temp/*/alignDir/concat.codon.updated.1.fasta>;
open(OUT,">$dir/$core_name.concat.nucl.fa");
foreach my $file(@files){
    print "Processing $file\n";
    $file =~ m/\/([^\/]+)\/alignDir/;
    my $index = $1;
    open(IN,$file);
    while(<IN>){
        chomp($_);
        if ($_=~ m/^>/){
            print OUT ">$index"."\n";
        }else{
            print OUT $_."\n";
        }
    }
}
close(OUT);

print STDERR "Gathering concat AA sequences\n";
@files = <$dir/PS_temp/*/alignDir/concat.updated.1.fasta>;
open(OUT,">$dir/$core_name.concat.AA.fa");
foreach my $file(@files){
    print "Processing $file\n";
    $file =~ m/\/([^\/]+)\/alignDir/;
    my $index = $1;
    open(IN,$file);
    while(<IN>){
        chomp($_);
        if ($_=~ m/^>/){
            print OUT ">$index"."\n";
        }else{
            print OUT $_."\n";
        }
    }
}
close(OUT);

print STDERR "Gathering 16S sequences\n";

@files = <$dir/PS_temp/*/alignDir/16s_reps_bac.long.1.fasta>;
my %strains = ();
open(OUT,">$dir/$core_name.16S.fa");
foreach my $file(@files){
    print "Processing $file\n";
    $file =~ m/\/([^\/]+)\/alignDir/;
    my $index = $1;
   
   
    open(IN,$file);
    while(<IN>){
        chomp($_);
        if ($_=~ m/^>/){
            next if exists $strains{$index};
            print OUT ">$index"."\n";
            $strains{$index}=0;
        }else{
            next if $strains{$index}==1;
            print OUT $_."\n";
            $strains{$index}++;
        }
    }
}
close(OUT);

print STDERR "Building concat AA tree\n";
`FastTree $dir/$core_name.concat.AA.fa > $dir/$core_name.AA.tre` unless -e "$dir/$core_name.AA.tre";
print STDERR "Building concat nucleotide tree\n";
`FastTree -nt -gtr < $dir/$core_name.concat.nucl.fa > $dir/$core_name.nucl.tre` unless -e "$dir/$core_name.nucl.tre";
print STDERR "Building 16S tree\n";
`FastTree -nt -gtr < $dir/$core_name.16S.fa > $dir/$core_name.16S.tre` unless -e "$dir/$core_name.16S.tre";
