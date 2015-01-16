# demultiplex_dualBC.pl

This script will demultiplex reads according to a barcode file. It will output
interleaved files (zipped), non interleaved files (I think not zippped),
trimmed files (adaptor contamination removed) and merged files (if the reads
overlap, they will get merged and modify the headers to be qiime ready). It
automatically detects the illumina nomenclature for index files.

## required software

* [Anthony Bolger's `trimmomatic.jar`](http://www.usadellab.org/cms/?page=trimmomatic)
* [FLASH, a tool for merging paired-end reads](http://ccb.jhu.edu/software/FLASH/)

## usage

perl demultiplex_dualBC.pl **<options> <illumina_directory> <mapping_file> <output_directory> <filename_core>**

**illumina_directory** The directory where the raw Illumina output lives.

**mapping_file** A tab delimited, QIIME-style mapping file, like so :

```
#SampleID       BarcodeSequence ReverseBarcode
SAMPLE1         AACGCTAA        GACTGGTT
SAMPLE2         AACGCTAA        TTAGCGTT
SAMPLE3         AACGCTAA        GTAGTCTT
```

**output_directory** The name of the directory where output will be written.

**filename_core** Some sort of identifier added to the files before the sample
names, which will look like `<core>_<samplename>.faa`.

### Options

**--trim-file** : Adapter sequences that need to be removed from the raw data in
the case of contamination.  An example is found in `adapter_all_16s.txt`.

**--trim-tool** The path for the trimming tool that will be used for removing
adaptor contamination. The file it is looking for is Anthony Bolger's
`trimmomatic.jar`, which can be obtained from the [Usadel Lab software
page](http://www.usadellab.org/cms/?page=trimmomatic).

**--reverse** For obscure reasons, this option is always required for demultiplexing
16S data.

**--q** [n] The minimum quality score that will be kept during the QC. 20 is the default.

**--read-len** : Read length. 250 or 300 these days.

**--frag-len** Length of the fragment you are targeting. For standard 16S you
will want something around 253.

**--min-overlap**  This will be used to merge the reads.  Won't merge reads
that overlap less than 10 residues (default 10).

**--max-overlap** This will be used to merge the reads. It won't
merge more than 70 residues.  This should be adjusted depending on fragment
length Vs read length. I use 120 often for this. Default is 70.

**--mismatch-ratio** Do not merge reads if the overlapped residues have a
mismatch ratio over this number. Default 0.25

**--skip-merge** Use this flag if you do not want to merge the reads.

**--phred** Phred alphabet used. 33 is default and what recent sequencing
machines use.

**--skip_interleaved** [1|0] Set this flag if you want to skip over the interleave file
generation.

**--no-mismatch** Skips over the mismatch gadgetery added in this script. By
default it will allow 1 mismatch per barcode half. Barcodes are set by 2
index reads.  It allows 1 mismatch for each index read.

**--min-read-length** [n] Discards reads that end up being shorter than the number
after the QC (and maybe merging).

