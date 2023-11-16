# demux_read_index
Demultiplex samples with in-read barcodes with spacers

# Requirements
Install provided conda environment, otherwise requires
- python3
- HTSeq

# Installation example
```
$ git clone git@github.com:cjowf/demux_read_index.git
$ cd demux_read_index
$ conda create -n demux_read_index -f environment.yaml
$ conda activate demux_read_index
$ pip install .
```

# Usage

```
$ demux_read_index.py -h
usage: demux_read_index.py [-h] -ss SS -r1 R1 -r2 R2 [-mm MM] [-fmmc] [-outdir OUTDIR]

This script performs demultiplexing of MiSeq paired-end reads

options:
  -h, --help      show this help message and exit
  -ss SS          SampleSheetFile Tab-delimited (headings:
                  Plate,Well,Index1,Spacer1,Index2,Spacer2,LabID,SampleID)
  -r1 R1          R1 fastq file
  -r2 R2          R2 fastq file
  -mm MM          Allow n barcode mismatches, default n=1
  -fmmc           Filter colliding barcode mismatches, otherwise disallow collisions
  -outdir OUTDIR  Output path, default current directory
```
## Details
- See test/sample_sheet.tsv for example of sample sheet format
- All output files written to OUTDIR
- Each sample in sample sheet will have two files, both with indexes and spacers **trimmed**
  - `<SampleID>_<Index1>_<Index2>_R1.fastq.gz`
  - `<SampleID>_<Index1>_<Index2>_R2.fastq.gz`
- Reads with one or both barcodes not found in sample sheet are **not trimmed**, found in files
  - `Undetermined_R1.fastq.gz`
  - `Undetermined_R2.fastq.gz`
- Reads with both barcodes in sample sheet, but not in a combination found in the sample sheet will be **trimmed**, found in files
  - `Undetermined_<Index1>_<Index2>_R1.fastq.gz`
  - `Undetermined_<Index1>_<Index2>_R2.fastq.gz`
- Demultiplexed statistics written to `OUTDIR/demux_stats.txt`

# Example workflow
1. bcl_convert for fastq generation.  For the sample sheet it should be a single sample entry for all 16S samples, barcodes could be:
     - Dummy barcodes, eg AAAAAAAAA, desired reads in "Undetermined" fastq files.
     - Primers with index cycle primer hybridization sites will have a defined barcode sequence. Using these for bcl_convert can be a quality filter and/or could allow multiplexing with other types of samples.  Desired reads in the sample_id fastq files.
2. demux_read_index.py for demultiplexing 16S samples.
      - All output files are primer and spacer trimmed, except Undetermined_R1 and Undetermined_R2.
      - Single threaded implementation, will take some time to run
      - No eta provided, will print '.' for every 10,000 reads
3. Check <out_dir>/demux_stats.txt for demuxed read counts. Note
      - Recognized barcodes in unrecognized combinations are put in Undetermined_<bc1>_<bc2>...
      - These are combinatorial barcodes; expect to see all possible combinations of forward and reverse barcodes, unrecognized combinations should have a small fraction of reads.
  
# TODO
- Automate tests
- Add index cycle demultiplexing
- Performance improvements
