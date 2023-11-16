#!/usr/bin/env python

################################################################
##################         Description        ##################
################################################################
# Author: Hardik I. Parikh
# Date: 13 Nov 2015
# Demultiplexing MiSeq Runs
#  
# This script performs demultiplexing of MiSeq paired-end reads
#
# Input parameters - 
# 
# SampleSheet - Tab-delimited file 
# 		Plate	Well	Index1	Spacer	Index2	LabID	SampleName
# 
# Index fastq file - I1.fastq.gz 
# 		Index1 file from MiSeq, actually contains index2
#		sequences, but read first by MiSeq, so labelled
#		index1
#
# Forward Reads fastq file - R1.fastq.gz (All forward reads)
# 
# Reverse Reads fastq file - R2.fastq.gz (All reverse reads)
# 
################################################################

import sys
import os
import errno
import re
import argparse
import HTSeq
from collections import namedtuple, defaultdict
from itertools import product
import warnings

sample_sheet_header = "Plate  Well  Index1  Spacer1  Index2 Spacer2  LabID  SampleID".split()

def check_seq(seq):
    return len(seq) == sum([seq.count(n) for n in "ACGT"])

def check_duplicates(lst):
    count = defaultdict(lambda: 0)
    for i in lst:
        count[i] += 1
    return [k for k, v in count.items() if v > 1]

Index = namedtuple("Index", "seq spacer trim")
def create_idx(seq, spacer):
    for f, s in zip(("index", "spacer"), (seq, spacer)):
        if not check_seq(s):
            raise ValueError("Bad sequence found for {}: {}".format(f, s))
    return Index(seq, spacer, len(seq+spacer))
UNK_IDX = create_idx("", "")

Sample = namedtuple("Sample", "plate well f_idx r_idx lab_id sample_id f_path r_path")
def create_sample(plate, well, f_idx, r_idx, lab_id, sample_id, outdir):
    filename = "_".join([sample_id, f_idx.seq, r_idx.seq])
    #filename = "_".join([lab_id, f_idx.seq, r_idx.seq])
    return Sample(
        plate, well,
        f_idx,
        r_idx,
        lab_id,
        sample_id,
        os.path.join(outdir,filename + "_R1.fastq"),
        os.path.join(outdir,filename + "_R2.fastq")
    )


def load_sample_sheet(sample_sheet_name, outdir, validate_header=True):
    fh = open(sample_sheet_name)

    header = fh.readline()
    if validate_header:
        if header.strip().split('\t') != sample_sheet_header:
            issues = ""
            header_lst = header.strip().split('\t')[:len(sample_sheet_header)]
            if len(header_lst) < len(sample_sheet_header):
                issues += "Expected at least {} columns, found {}".format(len(sample_sheet_header), len(header_lst))
            issues += "\n".join([ 
                "Read '{}' expected '{}'".format(h, s) 
                for h, s in zip(header_lst, sample_sheet_header) if h!=s
            ])
            raise ValueError("Sample sheet header is not as expeceted:\n" + issues)

    samples = []
    for i, line in enumerate(fh):
        fields = line.strip().split('\t')[:len(sample_sheet_header)]
        if len(fields) != len(sample_sheet_header):
            raise ValueError("Bad number of fields in sample sheet on line {}: got {}".format(i+1, len(fields)))
        plate, well, f_seq, f_space, r_seq, r_space, lab_id, sample_id = fields
        f_seq, f_space, r_seq, r_space = map(str.upper, (f_seq, f_space, r_seq, r_space))
        samples.append(create_sample(
            plate, well,
            create_idx(f_seq, f_space),
            create_idx(r_seq, r_space),
            lab_id,
            sample_id,
            outdir
        ))
    return samples

def validate_samples(samples):
    sample_idxs = [(s.f_idx.seq, s.r_idx.seq) for s in samples]
    if check_duplicates(sample_idxs):
        raise ValueError(
            "Duplicate sample indexes found: " + ", ".join(
                map("-".join, check_duplicates(sample_idxs))
            )
        )

    return samples

def validate_idxs(idxs):
    spacer_dict = defaultdict(set)
    for idx in idxs:
        spacer_dict[idx.seq].add(idx.spacer)
    issues = ""
    for idx, spacers in spacer_dict.items():
        if len(spacers) > 1:
            issues += "\nbarcode {}, spacers {}".format(idx, ", ".join(spacers))
    if issues:
        raise ValueError("Inconsistent spacers found for barcode(s):"+issues)

    l = len(idxs[0].seq)
    if not all([len(i.seq)==l for i in idxs[1:]]):
        raise ValueError("Inconsistent index length found")

    return l

def add_bc_mismatches(bc_map, filter_collisions=False):
    remove_x = set()
    for mm, bc in list(bc_map.items()):
        for i, n in product(range(len(mm)), [b"A", b"C", b"G", b"T"]):
            x = mm[:i] + n + mm[i+1:]
            if x in bc_map and bc_map[x] != bc:
                if filter_collisions:
                    warnings.warn("Barcode collision found, ignoring colliding mismatch")
                    remove_x.add(x)
                else:
                    raise ValueError("Barcode collision found: {} and {}".format(bc_map[x], bc))
            bc_map[x] = bc
    for x in remove_x:
        del bc_map[x]

# Function to write reads to output files
def printSamples(sample_reads_dict):
    for sample, sample_reads in sample_reads_dict.items():
        with open(sample.f_path, 'a') as outR1:
            with open(sample.r_path, 'a') as outR2:
                for r1, r2 in sample_reads:
                    r1.write_to_fastq_file(outR1)
                    r2.write_to_fastq_file(outR2)
    sample_reads_dict.clear()

SAMPLE_INFO = {'numReads' : 0}

#if __name__ == "__main__":
def main():
    # Argparse to get input from command line
    parser = argparse.ArgumentParser(description='This script performs demultiplexing of MiSeq paired-end reads')
    parser.add_argument('-ss', type=str, help=f"SampleSheetFile Tab-delimited (headings: {','.join(sample_sheet_header)})", required=True)
    parser.add_argument('-r1', type=str, help='R1 fastq file', required=True)
    parser.add_argument('-r2', type=str, help='R2 fastq file', required=True)
    parser.add_argument('-mm', type=int, help='Allow n barcode mismatches, default n=1', default=1)
    parser.add_argument('-fmmc', action='store_true', help='Filter colliding barcode mismatches, otherwise disallow collisions')
    parser.add_argument('-outdir', type=str, help='Output path, default current directory', default='.')
    args = parser.parse_args()

    # create output directory if it does not exist
    try:
        os.makedirs(args.outdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    samples = load_sample_sheet(args.ss, args.outdir)
    validate_samples(samples)
    sample_map = {(s.f_idx, s.r_idx): s for s in samples}
    sample_info = {s: dict(SAMPLE_INFO) for s in samples}

    f_idx_map = {s.f_idx.seq.encode(): s.f_idx for s in samples}
    r_idx_map = {s.r_idx.seq.encode(): s.r_idx for s in samples}

    # Index lengths
    forIndexLen = validate_idxs([s.f_idx for s in samples])
    revIndexLen = validate_idxs([s.r_idx for s in samples])

    for i in range(args.mm):
        add_bc_mismatches(f_idx_map, args.fmmc)
        add_bc_mismatches(r_idx_map, args.fmmc)
    print(f"Forward bc map size {len(f_idx_map)}")
    print(f"Reverse bc map size {len(r_idx_map)}")

    UNK_SAMPLE = Sample(
        "UNK", "UNK", UNK_IDX, UNK_IDX, "Undetermined", "Undetermined",
        os.path.join(args.outdir,"Undetermined_R1.fastq"),
        os.path.join(args.outdir,"Undetermined_R2.fastq")
    )
    sample_info[UNK_SAMPLE] = dict(SAMPLE_INFO)

    # Open Index1, R1, R2 fastq files - perform demultiplexing
    totReads = 0
    sample_reads_dict = defaultdict(list)

    r1fh = iter(HTSeq.FastqReader(args.r1))
    r2fh = iter(HTSeq.FastqReader(args.r2))

    #for i1,r1,r2 in zip(HTSeq.FastqReader(args.i1),HTSeq.FastqReader(args.r1),HTSeq.FastqReader(args.r2)):
    for r1, r2 in zip(r1fh, r2fh):

        totReads += 1

        if totReads % 10000 == 0:
            print("Processed - ", totReads)
            
        # Print Sample R1, R2 file; clear the dictionary 
        if totReads % 250000 == 0:
            printSamples(sample_reads_dict)

        # Map indexes
        f_idx = f_idx_map.get(r1.seq[:forIndexLen], 'unk')
        r_idx = r_idx_map.get(r2.seq[:revIndexLen], 'unk')
        sample = sample_map.get((f_idx, r_idx), UNK_SAMPLE)
        if f_idx == 'unk' or r_idx == 'unk':
            # unknown idxs, lump together to untrimmed set
            f_idx = r_idx = UNK_IDX
        elif sample == UNK_SAMPLE:
            #known idxs, unknown samples, place in new sample
            sample =create_sample(sample.plate, sample.well, f_idx, r_idx, sample.lab_id, sample.sample_id, args.outdir)
            sample_map[(f_idx, r_idx)] = sample
            samples.append(sample)
            sample_info[sample] = dict(SAMPLE_INFO)
        sample_info[sample]['numReads'] += 1

        # New read descriptions 
        forReadDesc = ":".join([r1.name, f_idx.seq, r_idx.seq, sample.sample_id])
        revReadDesc = ":".join([r2.name, f_idx.seq, r_idx.seq, sample.sample_id])

        # Create new HTSeq Objects for forward, reverse reads
        #remove forward index + spacer for r1 read
        newForRead = HTSeq.SequenceWithQualities(
            r1.seq[sample.f_idx.trim:], 
            forReadDesc, 
            r1.qualstr[sample.f_idx.trim:]
        )
        newRevRead = HTSeq.SequenceWithQualities(
            r2.seq[sample.r_idx.trim:], 
            revReadDesc, 
            r2.qualstr[sample.r_idx.trim:]
        )

        # Store 250,000 seq objects in sampleR1, sampleR2 dict
        sample_reads_dict[sample].append((newForRead, newRevRead))
            

    # For last reads < 250000
    print("Processing final reads ...")
    printSamples(sample_reads_dict)

    # Print Sample Info Dict
    samples.append(UNK_SAMPLE)
    with open(os.path.join(args.outdir, "demux_stats.txt"), "w") as ofh:
        for s in samples:
            percReads = (sample_info[s]["numReads"] / float(totReads) ) * 100
            outlist = [s.sample_id, str(sample_info[s]["numReads"]), "%.3f" % percReads, s.f_path, s.r_path]
            print("\t".join(outlist), file=ofh)

    print("Done!")

