#!/usr/bin/env python3

"""Convert a GenBank file into FFN format.

Usage:
    genbank_to_ffn.py <genbank_file>
"""
import sys
import os
import argparse
import textwrap

from Bio import SeqIO
from Bio import Seq

def collapse_fasta(fasta_file, input_dir, output_dir):
#    sys.stdout.write("[STATUS] Collapsing {}".format(base))
    base = os.path.splitext(os.path.basename(fasta_file))[0]
    sys.stdout.write("[STATUS]\tGenome: {}\n".format(base))
    input_handle = open(os.path.join(input_dir, fasta_file), "r")
    output_handle = open(os.path.join(output_dir, base + "_collapsed.fna"), "w")
    output_handle.write(">{}\n".format(base))
    sequence = str()
    first = True
    for record in SeqIO.parse(input_handle, "fasta") :
#        print "Dealing with GenBank record %s" % seq_record.id
        sys.stdout.write("[STATUS]\t\tProcessing contig {}\n".format(record.id))
        if first == True:
            sequence = str(record.seq)
            first = False
        else:
            sequence = sequence + "NNNNNNNNNN" + str(record.seq)
    output_handle.write("{}\n".format(textwrap.fill(sequence, width = 60)))
    input_handle.close()
    output_handle.close()
#        for record in seq_record.features:
#            if seq_feature.type=="CDS":
#                output_handle.write(">{} {}\n{}\n".format(seq_feature.qualifiers['locus_tag'][0], seq_record.name, seq_feature.extract(seq_record.seq)))
#    input_handle.close()
#    output_handle.close()

#    out_file = "%s.fna" % os.path.splitext(gb_file)[0]
#    with open(out_file, "w") as out_handle:
#        GFF.write(SeqIO.parse(gb_file, "genbank"), out_handle)

parser = argparse.ArgumentParser(description='Process mash results for kallisto index creation. Run script from within directory containing all genome fasta files.')
parser.add_argument('-d', '--directory', dest="input_dir", default="./", help="How many strains of each species to keep for the quantification step")
parser.add_argument('-o', '--output', dest="output_dir", default="collapsed", help="Directory to put files for kallisto index creation. Default is moving them one directory up.")
parser.add_argument('-e', '--extension', dest="ext", default=".fna", help="If set, lists files that would be moved, but does not create or move any files.")
args = parser.parse_args()

if __name__ == "__main__":
    files = os.listdir(args.input_dir)
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    for f in files:
#        sys.stdout.write("{}".format(os.path.splitext(f)[1]))
        if os.path.splitext(f)[1] == args.ext:
            collapse_fasta(f, args.input_dir, args.output_dir)


#Short version:
#SeqIO.write(SeqIO.parse(input_handle, "genbank"), output_handle, "fasta")
