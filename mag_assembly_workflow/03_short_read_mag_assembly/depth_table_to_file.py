#!/usr/bin/env python3

import sys
import os
import errno
import argparse
import hashlib
import subprocess
import shutil
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="[STRING] Input depth file", dest="input", required=True)
parser.add_argument("-p", help="[STRING] Output prefix", dest="prefix", default = ".")
args = parser.parse_args()

dat = dict()
lines = open( args.input, 'r')
first = True
for line in lines:
    stripline = line.strip()
    splitline = stripline.split('\t')
    if (first):
        headers = splitline
        for i,e in enumerate(headers):
            if headers[i] == "contigName" or headers[i].endswith(".bam"):
                dat[headers[i]] = []

        first = False
    else:
        line_array = splitline
        for i,e in enumerate(line_array):
            if headers[i] == "contigName" or headers[i].endswith(".bam"):
                dat[headers[i]].append(line_array[i])

output_list = open(args.prefix + "/" + re.sub('\.depth', '', os.path.basename(args.input)) + ".depth_list", 'w')
for key in dat:
#    sys.stdout.write("{}\n".format(key))
    if key != "contigName":
        base = re.sub('_nodup\.bam', '', key)
        outfile = args.prefix + "/" + base + ".depth"
        output = open(outfile, 'w')
        for i,e in enumerate(dat[key]):
            output.write("{}\t{}\n".format(dat["contigName"][i], dat[key][i]))
        output.close()
        output_list.write("{}\n".format(outfile))

output_list.close()
