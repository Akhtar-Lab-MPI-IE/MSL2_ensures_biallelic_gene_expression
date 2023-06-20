import os
import argparse
import sys
import gzip

parser = argparse.ArgumentParser(description='subset fastq file based on ixList')
parser.add_argument('--fq', type=str, required=True,
                    help='fastq file to process.')
parser.add_argument('--ix', type=str, required=True,
                    help='ixFile containing headers to include.')

args = parser.parse_args()

if not os.path.exists(args.fq):
    print("fastq file not found. exiting.")
    sys.exit()

if not os.path.exists(args.ix):
    print("ix file not found. exiting.")
    sys.exit()

ixSet = set()
counter = 0
with gzip.open(args.ix, 'rt') as f:
    for line in f:
        counter += 1
        ixSet.add(line.strip())
    if counter % 1000000 == 0:
        sys.stderr.write("Processed {} index lines.".format(counter))

with gzip.open(args.fq, 'rt') as f:
    fqCounter = 0
    counter = 0
    appendStatus = False
    for line in f:
        counter += 1
        fqCounter += 1
        if fqCounter % 1000000 == 0:
            sys.stderr.write("Processed {} lines in fastq files.".format(fqCounter))
        if counter == 1:
            head = line.strip().split(' ')[0].replace('@', '')
            if head in ixSet:
                appendStatus = True
                print(line.strip() )
            else:
                appendStatus = False
        elif appendStatus == True and counter != 4:
            print(line.strip() )
        elif appendStatus == False and counter != 4:
           continue
        elif counter == 4:
            if appendStatus == True:
                print(line.strip() )
                counter = 0
                appendStatus = False
            else:
                counter = 0
                appendStatus = False
        else:
            sys.stderr.write("Reached an exception state. Exit.")
            sys.stderr.write("counter {} , appendStatus {} ".format(counter, appendStatus))
            sys.exit()
 