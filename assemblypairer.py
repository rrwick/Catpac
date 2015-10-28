#!/usr/bin/env python


# Copyright 2015 Ryan Wick

# This file is part of AssemblyPairer.

# AssemblyPairer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.

# AssemblyPairer is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# AssemblyPairer.  If not, see <http:# www.gnu.org/licenses/>.


from __future__ import division
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import datetime
import shutil


def main():
    startTime = datetime.datetime.now()
    args = getArguments()

    print("\nAssemblyPairer\n--------------\n")

    # Load in the contigs from each assembly.
    print("Loading assemblies... ", end="")
    sys.stdout.flush()
    contigs1 = loadContigs(args.assembly1)
    contigs2 = loadContigs(args.assembly2)
    contigs1TotalLength = getTotalContigLength(contigs1)
    contigs2TotalLength = getTotalContigLength(contigs2)
    print("done\n")
    print("Loaded assembly 1: " + str(len(contigs1)) + " contigs, " + str(contigs1TotalLength) + " bp")
    print("Loaded assembly 2: " + str(len(contigs2)) + " contigs, " + str(contigs2TotalLength) + " bp\n")

    # Remove contigs below the length threshold.
    print("Filtering out contigs less than " + str(args.length) + " bp... ", end="")
    sys.stdout.flush()
    contigs1 = filterContigsByLength(contigs1, args.length)
    contigs2 = filterContigsByLength(contigs2, args.length)

    # Make a temporary directory for the alignment files.
    tempdir = os.getcwd() + '/temp'
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # Save the reduced contig sets to file.
    saveContigsToFile(contigs1, tempdir + "/contigs1.fasta")
    saveContigsToFile(contigs2, tempdir + "/contigs2.fasta")
    print("done\n")
    print("Filtered assembly 1: " + str(len(contigs1)) + " contigs, " + str(getTotalContigLength(contigs1)) + " bp")
    print("Filtered assembly 2: " + str(len(contigs2)) + " contigs, " + str(getTotalContigLength(contigs2)) + " bp\n")

    # Build a BLAST database using the first assembly.
    print("Building BLAST database... ", end="")
    sys.stdout.flush()
    makeblastdbCommand = ["makeblastdb", "-dbtype", "nucl", "-in", tempdir + "/contigs1.fasta"]
    p = subprocess.Popen(makeblastdbCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    print("done\n")

    # BLAST the second assembly against the first.
    print("Running BLAST search... ", end="")
    sys.stdout.flush()
    blastnCommand = ["blastn", "-db", tempdir + "/contigs1.fasta", "-query", tempdir + "/contigs2.fasta", "-outfmt", "6 length pident qseqid qstart qend qseq sseqid sstart send sseq mismatch gaps"]
    p = subprocess.Popen(blastnCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # Save the alignments in Python objects.
    alignmentStrings = out.splitlines()
    blastAlignments = []
    for alignmentString in alignmentStrings:
        alignment = BlastAlignment(alignmentString)
        blastAlignments.append(alignment)
    print("done\n")

    # Filter the alignments for length and identity.
    print("BLAST alignments before filtering:", len(blastAlignments))
    blastAlignments = filterBlastAlignments(blastAlignments, int(args.length), float(args.identity))
    print("BLAST alignments after filtering: ", len(blastAlignments))




    # CHECK TO ENSURE THERE ISN'T OVERLAP!!!!!
    # CHECK TO ENSURE THERE ISN'T OVERLAP!!!!!
    # CHECK TO ENSURE THERE ISN'T OVERLAP!!!!!
    # CHECK TO ENSURE THERE ISN'T OVERLAP!!!!!
    # CHECK TO ENSURE THERE ISN'T OVERLAP!!!!!




    # Display some summary information about the alignments.
    mismatches, gaps, length = totalMismatchesGapsAndLength(blastAlignments)
    print("\nTotal alignment mismatches:", mismatches)
    print("Total alignment gaps:      ", gaps)

    print("\nTotal alignment length:", length)
    contigs1Percent = 100.0 * length / contigs1TotalLength
    contigs2Percent = 100.0 * length / contigs2TotalLength
    print("{0:.3f}".format(contigs1Percent) + "% of assembly 1")
    print("{0:.3f}".format(contigs2Percent) + "% of assembly 2\n ")


    # Save the results to file
    print("Saving alignments to file... ", end="")
    sys.stdout.flush()
    saveAlignmentsToFile(blastAlignments, args.out1, True)
    saveAlignmentsToFile(blastAlignments, args.out2, False)
    print("done\n")

    # Delete the temporary files.
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    # Print a final message.
    endTime = datetime.datetime.now()
    duration = endTime - startTime
    print('Finished!')
    print('Total time to complete:', convertTimeDeltaToReadableString(duration))






def getArguments():
    parser = argparse.ArgumentParser(description='AssemblyPairer')
    parser.add_argument('assembly1', help='The first set of assembled contigs')
    parser.add_argument('assembly2', help='The second set of assembled contigs')
    parser.add_argument('out1', help='The filename for the first set of paired contigs')
    parser.add_argument('out2', help='The filename for the second set of paired contigs')
    parser.add_argument('-l', '--length', action='store', help='Minimum alignment length', default=1000)
    parser.add_argument('-i', '--identity', action='store', help='Minimum alignment percent identity', default=99.0)

    return parser.parse_args()



# This function takes a contig filename and returns a list of Contig objects.
def loadContigs(contigFilename):

    contigs = []

    contigFile = open(contigFilename, 'r')

    name = ''
    sequence = ''
    for line in contigFile:

        strippedLine = line.strip()

        # Skip empty lines.
        if len(strippedLine) == 0:
            continue

        # Header lines indicate the start of a new contig.
        if strippedLine[0] == '>':

            # If a contig is currently stored, save it now.
            if len(name) > 0:
                contig = Contig(name, sequence)
                contigs.append(contig)
                name = ''
                sequence = ''

            name = strippedLine[1:]

        # If not a header line, we assume this is a sequence line.
        else:
            sequence += strippedLine

    # Save the last contig.
    if len(name) > 0:
        contig = Contig(name, sequence)
        contigs.append(contig)

    return contigs



# This function takes a list of Contig objects and returns another list of
# Contig objects.  Only contigs with a length greater than or equal to the
# threshold will make it into the returned list.
def filterContigsByLength(contigs, lengthThreshold):

    lengthThreshold = int(lengthThreshold)
    filteredContigs = []

    for contig in contigs:
        if contig.length >= lengthThreshold:
            filteredContigs.append(contig)

    return filteredContigs



# This class holds a contig: its name, sequence and length.
class Contig:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.length = len(sequence)

    def __lt__(self, other):
        return self.length < other.length


# This class holds a BLAST alignment
class BlastAlignment:
    def __init__(self, blastString):
        blastStringParts = blastString.split("\t")

        self.length = int(blastStringParts[0])
        self.percentIdentity = float(blastStringParts[1])

        self.contig1Name = blastStringParts[2]
        self.contig1Start = int(blastStringParts[3])
        self.contig1End = int(blastStringParts[4])
        self.contig1Sequence = blastStringParts[5]

        self.contig2Name = blastStringParts[6]
        self.contig2Start = int(blastStringParts[7])
        self.contig2End = int(blastStringParts[8])
        self.contig2Sequence = blastStringParts[9]

        self.mismatches = int(blastStringParts[10])
        self.gaps = int(blastStringParts[11])



def convertTimeDeltaToReadableString(timeDelta):
    seconds = timeDelta.seconds
    hours = timeDelta.days * 24
    hours += seconds // 3600
    seconds = seconds % 3600
    minutes = seconds // 60
    seconds = seconds % 60
    seconds += timeDelta.microseconds / 1000000.0
    secondString = "{:.1f}".format(seconds)

    returnString = ""
    if hours > 0:
        return str(hours) + ' h, ' + str(minutes) + ' min, ' + secondString + ' s'
    if minutes > 0:
        return str(minutes) + ' min, ' + secondString + ' s'
    return secondString + ' s'



def saveContigsToFile(contigList, filename):
    outfile = open(filename, 'w')
    for contig in contigList:
        outfile.write('>' + contig.name + '\n')
        sequence = contig.sequence
        while len(sequence) > 60:
            outfile.write(sequence[0:60] + '\n')
            sequence = sequence[60:]
        outfile.write(sequence + '\n')



# If contig1 is True, then the alignments for the first set of contigs are
# saved to file.  If False, the alignments for the second set are saved.
def saveAlignmentsToFile(alignmentList, filename, contig1):
    outfile = open(filename, 'w')

    for alignment in alignmentList:

        alignmentName = ""
        if contig1:
            alignmentName += alignment.contig1Name
            alignmentName += "_" + str(alignment.contig1Start)
            alignmentName += "_to_" + str(alignment.contig1End)
        else:
            alignmentName += alignment.contig2Name
            alignmentName += "_" + str(alignment.contig2Start)
            alignmentName += "_to_" + str(alignment.contig2End)

        outfile.write('>' + alignmentName + '\n')

        sequence = ""
        if contig1:
            sequence = alignment.contig1Sequence
        else:
            sequence = alignment.contig2Sequence

        while len(sequence) > 60:
            outfile.write(sequence[0:60] + '\n')
            sequence = sequence[60:]
        outfile.write(sequence + '\n')



def getTotalContigLength(contigs):
    totalLength = 0
    for contig in contigs:
        totalLength += contig.length
    return totalLength



def filterBlastAlignments(alignments, minLength, minIdentity):
    filteredAlignments = []

    for alignment in alignments:
        if alignment.length >= minLength and alignment.percentIdentity >= minIdentity:
            filteredAlignments.append(alignment)

    return filteredAlignments


def totalMismatchesGapsAndLength(alignments):

    mismatches = 0
    gaps = 0
    length = 0

    for alignment in alignments:
        mismatches += alignment.mismatches
        gaps += alignment.gaps
        length += alignment.length

    return (mismatches, gaps, length)



# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()