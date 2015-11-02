#!/usr/bin/env python


# Copyright 2015 Ryan Wick

# This file is part of Catpac.

# Catpac is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.

# Catpac is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# Catpac.  If not, see <http:# www.gnu.org/licenses/>.


from __future__ import division
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import datetime
import shutil
import itertools


def main():
    startTime = datetime.datetime.now()
    args = getArguments()

    print("\nCatpac: a Contig Alignment Tool for Pairwise Assembly Comparison\n----------------------------------------------------------------\n")

    # Load in the contigs from each assembly.
    print("Loading assemblies... ", end="")
    sys.stdout.flush()
    contigs1 = loadContigs(args.assembly1)
    contigs2 = loadContigs(args.assembly2)
    contigs1TotalLength = getTotalContigLength(contigs1)
    contigs2TotalLength = getTotalContigLength(contigs2)
    contigs1MedianReadDepth, contigs1MedianAbsoluteDeviation = getMedianReadDepthByBaseAndMedianAbsoluteDeviation(contigs1)
    contigs2MedianReadDepth, contigs2MedianAbsoluteDeviation = getMedianReadDepthByBaseAndMedianAbsoluteDeviation(contigs2)
    calculateRelativeDepthAndZScore(contigs1, contigs1MedianReadDepth, contigs1MedianAbsoluteDeviation)
    calculateRelativeDepthAndZScore(contigs2, contigs2MedianReadDepth, contigs2MedianAbsoluteDeviation)
    print("done\n")
    print("Loaded assembly 1: ")
    print("   " + str(len(contigs1)) + " contigs, " + str(contigs1TotalLength) + " bp")
    print("   median read depth (by base): " + str(contigs1MedianReadDepth))
    print("   median absolute deviation:   " + str(contigs1MedianAbsoluteDeviation))
    print("\nLoaded assembly 2: ")
    print("   " + str(len(contigs2)) + " contigs, " + str(contigs2TotalLength) + " bp")
    print("   median read depth (by base): " + str(contigs2MedianReadDepth))
    print("   median absolute deviation:   " + str(contigs2MedianAbsoluteDeviation) + "\n")

    # Build dictionaries for the contigs, so we can later use a contig's name
    # to get the rest of the contig's details.
    contigs1Dict = {}
    for contig in contigs1:
        contigs1Dict[contig.fullname] = contig
    contigs2Dict = {}
    for contig in contigs2:
        contigs2Dict[contig.fullname] = contig

    # Make a temporary directory for the alignment files.
    tempdir = os.getcwd() + '/temp'
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # Remove contigs below the length threshold.
    if int(args.length) > 0:
        print("Filtering out contigs less than " + str(args.length) + " bp... ", end="")
        sys.stdout.flush()
        contigs1 = filterContigsByLength(contigs1, args.length)
        contigs2 = filterContigsByLength(contigs2, args.length)
        print("done")
        print("   Filtered assembly 1: " + str(len(contigs1)) + " contigs, " + str(getTotalContigLength(contigs1)) + " bp")
        print("   Filtered assembly 2: " + str(len(contigs2)) + " contigs, " + str(getTotalContigLength(contigs2)) + " bp\n")

    # Remove contigs outside the relative read depth threshold.
    if float(args.minreaddepth) > 0.0 or float(args.maxreaddepth) < float("inf"):
        if float(args.minreaddepth) > 0.0 and float(args.maxreaddepth) < float("inf"):
            print("Filtering out contigs with a relative read depth less than " + str(args.minreaddepth) + " or greater than " + str(args.maxreaddepth) + "... ", end="")
        elif float(args.minreaddepth) > 0.0:
            print("Filtering out contigs with a relative read depth less than " + str(args.minreaddepth) + "... ", end="")
        elif float(args.maxreaddepth) < float("inf"):
            print("Filtering out contigs with a relative read depth greater than " + str(args.maxreaddepth) + "... ", end="")
        sys.stdout.flush()
        contigs1 = filterContigsByReadDepth(contigs1, float(args.minreaddepth) * contigs1MedianReadDepth, float(args.maxreaddepth) * contigs1MedianReadDepth)
        contigs2 = filterContigsByReadDepth(contigs2, float(args.minreaddepth) * contigs2MedianReadDepth, float(args.maxreaddepth) * contigs2MedianReadDepth)
        print("done")
        print("   Filtered assembly 1: " + str(len(contigs1)) + " contigs, " + str(getTotalContigLength(contigs1)) + " bp")
        print("   Filtered assembly 2: " + str(len(contigs2)) + " contigs, " + str(getTotalContigLength(contigs2)) + " bp\n")

    # Save the reduced contig sets to file.
    saveContigsToFile(contigs1, tempdir + "/contigs1.fasta")
    saveContigsToFile(contigs2, tempdir + "/contigs2.fasta")

    # Build a BLAST database using the first assembly.
    print("Building BLAST database... ", end="")
    sys.stdout.flush()
    makeblastdbCommand = ["makeblastdb", "-dbtype", "nucl", "-in", tempdir + "/contigs1.fasta"]
    p = subprocess.Popen(makeblastdbCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    print("done")

    # BLAST the second assembly against the first.
    print("Running BLAST search... ", end="")
    sys.stdout.flush()
    blastnCommand = ["blastn", "-db", tempdir + "/contigs1.fasta", "-query", tempdir + "/contigs2.fasta", "-outfmt", "6 length pident sseqid sstart send sseq qseqid qstart qend qseq mismatch gaps gapopen"]
    p = subprocess.Popen(blastnCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # Save the alignments in Python objects.
    alignmentStrings = out.splitlines()
    blastAlignments = []
    for alignmentString in alignmentStrings:
        alignment = BlastAlignment(alignmentString, contigs1Dict, contigs2Dict)
        blastAlignments.append(alignment)
    print("done\n")

    # Filter the alignments for length, identity and overlap.
    print("BLAST alignments before filtering:        ", len(blastAlignments))
    blastAlignments = filterBlastAlignmentsByLength(blastAlignments, int(args.length))
    print("BLAST alignments after length filtering:  ", len(blastAlignments))
    blastAlignments = filterBlastAlignmentsByIdentity(blastAlignments, float(args.identity))
    print("BLAST alignments after identity filtering:", len(blastAlignments))
    blastAlignments = filterBlastAlignmentsByOverlap(blastAlignments, contigs1, contigs2, int(args.maxoverlap))
    print("BLAST alignments after overlap filtering: ", len(blastAlignments))

    # Display some summary information about the alignments.
    mismatches, gaps, gapopens, length = totalMismatchesGapsAndLength(blastAlignments)
    print("\nTotal alignment mismatches:", mismatches)
    print("Total alignment gap bases: ", gaps)
    print("Total alignment gap opens: ", gapopens)
    print("\nTotal alignment length:", length)
    contigs1Percent = 100.0 * length / contigs1TotalLength
    contigs2Percent = 100.0 * length / contigs2TotalLength
    print("  " + "{0:.3f}".format(contigs1Percent) + "% of assembly 1")
    print("  " + "{0:.3f}".format(contigs2Percent) + "% of assembly 2\n ")

    # Save the SNP table to file.
    if args.variants != "":
        print("Saving variants to file... ", end="")
        sys.stdout.flush()
        saveVariantsToCsvFile(blastAlignments, args.variants)
        print("done")

    # Save the alignments to file.
    if args.alignment1 != "" or args.alignment2 != "":
        print("Saving alignments to file... ", end="")
        sys.stdout.flush()
        if args.alignment1 != "":
            saveAlignmentsToFastaFile(blastAlignments, args.alignment1, True)
        if args.alignment2 != "":
            saveAlignmentsToFastaFile(blastAlignments, args.alignment2, False)
        print("done")

    # Delete the temporary files.
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    # Print a final message.
    endTime = datetime.datetime.now()
    duration = endTime - startTime
    print('\nFinished! Total time to complete:', convertTimeDeltaToReadableString(duration))






def getArguments():
    parser = argparse.ArgumentParser(description='Catpac: a Contig Alignment Tool for Pairwise Assembly Comparison')
    parser.add_argument('assembly1', help='The first set of assembled contigs')
    parser.add_argument('assembly2', help='The second set of assembled contigs')
    parser.add_argument('-1', '--alignment1', action='store', help='Save the alignments for the first assembly to this FASTA file', default="")
    parser.add_argument('-2', '--alignment2', action='store', help='Save the alignments for the second assembly to this FASTA file', default="")
    parser.add_argument('-v', '--variants', action='store', help='Save a table of variants to this CSV file', default="")
    parser.add_argument('-l', '--length', action='store', help='Minimum alignment length', default=100)
    parser.add_argument('-i', '--identity', action='store', help='Minimum alignment percent identity', default=99.0)
    parser.add_argument('-o', '--maxoverlap', action='store', help='Maximum overlap between alignments', default=0)
    parser.add_argument('-n', '--minreaddepth', action='store', help='Minimum contig read depth relative to median', default=0.0)
    parser.add_argument('-x', '--maxreaddepth', action='store', help='Maximum contig read depth relative to median', default=float("inf"))

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



# This function takes a list of Contig objects and returns another list which
# only contains the contigs within the given read depth range.
def filterContigsByReadDepth(contigs, minReadDepth, maxReadDepth):

    filteredContigs = []

    for contig in contigs:
        if contig.depth >= minReadDepth and contig.depth <= maxReadDepth:
            filteredContigs.append(contig)

    return filteredContigs



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
        outfile.write('>' + contig.fullname + '\n')
        sequence = contig.sequence
        while len(sequence) > 60:
            outfile.write(sequence[0:60] + '\n')
            sequence = sequence[60:]
        outfile.write(sequence + '\n')



# If contig1 is True, then the alignments for the first set of contigs are
# saved to file.  If False, the alignments for the second set are saved.
def saveAlignmentsToFastaFile(alignments, filename, contig1):
    outfile = open(filename, 'w')

    for alignment in alignments:

        alignmentName = ""
        if contig1:
            alignmentName += alignment.contig1.fullname
            alignmentName += "_" + str(alignment.contig1Start)
            alignmentName += "_to_" + str(alignment.contig1End)
        else:
            alignmentName += alignment.contig2.fullname
            alignmentName += "_" + str(alignment.contig2Start)
            alignmentName += "_to_" + str(alignment.contig2End)

        outfile.write('>' + alignmentName + '\n')

        sequence = ""
        if contig1:
            sequence = alignment.contig1Sequence.replace("-","")
        else:
            sequence = alignment.contig2Sequence.replace("-","")

        while len(sequence) > 60:
            outfile.write(sequence[0:60] + '\n')
            sequence = sequence[60:]
        outfile.write(sequence + '\n')



def saveVariantsToCsvFile(alignments, filename):
    outfile = open(filename, 'w')
    
    # Create the CSV header line:
    outfile.write("Type,")
    outfile.write("Contig 1 sequence,")
    outfile.write("Contig 2 sequence,")
    outfile.write("Contig 1 name,")
    outfile.write("Contig 1 position,")
    outfile.write("Contig 1 depth,")
    outfile.write("Contig 1 depth relative to median depth,")
    outfile.write("Contig 1 depth robust z-score,")
    outfile.write("Contig 2 name,")
    outfile.write("Contig 2 position,")
    outfile.write("Contig 2 depth,")
    outfile.write("Contig 2 depth relative to median depth,")
    outfile.write("Contig 2 depth robust z-score")
    outfile.write("\n")

    variants = []
    for alignment in alignments:
        variants.extend(alignment.getVariants())

    for variant in variants:
        outfile.write(variant.getCsvString())
        outfile.write("\n")



def getTotalContigLength(contigs):
    totalLength = 0
    for contig in contigs:
        totalLength += contig.length
    return totalLength



def filterBlastAlignmentsByLength(alignments, minLength):
    filteredAlignments = []

    for alignment in alignments:
        if alignment.length >= minLength:
            filteredAlignments.append(alignment)

    return filteredAlignments


def filterBlastAlignmentsByIdentity(alignments, minIdentity):
    filteredAlignments = []

    for alignment in alignments:
        if alignment.percentIdentity >= minIdentity:
            filteredAlignments.append(alignment)

    return filteredAlignments



def totalMismatchesGapsAndLength(alignments):

    mismatches = 0
    gaps = 0
    gapopens = 0
    length = 0

    for alignment in alignments:
        mismatches += alignment.mismatches
        gaps += alignment.gaps
        gapopens += alignment.gapopens
        length += alignment.length

    return (mismatches, gaps, gapopens, length)



def filterBlastAlignmentsByOverlap(alignments, contigs1, contigs2, maxOverlap):

    overlappingAlignmentPairs = []
    alignmentPairs = list(itertools.combinations(alignments, 2))
    for alignmentPair in alignmentPairs:
        if doesAlignmentPairOverlap(alignmentPair, maxOverlap):
            overlappingAlignmentPairs.append(alignmentPair)

    filteredAlignments = []
    for alignment in alignments:
        if alignmentPassesOverlapFilter(alignment, overlappingAlignmentPairs):
            filteredAlignments.append(alignment)

    return filteredAlignments



# This function looks at two alignments and determines if they overlap, either
# in contig 1 or in contig 2.
# It uses sets of positions, which probably isn't a very efficient way to do
# this, but it works well enough.
def doesAlignmentPairOverlap(alignmentPair, maxOverlap):
    alignment1 = alignmentPair[0]
    alignment2 = alignmentPair[1]

    if alignment1 == alignment2:
        return False

    if alignment1.contig1 == alignment2.contig1:

        a1c1Start = alignment1.contig1Start
        a1c1End = alignment1.contig1End
        a2c1Start = alignment2.contig1Start
        a2c1End = alignment2.contig1End

        a1c1Positions = set(range(a1c1Start, a1c1End + 1))
        a2c1Positions = set(range(a2c1Start, a2c1End + 1))
        c1OverlapLength = len(a1c1Positions & a2c1Positions)
        if c1OverlapLength > maxOverlap:
            return True

    if alignment1.contig2 == alignment2.contig2:

        a1c2Start = alignment1.contig2Start
        a1c2End = alignment1.contig2End
        a2c2Start = alignment2.contig2Start
        a2c2End = alignment2.contig2End

        # Swap the positions, if necessary (happens when a hit is on the
        # reverse complement strand)
        if (a1c2Start > a1c2End):
            a1c2Start, a1c2End = a1c2End, a1c2Start
        if (a2c2Start > a2c2End):
            a2c2Start, a2c2End = a2c2End, a2c2Start

        a1c2Positions = set(range(a1c2Start, a1c2End + 1))
        a2c2Positions = set(range(a2c2Start, a2c2End + 1))
        c2OverlapLength = len(a1c2Positions & a2c2Positions)
        if c2OverlapLength > maxOverlap:
            return True


# An alignment is said to pass the overlap filter is one of the two conditions
# is true:
#    1) it is not in any overlapping pairs
#    2) it is in overlapping pairs, but it is always the longer alignment in
#       the pair
def alignmentPassesOverlapFilter(alignment, overlappingAlignmentPairs):

    for overlappingAlignmentPair in overlappingAlignmentPairs:
        alignment1 = overlappingAlignmentPair[0]
        alignment2 = overlappingAlignmentPair[1]

        # If an alignment is in an overlapping pair where they are the same
        # length, it fails the filter.
        if alignment1.length == alignment2.length and \
               (alignment == alignment1 or alignment == alignment2):
           return False

        # If an alignment is in an overlapping pair where it is the shorter
        # one, it fails the filter.
        shorterAlignment = alignment1
        if alignment2.length < alignment1.length:
            shorterAlignment = alignment2
        if alignment == shorterAlignment:
            return False

    return True



# This function finds the median read depth by base.
def getMedianReadDepthByBaseAndMedianAbsoluteDeviation(contigs):

    readDepths = []
    for contig in contigs:
        for i in range(contig.length):
            readDepths.append(contig.depth)

    sortedReadDepths = sorted(readDepths)
    medianReadDepthByBase = getMedian(sortedReadDepths)

    absoluteDeviations = []
    for baseDepth in readDepths:
        absoluteDeviations.append(abs(baseDepth - medianReadDepthByBase))

    sortedAbsoluteDeviations = sorted(absoluteDeviations)
    medianAbsoluteDeviation = 1.4826 * getMedian(sortedAbsoluteDeviations)

    return (medianReadDepthByBase, medianAbsoluteDeviation)



def getMedian(sortedLst):

    count = len(sortedLst)
    index = (count - 1) // 2

    if (count % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1]) / 2.0



def calculateRelativeDepthAndZScore(contigs, medianReadDepth, medianAbsoluteDeviation):

    for contig in contigs:
        contig.relativeDepth = contig.depth / medianReadDepth

    for contig in contigs:
        contig.robustZScore = (contig.depth - medianReadDepth) / medianAbsoluteDeviation










# This class holds a contig: its name, sequence and length.
class Contig:

    def __init__(self, name, sequence):
        self.fullname = name
        nameParts = name.split("_")
        self.shortname = nameParts[0] + "_" + nameParts[1]
        self.number = int(nameParts[1])
        self.depth = float(nameParts[5])
        self.sequence = sequence
        self.length = len(sequence)

    def __lt__(self, other):
        return self.length < other.length

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def __eq__(self, other):
        return (self.fullname == other.fullname)














# This class holds a variant between two sequences.  It can be either a SNP or a small indel.
class Variant:

    def __init__(self, variantType, contig1, contig1Position, contig1Sequence, contig2, contig2Position, contig2Sequence):
        self.variantType = variantType

        self.contig1 = contig1
        self.contig1Position = contig1Position
        self.contig1Sequence = contig1Sequence

        self.contig2 = contig2
        self.contig2Position = contig2Position
        self.contig2Sequence = contig2Sequence

    def getCsvString(self):
        csvString = self.variantType
        csvString += ","
        csvString += self.contig1Sequence
        csvString += ","
        csvString += self.contig2Sequence
        csvString += ","
        csvString += self.contig1.shortname
        csvString += ","
        csvString += str(self.contig1Position)
        csvString += ","
        csvString += str(self.contig1.depth)
        csvString += ","
        csvString += str(self.contig1.relativeDepth)
        csvString += ","
        csvString += str(self.contig1.robustZScore)
        csvString += ","
        csvString += self.contig2.shortname
        csvString += ","
        csvString += str(self.contig2Position)
        csvString += ","
        csvString += str(self.contig2.depth)
        csvString += ","
        csvString += str(self.contig2.relativeDepth)
        csvString += ","
        csvString += str(self.contig2.robustZScore)

        return csvString















# This class holds a BLAST alignment
class BlastAlignment:

    def __init__(self, blastString, contigs1Dict, contigs2Dict):
        blastStringParts = blastString.split("\t")

        self.length = int(blastStringParts[0])
        self.percentIdentity = float(blastStringParts[1])

        contig1Name = blastStringParts[2]
        self.contig1 = contigs1Dict[contig1Name]

        self.contig1Start = int(blastStringParts[3])
        self.contig1End = int(blastStringParts[4])
        self.contig1Sequence = blastStringParts[5]

        contig2Name = blastStringParts[6]
        self.contig2 = contigs2Dict[contig2Name]

        self.contig2Start = int(blastStringParts[7])
        self.contig2End = int(blastStringParts[8])
        self.contig2Sequence = blastStringParts[9]

        self.mismatches = int(blastStringParts[10])
        self.gaps = int(blastStringParts[11])
        self.gapopens = int(blastStringParts[12])
    
    def __eq__(self, other):
        return (self.contig1 == other.contig1 and
            self.contig1Start == other.contig1Start and
            self.contig1End == other.contig1End and
            self.contig2 == other.contig2 and
            self.contig2Start == other.contig2Start and
            self.contig2End == other.contig2End)

    def __str__(self):
        return self.contig1Name + ": " + str(self.contig1Start) + " to " + str(self.contig1End) + ", " + \
               self.contig2Name + ": " + str(self.contig2Start) + " to " + str(self.contig2End)

    def __repr__(self):
        return self.contig1Name + "_" + str(self.contig1Start) + "_" + str(self.contig1End) + "_" + \
               self.contig2Name + "_" + str(self.contig2Start) + "_" + str(self.contig2End)

    # This function returns a list of all the variants within the alignment.
    def getVariants(self):

        # First loop through the alignment, pulling out all variants on a
        # base-by-base basis.
        singleNucleotideVariants = []
        for i in range(len(self.contig1Sequence)):
            base1 = self.contig1Sequence[i]
            base2 = self.contig2Sequence[i]

            if base1 != base2:

                contig1Dashes, contig2Dashes = self.countDashesUpToPosition(i)

                contig1Position = self.contig1Start + i - contig1Dashes
                contig2Position = self.contig2Start + i - contig2Dashes

                if base1 != "-" and base2 != "-":
                    variant = Variant("SNP", self.contig1, contig1Position, base1, self.contig2, contig2Position, base2)
                    singleNucleotideVariants.append(variant)
                else:
                    variant = Variant("indel", self.contig1, contig1Position, base1, self.contig2, contig2Position, base2)
                    singleNucleotideVariants.append(variant)

        # Now we want to collapse multi-base indels into single variants.
        variants = []

        indelInProgress = False
        indelContig1Sequence = ""
        indelContig2Sequence = ""
        indelContig1Position = 0
        indelContig2Position = 0

        for singleNucleotideVariant in singleNucleotideVariants:

            # If the variant is a SNP, then save the SNP variant and complete
            # the indel in progress (if there is one).
            if singleNucleotideVariant.variantType == "SNP":
                if indelInProgress:
                    variant = Variant("indel", self.contig1, indelContig1Position, indelContig1Sequence, self.contig2, indelContig2Position, indelContig2Sequence)
                    variants.append(variant)

                    indelInProgress = False
                    indelContig1Sequence = ""
                    indelContig2Sequence = ""
                    indelContig1Position = 0
                    indelContig2Position = 0

                variants.append(singleNucleotideVariant)

            # If the variant is an indel, either start a new indel in progress
            # or continue an existing indel, as appropriate.
            else:
                if not indelInProgress:
                    indelInProgress = True
                    indelContig1Sequence = singleNucleotideVariant.contig1Sequence
                    indelContig2Sequence = singleNucleotideVariant.contig2Sequence
                    indelContig1Position = singleNucleotideVariant.contig1Position
                    indelContig2Position = singleNucleotideVariant.contig2Position
                
                # If an indel is in progress...
                else:
                    # If the current indel matches the one in progress, extend
                    # the one in progress.
                    if (indelContig1Sequence[0] == '-' and singleNucleotideVariant.contig1Sequence == '-') or (indelContig2Sequence[0] == '-' and singleNucleotideVariant.contig2Sequence == '-'):
                        indelContig1Sequence += singleNucleotideVariant.contig1Sequence
                        indelContig2Sequence += singleNucleotideVariant.contig2Sequence

                    # If the current indel does not match the one in progress,
                    # finish the current one and start a new one.
                    else:
                        variant = Variant("indel", self.contig1, indelContig1Position, indelContig1Sequence, self.contig2, indelContig2Position, indelContig2Sequence)
                        variants.append(variant)

                        indelInProgress = True
                        indelContig1Sequence = singleNucleotideVariant.contig1Sequence
                        indelContig2Sequence = singleNucleotideVariant.contig2Sequence
                        indelContig1Position = singleNucleotideVariant.contig1Position
                        indelContig2Position = singleNucleotideVariant.contig2Position

        # Check to see if an indel is in progress at the end, and save it if
        # so.
        if indelInProgress:
            variant = Variant("indel", self.contig1, indelContig1Position, indelContig1Sequence, self.contig2, indelContig2Position, indelContig2Sequence)
            variants.append(variant)

        return variants

    # This function counts all occurrences of a dash up to the given index
    # in the alignment.  Counts for both the contig1 and contig2 sequences
    # are returned.  This function is used to help translate an alignment
    # position to a contig position.
    def countDashesUpToPosition(self, position):
        contig1Dashes = 0
        contig2Dashes = 0

        for i in range(position):
            contig1Base = self.contig1Sequence[i]
            contig2Base = self.contig2Sequence[i]

            if contig1Base == "-":
                contig1Dashes += 1
            if contig2Base == "-":
                contig2Dashes += 1

        return (contig1Dashes, contig2Dashes)











# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()