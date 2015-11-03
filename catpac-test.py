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
import shutil
import random



def main():

    args = getArguments()

    # Run the tests
    for i in range(args.number):
        runSingleCatpacSnpTest()



def getArguments():
    
    parser = argparse.ArgumentParser(description='Catpac tester', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--number', type=int, action='store', help='The number of random tests to conduct', default=50)

    return parser.parse_args()



def runSingleCatpacSnpTest():
    
    # Make a temporary directory for the alignment files.
    testdir = os.getcwd() + '/test'
    if not os.path.exists(testdir):
        os.makedirs(testdir)

    # Create two random sequences which share a region in common.
    seq1UniqueLength1 = random.randint(0, 200)
    seq2UniqueLength1 = random.randint(0, 200)
    sharedLength = random.randint(100, 1000)
    seq1UniqueLength2 = random.randint(0, 200)
    seq2UniqueLength2 = random.randint(0, 200)
    sharedSequence = getRandomDNASequence(sharedLength)
    sequence1 = getRandomDNASequence(seq1UniqueLength1) + sharedSequence + getRandomDNASequence(seq1UniqueLength2)
    sequence2 = getRandomDNASequence(seq2UniqueLength1) + sharedSequence + getRandomDNASequence(seq2UniqueLength2)

    # Create 5 random SNPs in the shared region of sequence 2
    startOfSnpRegion = seq2UniqueLength1 + 10
    endOfSnpRegion = seq2UniqueLength1 + sharedLength - 10
    snpLocationsSeq2 = getUniqueRandomNumbers(startOfSnpRegion, endOfSnpRegion, 5)
    for snpLocationSeq2 in snpLocationsSeq2:
        sequence2 = createSnp(sequence2, snpLocationSeq2)
    additionalBasesAtStartOfSeq1 = seq1UniqueLength1 - seq2UniqueLength1
    snpLocationsSeq1 = []
    for snpLocationSeq2 in snpLocationsSeq2:
        snpLocationsSeq1.append(snpLocationSeq2 + additionalBasesAtStartOfSeq1)

    # Save the sequences to FASTA files
    sequence1FilePath = testdir + "/seq1.fasta"
    sequence2FilePath = testdir + "/seq2.fasta"
    saveSequenceToFile("NODE_1_length_" + str(len(sequence1)) + "_cov_100.0", sequence1, sequence1FilePath)
    saveSequenceToFile("NODE_2_length_" + str(len(sequence2)) + "_cov_100.0", sequence2, sequence2FilePath)

    # Run Catpac on the two sequences and save the variants to file.
    variantsFilePath = testdir + "/variants.csv"
    catpacCommand = ["./catpac.py", sequence1FilePath, sequence2FilePath, "-l", "50", "-i", "90", "-v", variantsFilePath]
    p = subprocess.Popen(catpacCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # We expect to find the SNPs where we put them (plus 1 due to 0 vs 1 based
    # indexing).
    expectedSnpLocations = []
    for i in range(len(snpLocationsSeq1)):
        expectedSnpLocations.append((snpLocationsSeq1[i] + 1, snpLocationsSeq2[i] + 1))
    expectedSnpLocations.sort()

    # Look in the variants file for where the SNPs were actually found.
    actualSnpLocations = []
    variantsFile = open(variantsFilePath, 'r')
    for line in variantsFile:
        if line[0:3] == "SNP":
            lineParts = line.split(",")
            actualSnpLocations.append((int(lineParts[4]), int(lineParts[9])))

    # Exclude SNPs that aren't in the shared region.
    filteredActualSnpLocations = []
    for actualSnpLocation in actualSnpLocations:
        seq1Location = actualSnpLocation[0]
        seq2Location = actualSnpLocation[1]
        if seq1Location > seq1UniqueLength1 and seq1Location <= seq1UniqueLength1 + sharedLength and seq2Location > seq2UniqueLength1 and seq2Location <= seq2UniqueLength1 + sharedLength:
            filteredActualSnpLocations.append(actualSnpLocation)
    actualSnpLocations = filteredActualSnpLocations
    actualSnpLocations.sort()

    # Make sure each of the expected SNPs is in the actual SNPs.
    testPassed = True
    for expectedSnpLocation in expectedSnpLocations:
        if expectedSnpLocation not in actualSnpLocations:
            testPassed = False
            break

    # Make sure the number of found SNPs is the number of expected SNPs:
    if len(actualSnpLocations) != len(expectedSnpLocations):
        testPassed = False

    print("\nExpected SNP locations:", expectedSnpLocations)
    print("Actual SNP locations:  ", actualSnpLocations)

    if testPassed:
        print("PASS")
    else:
        print("FAIL")
        quit()

    # Delete the temporary files.
    if os.path.exists(testdir):
        shutil.rmtree(testdir)



def getRandomDNASequence(length):
    randomSequence = ""
    for i in range(length):
        randomSequence += getRandomBase()
    return randomSequence


def getRandomBase():
    randomNum = random.randint(1, 4)
    if randomNum == 1:
        return "A"
    elif randomNum == 2:
        return "C"
    elif randomNum == 3:
        return "G"
    else:
        return "T"


def getUniqueRandomNumbers(rangeStart, rangeEnd, count):
    randomNumbers = []
    for i in range(count):
        randomNumber = random.randint(rangeStart, rangeEnd)
        while randomNumber in randomNumbers:
            randomNumber = random.randint(rangeStart, rangeEnd)
        randomNumbers.append(randomNumber)
    return randomNumbers



def createSnp(sequence, location):
    oldBase = sequence[location]
    newBase = getRandomBase()
    while newBase == oldBase:
        newBase = getRandomBase()
    return sequence[0:location] + newBase + sequence[location+1:]



def saveSequenceToFile(sequenceName, sequence, filename):
    outfile = open(filename, 'w')
    outfile.write('>' + sequenceName + '\n')
    while len(sequence) > 60:
        outfile.write(sequence[0:60] + '\n')
        sequence = sequence[60:]
    outfile.write(sequence + '\n')


# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()