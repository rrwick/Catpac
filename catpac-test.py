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

    # Make a temporary directory for the alignment files.
    testdir = os.getcwd() + '/test'
    if not os.path.exists(testdir):
        os.makedirs(testdir)

    # Run the tests
    for i in range(args.testcount):
        runSingleCatpacTest(testdir)

    # Delete the temporary files.
    if os.path.exists(testdir):
        shutil.rmtree(testdir)



def getArguments():
    
    parser = argparse.ArgumentParser(description='Catpac tester', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('testcount', type=int, help='The number of random tests to conduct')

    return parser.parse_args()



def runSingleCatpacTest(testdir):
    

    print(getRandomDNASequence(10))



def getRandomDNASequence(length):

    randomSequence = ""

    for i in range(length):
        randomNum = random.randint(1, 4)
        if randomNum == 1:
            randomSequence += "A"
        elif randomNum == 2:
            randomSequence += "C"
        elif randomNum == 3:
            randomSequence += "G"
        elif randomNum == 4:
            randomSequence += "T"

    return randomSequence



# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()