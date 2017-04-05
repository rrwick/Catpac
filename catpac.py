#!/usr/bin/env python

"""
Copyright 2017 Ryan Wick

This file is part of Catpac.

Catpac is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any
later version.

Catpac is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
Catpac.  If not, see <http:# www.gnu.org/licenses/>.
"""


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
    start_time = datetime.datetime.now()
    args = get_arguments()

    print("\nCatpac: a Contig Alignment Tool for Pairwise Assembly Comparison"
          "\n----------------------------------------------------------------\n")

    # Load in the contigs from each assembly.
    print("Loading assemblies... ", end="")
    sys.stdout.flush()
    contigs1 = load_contigs(args.assembly1)
    contigs2 = load_contigs(args.assembly2)
    contigs1_total_length = get_total_contig_length(contigs1)
    contigs2_total_length = get_total_contig_length(contigs2)
    contigs1_median_read_depth, contigs1_median_absolute_deviation = \
        get_median_read_depth_by_base_and_median_absolute_deviation(contigs1)
    contigs2_median_read_depth, contigs2_median_absolute_deviation = \
        get_median_read_depth_by_base_and_median_absolute_deviation(contigs2)
    calculate_relative_depth_and_z_score(contigs1, contigs1_median_read_depth,
                                         contigs1_median_absolute_deviation)
    calculate_relative_depth_and_z_score(contigs2, contigs2_median_read_depth,
                                         contigs2_median_absolute_deviation)
    print("done\n")
    print("Loaded assembly 1: ")
    print("   " + str(len(contigs1)) + " contigs, " + str(contigs1_total_length) + " bp")
    print("   median read depth (by base): " + str(contigs1_median_read_depth))
    print("   median absolute deviation:   " + str(contigs1_median_absolute_deviation))
    print("\nLoaded assembly 2: ")
    print("   " + str(len(contigs2)) + " contigs, " + str(contigs2_total_length) + " bp")
    print("   median read depth (by base): " + str(contigs2_median_read_depth))
    print("   median absolute deviation:   " + str(contigs2_median_absolute_deviation) + "\n")

    # Build dictionaries for the contigs, so we can later use a contig's name
    # to get the rest of the contig's details.
    contigs1_dict = {}
    for contig in contigs1:
        contigs1_dict[contig.fullname] = contig
    contigs2_dict = {}
    for contig in contigs2:
        contigs2_dict[contig.fullname] = contig

    # Make a temporary directory for the alignment files.
    tempdir = os.getcwd() + '/temp'
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # Remove contigs below the length threshold.
    if int(args.length) > 0:
        print("Filtering out contigs less than " + str(args.length) + " bp... ", end="")
        sys.stdout.flush()
        contigs1 = filter_contigs_by_length(contigs1, args.length)
        contigs2 = filter_contigs_by_length(contigs2, args.length)
        print("done")
        print("   Filtered assembly 1: " + str(len(contigs1)) + " contigs, " +
              str(get_total_contig_length(contigs1)) + " bp")
        print("   Filtered assembly 2: " + str(len(contigs2)) + " contigs, " +
              str(get_total_contig_length(contigs2)) + " bp\n")

    # Remove contigs outside the relative read depth thresholds.
    if float(args.minreldepth) > 0.0 or float(args.maxreldepth) < float("inf"):
        if float(args.minreldepth) > 0.0 and float(args.maxreldepth) < float("inf"):
            print("Filtering out contigs with a relative read depth less than " +
                  str(args.minreldepth) + " or greater than " + str(args.maxreldepth) + "... ",
                  end="")
        elif float(args.minreldepth) > 0.0:
            print("Filtering out contigs with a relative read depth less than " +
                  str(args.minreldepth) + "... ", end="")
        elif float(args.maxreldepth) < float("inf"):
            print("Filtering out contigs with a relative read depth greater than " +
                  str(args.maxreldepth) + "... ", end="")
        sys.stdout.flush()
        contigs1 = filter_contigs_by_read_depth(contigs1,
                                                args.minreldepth * contigs1_median_read_depth,
                                                args.maxreldepth * contigs1_median_read_depth)
        contigs2 = filter_contigs_by_read_depth(contigs2,
                                                args.minreldepth * contigs2_median_read_depth,
                                                args.maxreldepth * contigs2_median_read_depth)
        print("done")
        print("   Filtered assembly 1: " + str(len(contigs1)) + " contigs, " +
              str(get_total_contig_length(contigs1)) + " bp")
        print("   Filtered assembly 2: " + str(len(contigs2)) + " contigs, " +
              str(get_total_contig_length(contigs2)) + " bp\n")

    # Remove contigs outside the read depth robust z-score thresholds.
    if float(args.mindepthz) > float("-inf") or float(args.maxdepthz) < float("inf"):
        if float(args.mindepthz) > float("-inf") and float(args.maxdepthz) < float("inf"):
            print("Filtering out contigs with a read depth robust z-score less than " +
                  str(args.mindepthz) + " or greater than " + str(args.maxdepthz) + "... ", end="")
        elif float(args.mindepthz) > float("-inf"):
            print("Filtering out contigs with a read depth robust z-score less than " +
                  str(args.mindepthz) + "... ", end="")
        elif float(args.maxdepthz) < float("inf"):
            print("Filtering out contigs with a read depth robust z-score greater than " +
                  str(args.maxdepthz) + "... ", end="")
        sys.stdout.flush()
        contigs1 = filter_contigs_by_z_score(contigs1, float(args.mindepthz), float(args.maxdepthz))
        contigs2 = filter_contigs_by_z_score(contigs2, float(args.mindepthz), float(args.maxdepthz))
        print("done")
        print("   Filtered assembly 1: " + str(len(contigs1)) + " contigs, " +
              str(get_total_contig_length(contigs1)) + " bp")
        print("   Filtered assembly 2: " + str(len(contigs2)) + " contigs, " +
              str(get_total_contig_length(contigs2)) + " bp\n")

    # Save the reduced contig sets to file.
    save_contigs_to_file(contigs1, tempdir + "/contigs1.fasta")
    save_contigs_to_file(contigs2, tempdir + "/contigs2.fasta")

    # Build a BLAST database using the first assembly.
    print("Building BLAST database... ", end="")
    sys.stdout.flush()
    makeblastdb_command = ["makeblastdb", "-dbtype", "nucl", "-in", tempdir + "/contigs1.fasta"]
    p = subprocess.Popen(makeblastdb_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = p.communicate()
    print("done")

    # BLAST the second assembly against the first.
    print("Running BLAST search... ", end="")
    sys.stdout.flush()
    blastn_command = ["blastn"]
    if args.blastn:
        blastn_command += ["-task", "blastn"]
    blastn_command += ["-db", tempdir + "/contigs1.fasta", "-query", tempdir + "/contigs2.fasta",
                       "-outfmt", "6 length pident sseqid sstart send sseq qseqid "
                                  "qstart qend qseq mismatch gaps gapopen"]
    p = subprocess.Popen(blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # Save the alignments in Python objects.
    alignment_strings = out.splitlines()
    blast_alignments = []
    for alignmentString in alignment_strings:
        alignment = BlastAlignment(alignmentString, contigs1_dict, contigs2_dict)
        blast_alignments.append(alignment)
    print("done\n")

    # Filter the alignments for length, identity and overlap.
    print("BLAST alignments before filtering:        ", len(blast_alignments))
    blast_alignments = filter_blast_alignments_by_length(blast_alignments, args.length)
    print("BLAST alignments after length filtering:  ", len(blast_alignments))
    blast_alignments = filter_blast_alignments_by_identity(blast_alignments, args.identity)
    print("BLAST alignments after identity filtering:", len(blast_alignments))
    blast_alignments = filter_blast_alignments_by_overlap(blast_alignments, args.maxoverlap)
    print("BLAST alignments after overlap filtering: ", len(blast_alignments))

    # Display some summary information about the alignments.
    mismatches, gaps, gap_opens, length = total_mismatches_gaps_and_length(blast_alignments)
    print("\nTotal alignment mismatches:", mismatches)
    print("Total alignment gap bases: ", gaps)
    print("Total alignment gap opens: ", gap_opens)
    print("\nTotal alignment length:", length)
    contigs1_percent = 100.0 * length / contigs1_total_length
    contigs2_percent = 100.0 * length / contigs2_total_length
    print("  " + "{0:.3f}".format(contigs1_percent) + "% of assembly 1")
    print("  " + "{0:.3f}".format(contigs2_percent) + "% of assembly 2\n ")

    # Save the SNP table to file.
    if args.variants != "":
        print("Saving variants to file... ", end="")
        sys.stdout.flush()
        save_variants_to_csv_file(blast_alignments, args.variants)
        print("done")

    # Save the alignments to file.
    if args.alignment1 != "" or args.alignment2 != "":
        print("Saving alignments to file... ", end="")
        sys.stdout.flush()
        if args.alignment1 != "":
            save_alignments_to_fasta_file(blast_alignments, args.alignment1, True)
        if args.alignment2 != "":
            save_alignments_to_fasta_file(blast_alignments, args.alignment2, False)
        print("done")

    # Delete the temporary files.
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    # Print a final message.
    duration = datetime.datetime.now() - start_time
    print('\nFinished! Total time to complete:', convert_time_delta_to_readable_string(duration))


def get_arguments():
    """
    This little trick (adapted from stackoverflow.com/questions/9025204)
    adds a space before numbers that start with a negative sign.  It allows
    for more easily passing negative numbers via a command line argument.
    """
    for i, arg in enumerate(sys.argv):
        if len(arg) > 1 and arg[0] == '-' and arg[1].isdigit():
            sys.argv[i] = ' ' + arg

    parser = argparse.ArgumentParser(description='Catpac: a Contig Alignment Tool for Pairwise '
                                                 'Assembly Comparison',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('assembly1',
                        help='The first set of assembled contigs')
    parser.add_argument('assembly2',
                        help='The second set of assembled contigs')
    parser.add_argument('-a', '--alignment1', action='store', default="",
                        help='Save the alignments for the first assembly to this FASTA file')
    parser.add_argument('-b', '--alignment2', action='store', default="",
                        help='Save the alignments for the second assembly to this FASTA file')
    parser.add_argument('-v', '--variants', action='store', default="",
                        help='Save a table of variants to this CSV file')
    parser.add_argument('-l', '--length', action='store', type=int, default=100,
                        help='Minimum alignment length')
    parser.add_argument('-i', '--identity', action='store', type=float, default=99.0,
                        help='Minimum alignment percent identity')
    parser.add_argument('-o', '--maxoverlap', action='store', type=int, default=0,
                        help='Maximum overlap between alignments')
    parser.add_argument('--minreldepth', action='store', type=float, default=0.0,
                        help='Minimum contig read depth relative to median')
    parser.add_argument('--maxreldepth', action='store', type=float, default=float("inf"),
                        help='Maximum contig read depth relative to median')
    parser.add_argument('--mindepthz', action='store', type=float, default=float("-inf"),
                        help='Minimum contig read depth robust z-score')
    parser.add_argument('--maxdepthz', action='store', type=float, default=float("inf"),
                        help='Maximum contig read depth robust z-score')
    parser.add_argument('--blastn', action='store_true',
                        help='Use blastn (instead of megablast)')
    return parser.parse_args()


def load_contigs(contig_filename):
    """
    This function takes a contig filename and returns a list of Contig objects.
    """
    contigs = []
    with open(contig_filename, 'r') as contig_file:
        name = ''
        sequence = ''
        for line in contig_file:

            stripped_line = line.strip()

            # Skip empty lines.
            if len(stripped_line) == 0:
                continue

            # Header lines indicate the start of a new contig.
            if stripped_line[0] == '>':

                # If a contig is currently stored, save it now.
                if len(name) > 0:
                    contig = Contig(name, sequence)
                    contigs.append(contig)
                    sequence = ''

                name = stripped_line[1:]

            # If not a header line, we assume this is a sequence line.
            else:
                sequence += stripped_line

        # Save the last contig.
        if len(name) > 0:
            contig = Contig(name, sequence)
            contigs.append(contig)

    return contigs


def filter_contigs_by_length(contigs, length_threshold):
    """
    This function takes a list of Contig objects and returns another list of
    Contig objects. Only contigs with a length greater than or equal to the
    threshold will make it into the returned list.
    """
    return [x for x in contigs if x.length >= length_threshold]


def filter_contigs_by_read_depth(contigs, min_read_depth, max_read_depth):
    """
    This function takes a list of Contig objects and returns another list which
    only contains the contigs within the given read depth range.
    """
    return [x for x in contigs if min_read_depth <= x.depth <= max_read_depth]


def filter_contigs_by_z_score(contigs, min_z_score, max_z_score):
    return [x for x in contigs if min_z_score <= x.robustZScore <= max_z_score]


def convert_time_delta_to_readable_string(time_delta):
    seconds = time_delta.seconds
    hours = time_delta.days * 24
    hours += seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    seconds += time_delta.microseconds / 1000000.0
    second_string = "{:.1f}".format(seconds)
    if hours > 0:
        return str(hours) + ' h, ' + str(minutes) + ' min, ' + second_string + ' s'
    if minutes > 0:
        return str(minutes) + ' min, ' + second_string + ' s'
    return second_string + ' s'


def save_contigs_to_file(contig_list, filename):
    with open(filename, 'w') as outfile:
        for contig in contig_list:
            outfile.write('>' + contig.fullname + '\n')
            sequence = contig.sequence
            while len(sequence) > 60:
                outfile.write(sequence[0:60] + '\n')
                sequence = sequence[60:]
            outfile.write(sequence + '\n')


def save_alignments_to_fasta_file(alignments, filename, contig1):
    """
    If contig1 is True, then the alignments for the first set of contigs are
    saved to file.  If False, the alignments for the second set are saved.
    """
    outfile = open(filename, 'w')

    for alignment in alignments:

        alignment_name = ""
        if contig1:
            alignment_name += alignment.contig1.fullname
            alignment_name += "_" + str(alignment.contig1Start)
            alignment_name += "_to_" + str(alignment.contig1End)
        else:
            alignment_name += alignment.contig2.fullname
            alignment_name += "_" + str(alignment.contig2Start)
            alignment_name += "_to_" + str(alignment.contig2End)

        outfile.write('>' + alignment_name + '\n')

        if contig1:
            sequence = alignment.contig_1_sequence.replace("-", "")
        else:
            sequence = alignment.contig_2_sequence.replace("-", "")

        while len(sequence) > 60:
            outfile.write(sequence[0:60] + '\n')
            sequence = sequence[60:]
        outfile.write(sequence + '\n')


def save_variants_to_csv_file(alignments, filename):
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
        variants.extend(alignment.get_variants())

    for variant in variants:
        outfile.write(variant.get_csv_string())
        outfile.write("\n")


def get_total_contig_length(contigs):
    total_length = 0
    for contig in contigs:
        total_length += contig.length
    return total_length


def filter_blast_alignments_by_length(alignments, min_length):
    filtered_alignments = []
    for alignment in alignments:
        if alignment.length >= min_length:
            filtered_alignments.append(alignment)
    return filtered_alignments


def filter_blast_alignments_by_identity(alignments, min_identity):
    filtered_alignments = []
    for alignment in alignments:
        if alignment.percentIdentity >= min_identity:
            filtered_alignments.append(alignment)
    return filtered_alignments


def total_mismatches_gaps_and_length(alignments):
    mismatches = 0
    gaps = 0
    gap_opens = 0
    length = 0
    for alignment in alignments:
        mismatches += alignment.mismatches
        gaps += alignment.gaps
        gap_opens += alignment.gap_opens
        length += alignment.length
    return mismatches, gaps, gap_opens, length


def filter_blast_alignments_by_overlap(alignments, max_overlap):
    overlapping_alignment_pairs = []
    alignment_pairs = list(itertools.combinations(alignments, 2))
    for alignmentPair in alignment_pairs:
        if does_alignment_pair_overlap(alignmentPair, max_overlap):
            overlapping_alignment_pairs.append(alignmentPair)

    filtered_alignments = []
    for alignment in alignments:
        if alignment_passes_overlap_filter(alignment, overlapping_alignment_pairs):
            filtered_alignments.append(alignment)

    return filtered_alignments


# This function looks at two alignments and determines if they overlap, either
# in contig 1 or in contig 2.
# It uses sets of positions, which probably isn't a very efficient way to do
# this, but it works well enough.
def does_alignment_pair_overlap(alignment_pair, max_overlap):
    alignment1 = alignment_pair[0]
    alignment2 = alignment_pair[1]

    if alignment1 == alignment2:
        return False

    if alignment1.contig1 == alignment2.contig1:

        a1c1_start = alignment1.contig1Start
        a1c1_end = alignment1.contig1End
        a2c1_start = alignment2.contig1Start
        a2c1_end = alignment2.contig1End

        a1c1_positions = set(range(a1c1_start, a1c1_end + 1))
        a2c1_positions = set(range(a2c1_start, a2c1_end + 1))
        c1_overlap_length = len(a1c1_positions & a2c1_positions)
        if c1_overlap_length > max_overlap:
            return True

    if alignment1.contig2 == alignment2.contig2:

        a1c2_start = alignment1.contig2Start
        a1c2_end = alignment1.contig2End
        a2c2_start = alignment2.contig2Start
        a2c2_end = alignment2.contig2End

        # Swap the positions, if necessary (happens when a hit is on the
        # reverse complement strand)
        if a1c2_start > a1c2_end:
            a1c2_start, a1c2_end = a1c2_end, a1c2_start
        if a2c2_start > a2c2_end:
            a2c2_start, a2c2_end = a2c2_end, a2c2_start

        a1c2_positions = set(range(a1c2_start, a1c2_end + 1))
        a2c2_positions = set(range(a2c2_start, a2c2_end + 1))
        c2_overlap_length = len(a1c2_positions & a2c2_positions)
        if c2_overlap_length > max_overlap:
            return True


def alignment_passes_overlap_filter(alignment, overlapping_alignment_pairs):
    """
    An alignment is said to pass the overlap filter is one of the two conditions
    is true:
   1) it is not in any overlapping pairs
   2) it is in overlapping pairs, but it is always the longer alignment in
      the pair
    """
    for overlappingAlignmentPair in overlapping_alignment_pairs:
        alignment1 = overlappingAlignmentPair[0]
        alignment2 = overlappingAlignmentPair[1]

        # If an alignment is in an overlapping pair where they are the same
        # length, it fails the filter.
        if (alignment1.length == alignment2.length and
                (alignment == alignment1 or alignment == alignment2)):
            return False

        # If an alignment is in an overlapping pair where it is the shorter
        # one, it fails the filter.
        shorter_alignment = alignment1
        if alignment2.length < alignment1.length:
            shorter_alignment = alignment2
        if alignment == shorter_alignment:
            return False

    return True


def get_median_read_depth_by_base_and_median_absolute_deviation(contigs):
    """
    This function finds the median read depth by base.
    """
    read_depths = []
    for contig in contigs:
        for i in range(contig.length):
            read_depths.append(contig.depth)

    sorted_read_depths = sorted(read_depths)
    median_read_depth_by_base = get_median(sorted_read_depths)

    absolute_deviations = []
    for baseDepth in read_depths:
        absolute_deviations.append(abs(baseDepth - median_read_depth_by_base))

    sorted_absolute_deviations = sorted(absolute_deviations)
    median_absolute_deviation = 1.4826 * get_median(sorted_absolute_deviations)

    return median_read_depth_by_base, median_absolute_deviation


def get_median(sorted_list):
    count = len(sorted_list)
    index = (count - 1) // 2
    if count % 2:
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1]) / 2.0


def calculate_relative_depth_and_z_score(contigs, median_read_depth, median_absolute_deviation):
    for contig in contigs:
        contig.relativeDepth = contig.depth / median_read_depth
    for contig in contigs:
        if median_absolute_deviation > 0.0:
            contig.robustZScore = (contig.depth - median_read_depth) / median_absolute_deviation
        else:
            contig.robustZScore = 0.0


class Contig:
    """
    This class holds a contig: its name, sequence and length.
    """
    def __init__(self, name, sequence):
        self.fullname = name
        name_parts = name.split("_")
        self.short_name = name_parts[0] + "_" + name_parts[1]
        self.number = int(name_parts[1])
        self.depth = float(name_parts[5])
        self.sequence = sequence
        self.length = len(sequence)

    def __lt__(self, other):
        return self.length < other.length

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def __eq__(self, other):
        return self.fullname == other.fullname


class Variant:
    """
    This class holds a variant between two sequences. It can be either a SNP or a small indel.
    """

    def __init__(self, variant_type, contig1, contig_1_position, contig_1_sequence,
                 contig_1_reverse_complement, contig2, contig_2_position, contig_2_sequence,
                 contig_2_reverse_complement):
        self.variantType = variant_type

        self.contig1 = contig1
        self.contig_1_position = contig_1_position
        self.contig_1_sequence = contig_1_sequence
        self.contig_1_reverse_complement = contig_1_reverse_complement

        self.contig2 = contig2
        self.contig_2_position = contig_2_position
        self.contig_2_sequence = contig_2_sequence
        self.contig_2_reverse_complement = contig_2_reverse_complement

    def get_csv_string(self):
        csv_string = self.variantType
        csv_string += ","
        csv_string += self.contig_1_sequence
        csv_string += ","
        csv_string += self.contig_2_sequence
        csv_string += ","
        csv_string += self.contig1.short_name
        if self.contig_1_reverse_complement:
            csv_string += "_rev_comp"
        csv_string += ","
        csv_string += str(self.contig_1_position)
        csv_string += ","
        csv_string += str(self.contig1.depth)
        csv_string += ","
        csv_string += str(self.contig1.relativeDepth)
        csv_string += ","
        csv_string += str(self.contig1.robustZScore)
        csv_string += ","
        csv_string += self.contig2.short_name
        if self.contig_2_reverse_complement:
            csv_string += "_rev_comp"
        csv_string += ","
        csv_string += str(self.contig_2_position)
        csv_string += ","
        csv_string += str(self.contig2.depth)
        csv_string += ","
        csv_string += str(self.contig2.relativeDepth)
        csv_string += ","
        csv_string += str(self.contig2.robustZScore)
        return csv_string

    def can_combine(self, other):
        """
        This function checks whether the other variant can be combined with this
        one into a continuous indel.
        """
        # Only indels can be combined.
        if self.variantType != "indel" or other.variantType != "indel":
            return False

        # Only indels on the same contig can be combined.
        if self.contig_1_sequence[0] == "-" and other.contig_1_sequence[0] != "-":
            return False
        if self.contig_2_sequence[0] == "-" and other.contig_2_sequence[0] != "-":
            return False

        # Only adjacent indels can be combined.  The contig with the sequence
        # will have its positions shifted forward, but the contig with the
        # gap will not.  So we need to specifically check whether the contig
        # with the sequence has its position shifted by the right amount.
        if self.contig_1_sequence[0] != "-":
            contig1_required_position = self.contig_1_position + len(self.contig_1_sequence)
            return contig1_required_position == other.contig_1_position
        else:
            contig2_required_position = self.contig_2_position + len(self.contig_2_sequence)
            return contig2_required_position == other.contig_2_position

    def combine(self, other):
        """
        This function assumes that combination is okay, and it returns a new
        variant that is a combination of self and other.
        """
        combined_contig1_sequence = self.contig_1_sequence + other.contig_1_sequence
        combined_contig2_sequence = self.contig_2_sequence + other.contig_2_sequence
        return Variant(self.variantType,
                       self.contig1, self.contig_1_position, combined_contig1_sequence,
                       self.contig_1_reverse_complement,
                       self.contig2, self.contig_2_position, combined_contig2_sequence,
                       self.contig_2_reverse_complement)


class BlastAlignment:
    def __init__(self, blast_string, contigs_1_dict, contigs_2_dict):
        blast_string_parts = blast_string.split("\t")

        self.length = int(blast_string_parts[0])
        self.percentIdentity = float(blast_string_parts[1])

        self.contig1Name = blast_string_parts[2]
        self.contig1 = contigs_1_dict[self.contig1Name]

        self.contig1Start = int(blast_string_parts[3])
        self.contig1End = int(blast_string_parts[4])
        self.contig_1_sequence = blast_string_parts[5]
        self.contig_1_reverse_complement = self.contig1Start > self.contig1End

        self.contig2Name = blast_string_parts[6]
        self.contig2 = contigs_2_dict[self.contig2Name]

        self.contig2Start = int(blast_string_parts[7])
        self.contig2End = int(blast_string_parts[8])
        self.contig_2_sequence = blast_string_parts[9]
        self.contig_2_reverse_complement = self.contig2Start > self.contig2End

        self.mismatches = int(blast_string_parts[10])
        self.gaps = int(blast_string_parts[11])
        self.gap_opens = int(blast_string_parts[12])

    def __eq__(self, other):
        return (self.contig1 == other.contig1 and
                self.contig1Start == other.contig1Start and
                self.contig1End == other.contig1End and
                self.contig2 == other.contig2 and
                self.contig2Start == other.contig2Start and
                self.contig2End == other.contig2End)

    def __str__(self):
        return self.contig1Name + ": " + str(self.contig1Start) + " to " + str(self.contig1End) + \
               ", " + self.contig2Name + ": " + str(self.contig2Start) + " to " + \
               str(self.contig2End)

    def __repr__(self):
        return self.contig1Name + "_" + str(self.contig1Start) + "_" + str(self.contig1End) + \
               "_" + self.contig2Name + "_" + str(self.contig2Start) + "_" + str(self.contig2End)

    def get_variants(self):
        """
        This function returns a list of all the variants within the alignment.
        """
        # First loop through the alignment, pulling out all variants on a
        # base-by-base basis.
        single_nucleotide_variants = []
        for i in range(len(self.contig_1_sequence)):
            base1 = self.contig_1_sequence[i]
            base2 = self.contig_2_sequence[i]

            if base1 != base2 and base1 != 'N' and base2 != 'N':

                contig1_dashes, contig2_dashes = self.count_dashes_up_to_position(i)

                steps1 = i - contig1_dashes
                steps2 = i - contig2_dashes

                if self.contig_1_reverse_complement:
                    contig1_position = self.contig1.length - (self.contig1Start - steps1) + 1
                else:
                    contig1_position = self.contig1Start + steps1

                if self.contig_2_reverse_complement:
                    contig2_position = self.contig2.length - (self.contig2Start - steps2) + 1
                else:
                    contig2_position = self.contig2Start + steps2

                if base1 != "-" and base2 != "-":
                    variant = Variant("SNP", self.contig1, contig1_position, base1,
                                      self.contig_1_reverse_complement, self.contig2,
                                      contig2_position, base2, self.contig_2_reverse_complement)
                    single_nucleotide_variants.append(variant)
                else:
                    variant = Variant("indel", self.contig1, contig1_position, base1,
                                      self.contig_1_reverse_complement, self.contig2,
                                      contig2_position, base2, self.contig_2_reverse_complement)
                    single_nucleotide_variants.append(variant)

        # Now we want to collapse multi-base indels into single variants.
        variants = []
        single_nucleotide_variant_count = len(single_nucleotide_variants)
        i = 0
        while i < single_nucleotide_variant_count:
            variant = single_nucleotide_variants[i]

            # Combine forward as much as possible.
            while True:

                # If there is no next variant, quit this loop.
                if i + 1 >= single_nucleotide_variant_count:
                    break

                next_variant = single_nucleotide_variants[i + 1]
                if variant.can_combine(next_variant):
                    variant = variant.combine(next_variant)
                    i += 1
                else:
                    break

            variants.append(variant)
            i += 1

        return variants

    def count_dashes_up_to_position(self, position):
        """
        This function counts all occurrences of a dash up to the given index
        in the alignment.  Counts for both the contig1 and contig2 sequences
        are returned.  This function is used to help translate an alignment
        position to a contig position.
        """
        contig1_dashes = 0
        contig2_dashes = 0
        for i in range(position):
            contig1_base = self.contig_1_sequence[i]
            contig2_base = self.contig_2_sequence[i]
            if contig1_base == "-":
                contig1_dashes += 1
            if contig2_base == "-":
                contig2_dashes += 1
        return contig1_dashes, contig2_dashes


if __name__ == '__main__':
    main()
