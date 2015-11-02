# Catpac: a Contig Alignment Tool for Pairwise Assembly Comparison

Catpac is a Python script that conducts a focused, iterative *de novo* assembly.

The user gives Catpac two files of assembled contigs from very closely related specimens.  Catpac will conduct a BLAST search between the two and output the alignments.

The purpose of Catpac is to produce two reduced sets of the assemblies that are more easily compared.  Small-scale variations like SNPs and small indels should be represented in these reduced sets.  However, larger scale structural variations will be filtered out.

Catpac accepts assemblies in FASTA format which follow the Velvet/SPAdes convention of headers: e.g. `>NODE_2_length_1382_cov_64.471779`

## Installation

No compilation or installation is required - just download/clone and run catpac.py.

## License

GNU General Public License, version 3
