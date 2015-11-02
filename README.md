# Catpac: a Contig Alignment Tool for Pairwise Assembly Comparison

Catpac is a Python script for aligning and comparing two assemblies from very closely related samples.  The user provides Catpac with two files of assembled contigs.  Catpac will conduct a BLAST search between the two assemblies, and it can output the alignments and variants.

Catpac was designed for two tasks:
* Producing a reduced set of the two assemblies which are more easily compared
* Locating small-scale variations between the two assemblies

Catpac was built with bacterial genome assemblies in mind.  It may work for larger assemblies, but this has not been tested.


## Usage

```
catpac.py [-h] [-a ALIGNMENT1] [-b ALIGNMENT2] [-v VARIANTS]
          [-l LENGTH] [-i IDENTITY] [-o MAXOVERLAP]
          [--minreldepth MINRELDEPTH] [--maxreldepth MAXRELDEPTH]
          [--mindepthz MINDEPTHZ] [--maxdepthz MAXDEPTHZ]
          assembly1 assembly2
```

##### Required arguments

Catpac must be given exactly two assembly files in FASTA format:

`catpac.py assembly1.fasta assembly2.fasta`

Catpac only accepts assemblies which follow the [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/)/[SPAdes](http://bioinf.spbau.ru/spades) convention of headers: e.g. `>NODE_2_length_1382_cov_64.471779`.  This is because Catpac needs the read depth for the contigs.  If you wish to use Catpac for assembled contigs that have a different header format, you will have to edit the initialisation method (`__init__`) of the Contig class.


##### Output alignments FASTA

The alignments for each of the assemblies can be outputted in FASTA format using the `-a` and `-b` arguments:

`catpac.py assembly1.fasta assembly2.fasta -a alignment1.fasta -b alignment2.fasta`

This usage of Catpac will produce two reduced sets of the assemblies which are easily compared.  The output files will be written in the same order and for the same strand of DNA.  For two identical assemblies, the alignments should be the same size as the original assemblies.  For less similar assemblies, the alignments will constitute the subset of the assemblies which correspond to each other.  The less similar the two assemblies are, the smaller the alignment files will be.


##### Output variant CSV file

The variants found in the alignments can be saved to a CSV file using the `-v` argument:

`catpac.py assembly1.fasta assembly2.fasta -v variants.csv`

This file will contain the variant sequences, their positions in the assembled contigs and information about the contigs' read depths.

Only variants contained within the BLAST alignments will be present here, which includes small scale varaitions like SNPs and small indels.  Larger scale structural variations will not be included.


##### Filter options

`-l` or `--length`: This specifies the minimum allowed alignment length.  BLAST alignments shorter than this will be excluded.  Default: 100 bp

`-i` or `--identity`: This specifies the minimum allowed percent identity for the alignments.  BLAST alignments with less identity will be excluded.  Default: 99.0%

`-o` or `--overlap`: To prevent assembled sequences from being duplicated in the alignment output, Catpac will exclude alignments which overlap each other by more than this amount.  Due to the nature of de Bruijn graph assemblers, it makes sense to allow a small amount of overlap equal to the k-mer size used for assembly.  Default: 0 (no overlap allowed)

`--minreldepth`, `--maxreldepth`, `--mindepthz`, `--maxdepthz`: These specify limits on the read depth of assembled contigs, relative to the assembly's median read depth.  Contigs which exceed these limits will be excluded from the analysis.  More information is below in the 'Contig read depth' section.  Default: no limits


## Contig read depth

Repeated regions of DNA can assemble into single contigs with higher than normal read depth.  Because the sequences in such contigs originate from more than one original section of DNA, it can be difficult to map reads to these contigs or generate reliable base calls.  For this reason, Catpac reports information about contig read depth and allows the user to filter based on this information (see 'Filter options' section above).

Assemblies given to Catpac have their median read depth calculated on a per base level.  This value should be close to the 'correct' read depth for a sequence which occurs only once.

##### Relative read depth

The read depth of a contig divided by the median indicates how abundant a contig is compared to the 'correct' value.  For example, a sequence which occurs twice (but assembled into a single contig) would be expected to have double the median read depth and therefore have a relative read depth of about 2.

Catpac users can filter using relative read depth.  For example, this command will only use contigs which have a relative read depth between 0.5 and 1.5 times the median read depth:

`catpac.py assembly1.fasta assembly2.fasta -v variants.csv --minreldepth 0.5 --maxreldepth 1.5`

##### Robust z-score

In addition to relative read depth, Catpac also calculates the robust z-score for each contig's read depth.  This value is similar to the more commonly used z-score, but it is less sensitive to outliers (e.g. contigs with very high read depth).  For more information, see [Birmingham, A., Selfors, L. M., Forster, T., Wrobel, D., Kennedy, C. J., Shanks, E., ... Shamu, C. E. (2009). Statistical methods for analysis of high-throughput RNA interference screens. _Nature Methods_, 6(8), 569-575.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2789971/)

The robust z-score can also be used to filter.  For example, this command will only use contigs which have a robust z-score between -2.0 and 2.0:

`catpac.py assembly1.fasta assembly2.fasta -v variants.csv --mindepthz -2.0 --mindepthz 2.0`


## Known issues

Because BLAST commands can have problems with file paths that contain spaces, Catpac can as well.  It is therefore recommended to only use Catpac in directories without spaces in their paths.

## Installation

No compilation or installation is required – just download/clone and run catpac.py from the command line.

You must have BLAST installed on the same machine as Catpac.  This can be done manually using the [instructions on the BLAST website](http://www.ncbi.nlm.nih.gov/books/NBK279671/) or via a package manager, e.g.:
* `brew install homebrew/science/blast` (using Homebrew on Mac)
* `sudo apt-get install ncbi-blast+` (using APT on Ubuntu)

It may be necessary to add executable permissions to catpac.py before using it:
`chmod u+x catpac.py`

## License

GNU General Public License, version 3
