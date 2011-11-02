FAVR - Filtering and Annotation of Variants that are Rare
---------------------------------------------------------

Version: 1.0

Authors: Tu Nguyen-Dumont (1), Fabrice Odefrey (1),
         Melissa C Southey (1), Daniel J Park (1),
         Bernard J Pope (2).

         (1) Genetic Epidemiology Laboratory, Department of Pathology.
         (2) Victorian Life Sciences Computation Initiative.

         The University of Melbourne, Australia.

Web:     https://github.com/bjpop/favr

License: ...

Requirements: Python 2.x (where x is >= 6), and the PySam library.

General description
-------------------

Characterizing genetic diversity through the analysis of massively parallel
sequence (MPS) data offers enormous potential in terms of our understanding of
predisposition to complex human disease. Great challenges remain, however,
regarding our ability to resolve those genetic variants that are genuinely
associated with disease from the millions of ‘bystanders’ and artefactual
signals. FAVR is designed to assist in the resolution of some of these issues
in the context of rare germline variants by facilitating ‘platform-steered’
artefact filtering. FAVR is a suite of tools that allow:

   (i) favr_family_annotate.py:

       Annotation of variants based on their presence/absence in
       samples from family members.

  (ii) favr_nonfamily_filter.py:

       Filtering of variants based on their presence/absence and
       abundance in samples from non-family members.

 (iii) favr_refgene.py:

       Annotation based on RefGene co-ordinates of genetic features.

  (iv) favr_filter_35s.py:

       Filtering of artefacts derived from ‘imbalanced’ paired end sequencing.

favr_family_annotate
--------------------

Annotate rare variants by comparing to samples from the same family.

Command line usage:

   ./favr_family_annotate.py
      [-h | --help]
      --variants=<variant list>
      reads1.bam reads2.bam ...

Explanation of the arguments:

   --variants=<variant list>

      List of variants, one per line. Lines that start with a
      "Residue Based Coordinate System" are considered variants,
      other lines are ignored. The Residue Based Coordinate System has
      the form:

          chromosome,coordinate,orientation,alleles

      for example:

          22,30163533,1,A/C

      that is: chromosome 22, position 30163533, orientation 1 (forward),
      reference base 'A', variant base 'C'.

      This format is supported by various genomics tools including SIFT.

      We assume the position is 1 based, which means that the first coordinate
      in a chromosome is 1 (not zero).

      The rest of the line after the coordinates can be any text, which
      will be preserved in the output.

   reads1.bam reads2.bam ...

      A list of bam files containing aligned sequence reads for
      particular samples. Each bam file must be accompanied by
      an index (.bai) file (but you don't mention those index
      files on the command line).

favr_nonfamily_filter
---------------------

Filter rare variants by comparing to samples from different families.

Command line usage:

   ./favr_nonfamily_filter
      [-h | --help]
      --variants=<variant list>
      reads1.bam reads2.bam ...

Explanation of the arguments:

   --variants=<variant list>
      same as the favr_family_annotate.py tool (described above).

   reads1.bam reads2.bam ...
      same as the favr_family_annotate.py tool (described above).

The "keep" variants will be saved in the file called "keepfile" and
the "binned" variants will be saved in the file called "binfile".

If there is already a "keepfile" or a "binfile" in the current directory,
then a number will be appended onto the end of the new filename, for
example "keepfile2" or "binfile2" to avoid overwriting existing files.

The keepfile is a CSV file containing the kept variants. It is sorted on
chromosome coordinates.

The binfile is a list of lines, where each line has the form:

   chr17:29684114 7/58 7/44 8/27

The first part of the line is the coordinate of the read. The rest of
the line contains a sequence of pairs of numbers: X/Y. The first number
indicates the number of times the variant was seen in a sample. The second
number indicates the coverage of the sample at that point.

For example "7/58" means 7 reads were the same as the variant out of a total
of 58 reads in the sample at that point.

The order of the X/Y numbers corresponds to the order of the BAM files on
the command line. Using the sample command line above, we would have:
   7/58 from S16final.bam
   7/44 from S17final.bam
   8/27 from S19final.bam

The bin file is also sorted on chromosome coordinates.

The rule for deciding whether a variant should be kept or binned is contained
in the file "bin.py". There you will find a function called "bin" which
returns True if the read is to be binned and False otherwise. You should
change this function if you want to change the rule. An example is given
to demonstrate how it works. It currently says:

   >= 5 reads of the variant must be seen in >= 20% of the samples.
