--------------------------------------------------------------------------------
FAVR - Filtering and Annotation of Variants that are Rare
--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------
General description
--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------
favr_family_annotate
--------------------------------------------------------------------------------

Annotate rare variants by comparing to samples from the same family.

Command line usage:

   ./favr_family_annotate.py
      [-h | --help]
      --variants=<variant list>
      --annotations=<output CSV file with annotations added>
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

   --annotations=<output CSV file with annotations added>

      Save the annotated variants to this file. Annotations are appended to
      each variant line. Variants which are found
      in at least one family member sample are annotated with 'IN RELATIVE',
      whereas variants which are not found in at least one family member
      are annotated with 'NOT IN RELATIVE'.

   reads1.bam reads2.bam ...

      A list of bam files containing aligned sequence reads for
      particular samples. Each bam file must be accompanied by
      an index (.bai) file (but you don't mention those index
      files on the command line).

--------------------------------------------------------------------------------
favr_nonfamily_filter
--------------------------------------------------------------------------------

Filter rare variants by comparing to samples from different families.

The rule for deciding whether a variant should be kept or binned is determined
by function favr_nonfamily_classify.classify(). The default rule is:

   readCountThreshold = 1
   samplesPercent = 30
   bin = 0
   for S in samples:
      if reads_same_as_variant(S) >= readCountThreshold:
         bin++
   if (bin * 100 / totalSamples) >= samplesPercent:
      BIN
   else:
      KEEP

You can choose a different rule by modifying the classify() function.

Command line usage:

   ./favr_nonfamily_filter.py
      [-h | --help]
      --variants=<variant list as CSV file>
      --bin=<bin filename>
      --keep=<keep filename>
      --log=<log filename>
      reads1.bam reads2.bam ...

Explanation of the arguments:

   --variants=<variant list>
      same as the favr_family_annotate.py tool (described above).

   --bin=<bin filename>

      The bin file is an output of the program that contains the variants
      which were filtered out (binned).

      The format of the bin file is zero or more groups of lines of the form:

      chr1:1897857
         <vars/coverage: 1/9>
         <vars/coverage: 0/4>

      The group starts with the coordinate of the variant. It is followed
      by N lines of the form: <vars/coverage: X/Y>. Each one of these
      corresponds to one of the sample bam files,
      in the same order they were specified on the command line. X equals
      the number of reads in the sample whose base was the same as the
      variant (at the same position as the variant). Y equals the total
      number of reads in the same which covered the same position as the
      variant. So a ratio of 1/9 means: at the same position as the
      variant, 1 base in the sample was the same as the variant base, and
      9 reads in the sample covered that position.

   --keep=<keep filename>

      The keep file is an output of the program that contains the variants
      which were kept (not binned).

      The keepfile is the same format as the input variants file, and
      is sorted on chromosome coordinates.

   --log=<log filename>

      The logfile records the reasons why each variant was either binned
      or kept, based on the output of the classify function.

   reads1.bam reads2.bam ...
      same as the favr_family_annotate.py tool (described above).

--------------------------------------------------------------------------------
favr_filter_35s
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
favr_refgene
--------------------------------------------------------------------------------
