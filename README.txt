Variant filter.
---------------

How to run on bruce:

    module load python-gcc/2.6.4
    ./filter.py S15variants.csv S16final.bam S17final.bam S19final.bam

The first argument is a CSV file containing the variant list. The remaining
arguments are BAM files containing reads of the samples being compared. The
BAM files must be accompanied by an index (.bai) file.

In Excel you can save a spreadsheet as a CSV file.

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

If you want to run it on another computer, such as your laptop, then you need
Python >= 2.6 (but probably not >= 3.0) and the pysam library installed.
