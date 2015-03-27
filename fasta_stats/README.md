fasta_stats
===========

Working on a class for reading and reporting on contents of Fasta-format sequences.  This class will be used in FRCbam, and here we initiate development the class and use it to power a simple tool for quickly reporting contig contents.

TODO
----

* check returns of `malloc` and `realloc` in `kseq.h`
* consider retyping `kseq.h` to avoid `int`
* `free` the `kseq` buffer(s?) after constructing `FastaSequenceStats` object
* do we store sequence in `FastaSequenceStats`?  only do that while testing
* methods for `FastaSequenceStats` to return contig table entry, and `FastaFileStats` to return contig table
* methods ... ... gap table entry ... ... gap table
