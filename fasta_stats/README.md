fasta_stats
===========

Working on a class for reading and reporting on contents of Fasta-format sequences.  This class will be used in FRCbam, and here we initiate development the class and use it to power a simple tool for quickly reporting contig contents.

Options
-------

Just a few so far.

Option  |  Meaning
------- | --------
`-q/--query NAME` | Only produce stats for sequence *NAME*
`-o/--output FILE` | Produce table output to *FILE*
`-/--stdout` | Produce table output to **stdout**
`-g/--gaps-bed FILE` | Produce BED file describing N-gaps to *FILE* (default **gaps.bed**)
`-G/--no-gaps-bed` | Do not produce BED file describing N-gaps to *FILE*
`-t/--total` | Produce a summary of total output as the first non-header output line (**default**)
`-T/--no-total` | Do not produce a total summary of the input
`-s/--sequences` | Produce a summary for each sequence (**default**)
`-S/--no-sequences` | Do not produce a summary for each sequence
`-d/--header` | Add a header to table output (**default**)
`-D/--no-header` | Do not add a header to the table output
`--comma` | Comma-separate table columns (**default**)
`--tab` | Tab-separate table columns
`--debug INT` | Debug level (0 for off)
`-h/-?/--help` | Very little help

Options processing uses Brodie Thiesfield's [SimpleOpt.h](https://github.com/brofield/simpleopt) which carries the MIT License.

TODO
----

* check returns of `malloc` and `realloc` in `kseq.h`
* `free` the `kseq` buffer(s?) after constructing `FastaSequenceStats` object

DONE
----
* Retyped `kseq.h` to use `int64_t` where possible
* Methods added to `FastaFileStats` to produce table containing total stats, and sequence stats
* Methods added to `FastaFileStats` to return `SequenceStats` object for sequence by name
* Method added to `FastaFileStats` produce gap BED
