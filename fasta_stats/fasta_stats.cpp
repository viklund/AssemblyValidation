extern "C" {
#include <zlib.h>
#include "kseq.h"
}

#include <cstdio>
#include <string>
#include <cctype>
#include <sstream>

KSEQ_INIT(gzFile, gzread)


// Class to hold sequence intervals, we use it for gaps.  Intervals have a
// sequence name, a start, a length, and an index that is the base (0 or 1)
// used to represent the first position of the sequence.  Note that this is
// *not* what is used for output... for BED and GFF formats, for example, the
// index is used to calculate the correct start position for either.

class SequenceInterval {
  public:
    const std::string& seq;
    const int64_t start;
    const int64_t length;
    const int64_t index;
    SequenceInterval(const std::string& s, const int64_t st, 
                     const int64_t l, const int64_t i = 1)
        : seq(s), start(st), length(l), index(i)
    {
        if (index != 0 && index != 1) {
            std::cerr << "SequenceInterval: index must be 0 or 1" << std::endl;
            exit(1);
        }
    }
    ~SequenceInterval();
    int64_t start_bed() { return start - index; }
    int64_t end_bed()   { return start - index + length; }
    int64_t start_gff() { return start - index + 1; }
    int64_t end_gff()   { return start - index + length; }
    const std::string& bed_entry()
    { 
        std::stringstream s;
        s << seq << '\t' << start_bed() << '\t' << end_bed() << std::endl;
        return s.str();
    }
};


class FastaSequenceStats {
  public:
    static bool           track_case = false;
    static bool           track_all_chars = false;
  public:
    const std::string     name;
    const std::string     comment;
    const int64_t        length;
    const std::string     seq;
    std::map<char, int64_t>      composition;
    std::vector<SequenceInterval> gaps;
  public:
    //FastaSequenceStats(const char* n, const char* c, int64_t l, const char* s)
    FastaSequenceStats(const kseq_t* k)
        : name(k->name.s), comment(k->comment.s), length(k->seq.l), seq(k->seq.s)
    {
        if (! length)  // length-0 sequence, because the kseq lib can return it
            return;
        if (! calc_composition(s)) {
            std::cerr << "composition could not be calculated" << std::endl;
            exit(1);
        }
    }
    ~FastaSequenceStats();
    bool calc_composition(const char* s)
    {
        const char* const start = s;
        int64_t current_gap_start = 0;
        int64_t current_gap_length = 0;
        while ((char c = *s++)) {
            // position with index 1 is s - start
            if (! track_case) c = toupper(c);
            if (c == 'N' || c == 'n') {
                composition[c]++;
                if (! current_gap_start) {
                    current_gap_start = s - start;
                    current_gap_length = 1;
                } else ++current_gap_length;
                continue;
            } else if (current_gap_start) {
                gaps.push_back(SequenceInterval(name, current_gap_start, 
                                                current_gap_length));
                current_gap_start = 0;
                current_gap_length = 0;
            }
            switch(c) {
                case 'A': case 'C': case 'G': case 'T':
                case 'a': case 'c': case 'g': case 't':
                    composition[c]++;
                    break;
                default:
                    composition[track_all_chars ? c : '*']++;
                    break;
            }
        }
        if (s - start != length) {
            std::cerr << seq << ": inconsistency between stated length " << 
                length << " and calculated length " << s - start << std::endl;
            exit(1);
        }
        if (current_gap_start)
            gaps.push_back(SequenceInterval(name, current_gap_start, 
                                            current_gap_length));
    }
};

class FastaFileStats {
  public:
    const std::string     filename;
    std::vector<FastaSequenceStats> seq_stats;
    FastaFileStats(const std::string& fn)
        : name(fn);
    {
        const char* fcstr = filename.c_str();
        gzFile fp = gzopen(fcstr, "r");
        if (! fp) {
            std::cerr << "could not open fasta file " << filename << std::endl;
            exit(1);
        }
        kseq_t *seq = kseq_init(fp);
        int l;  // here is one place we might start retyping kseq
        while ((l = kseq_read(seq)) >= 0) {
            seq_stats.push_back(FastaSequenceStats(seq));
            printf("name: %s\n", seq->name.s);
            if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            printf("seq: %s\n", seq->seq.s);
            if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        }
    }
};


int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
    std::vector<FastaFileStats> seqs;
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        SeqComp c = seq
	}
	printf("return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
