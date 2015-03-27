// TODO: add CpG islands


extern "C" {
#include <zlib.h>
#include "kseq.h"
}

#include <vector>
#include <map>
#include <cstdio>
#include <string>
#include <cctype>
#include <climits>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

KSEQ_INIT(gzFile, gzread)


// Class to hold sequence intervals, we use it for gaps.  Intervals have a
// sequence name, a start, a length, and an index that is the base (0 or 1)
// used to represent the first position of the sequence.  Note that this is
// *not* what is used for output... for BED and GFF formats, for example, the
// index is used to calculate the correct start position for either.

class SequenceInterval {
  public:
    std::string seq;
    int64_t start;
    int64_t length;
    int64_t index;
    SequenceInterval(const std::string& s, const int64_t st,
                     const int64_t l, const int64_t i = 1)
        : seq(s), start(st), length(l), index(i)
    {
        if (index != 0 && index != 1) {
            std::cerr << "SequenceInterval: index must be 0 or 1" << std::endl;
            exit(1);
        }
    }
    int64_t start_bed() const { return start - index; }
    int64_t end_bed()   const { return start - index + length; }
    int64_t start_gff() const { return start - index + 1; }
    int64_t end_gff()   const { return start - index + length; }
    const std::string bed_record() const
    {
        std::stringstream s;
        s << seq << '\t' << start_bed() << '\t' << end_bed();
        return s.str();
    }
};

class SizeStats {
  public:
    int64_t num;
    int64_t total_length;
    int64_t min_size;
    int64_t num_min_size;
    int64_t max_size;
    double  mean_size;
    double  median_size;
    SizeStats();
    SizeStats(std::vector<int64_t>& g)
    {
        fill(g);
    }
    void fill(std::vector<int64_t>& g)
    {
        num = int64_t(g.size());
        if (! num)
            return;
        size_t i;
        total_length = 0;
        for (i = 0; i < num; ++i)
            total_length += g[i];
        std::sort(g.begin(), g.end());
        min_size = g.front();
        for (i = 1; i < num && g[i] == min_size; ++i);
        num_min_size = int64_t(i);
        max_size = g.back();
        mean_size = double(total_length) / double(num);
        if ((num % 2) == 0)
            median_size = (g[num / 2] + g[(num / 2) - 1]) / 2.0;
        else
            median_size = g[num / 2];
    }
};

typedef std::map<char, int64_t>        char_count_map;
typedef char_count_map::iterator       char_count_map_I;
typedef char_count_map::const_iterator char_count_map_cI;

class SequenceStats {
  public:
    std::string         name;
    std::string         comment;
    int64_t             num;
    int64_t             length;
    char_count_map      composition;
    SizeStats           sequences;
    SizeStats           gaps;
};

class FastaSequenceStats {
  public:

    // stats.name     name of the Fasta sequence
    // stats.comment  sequence description
    // stats.length   length of the sequence

    SequenceStats                 stats;

    // N-gap sizes
    std::vector<int64_t>  g;

    // full N-gap descriptions
    std::vector<SequenceInterval> _gaps;

    // const because we don't want anyone to change these values yet
    static const bool           track_case = false;
    static const bool           track_all_chars = false;
    static const int64_t        min_gap_length = 1;

  public:

    // -------- c-tor, d-tor
    //
    FastaSequenceStats(const kseq_t* k)
    {
        stats.name.assign(k->name.s);
        stats.comment.assign(k->comment.s);
        stats.length = k->seq.l;
        stats.num = 1;
        if (! stats.length)  // length-0 sequence, because the kseq lib can return it
            return;
        calc_composition(k->seq.s);
    }

    // -------- c-tor accessory functions
    //
    void calc_composition(const char* s)
    {
        const char* const start = s;
        int64_t current_gap_start = 0;
        int64_t current_gap_length = 0;
        unsigned char c;
        while ((c = *s++)) {
            // position with index 1 is s - start
            if (! track_case) c = toupper(c);
            if (c == 'N' || c == 'n') {
                stats.composition[c]++;
                if (! current_gap_start) {
                    current_gap_start = s - start;
                    current_gap_length = 1;
                } else ++current_gap_length;
                continue;
            } else if (current_gap_start) {
                if (current_gap_length >= min_gap_length)
                    _gaps.push_back(SequenceInterval(stats.name,
                                                     current_gap_start,
                                                     current_gap_length));
                current_gap_start = 0;
                current_gap_length = 0;
            }
            switch(c) {
                case 'A': case 'C': case 'G': case 'T':
                case 'a': case 'c': case 'g': case 't':
                    stats.composition[c]++;
                    break;
                default:
                    stats.composition[track_all_chars ? c : '*']++;
                    break;
            }
        }
        if (s - start != stats.length) {
            std::cerr << stats.name << ": inconsistency between stated length " <<
                stats.length << " and calculated length " << s - start << std::endl;
            exit(1);
        }
        if (current_gap_start && current_gap_length >= min_gap_length)
            _gaps.push_back(SequenceInterval(stats.name,
                                             current_gap_start,
                                             current_gap_length));

        // Done reading the sequence for its composition.
        // Summarise gaps: fill size vector g and call stats.gaps.fill()
        g.reserve(_gaps.size());
        for (size_t i = 0; i < _gaps.size(); ++i)
            g[i] = _gaps[i].length;
        stats.gaps.fill(g);
    }

    // -------- extract info
    //
    void gaps_bed(std::ostream& os = std::cout) const
    {
        for (size_t i = 0; i < _gaps.size(); ++i)
           os << _gaps[i].bed_record() << std::endl;
    }

};

class FastaFileStats {
  public:
    const std::string               filename;

    // one entry for each sequence
    std::vector<FastaSequenceStats> seqs;

    // Summary stats for the whole file, filled after done with individual
    // sequences.
    SequenceStats                   stats;

    // -------- c-tor, d-tor
    //
    FastaFileStats(const char* fn)
        : filename(fn)
    {
        gzFile fp = gzopen(fn, "r");
        if (! fp) {
            std::cerr << "could not open fasta file " << filename << std::endl;
            exit(1);
        }
        kseq_t *seq = kseq_init(fp);
        int l;  // here is one place we might start retyping kseq
        while ((l = kseq_read(seq)) >= 0) {

            printf("name: %s\n", seq->name.s);
            if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            printf("seq: %s\n", seq->seq.s);
            if (seq->qual.l) printf("qual: %s\n", seq->qual.s);

            seqs.push_back(FastaSequenceStats(seq));

        }
        kseq_destroy(seq);
        gzclose(fp);

        // done reading file and calculating sequence-specific stats
        // now calculate summary stats for the whole genome
        stats.name.assign(fn);
        stats.num = int64_t(seqs.size());
        stats.length = 0;
        std::vector<int64_t>  s(seqs.size());  // for sequences
        std::vector<int64_t>  g;  // for gaps
        for (size_t i = 0; i < seqs.size(); ++i) {
            // add sequence length to s
            s[i] = seqs[i].stats.length;
            // concatenate gap sizes to g
            g.insert(g.end(), seqs[i].g.begin(), seqs[i].g.end());
            // fill stats.composition from seqs[i].stats.composition
            for (char_count_map_cI it = seqs[i].stats.composition.begin();
                 it != seqs[i].stats.composition.end();
                 ++it) {
                stats.composition[it->first] += it->second;
            }
            //
        }
        // do we want e.g., N50 when we generate stats for sequences?
        stats.sequences.fill(s);
        stats.gaps.fill(g);
    }

    // -------- extract info
    //
    // Statistics for the entire file
    SequenceStats stats_for_file() const
    {
        return stats;
    }

    // Statistics for a particular sequence
    SequenceStats stats_for_sequence(const std::string& s) const
    {
        size_t i;
        for (i = 0; i < seqs.size() && seqs[i].stats.name != s; ++i);
        if (i >= seqs.size()) {
            std::cerr << "could not find sequence named " << s << std::endl;
            exit(1);
        }
        return seqs[i].stats;
    }

    // TODO: Answer composition query for given sequence name, including gap summary

    // -------- produce BED file describing observed N-gaps
    //
    void create_gaps_bed(const std::string& fn) const
    {
        std::ofstream os(fn.c_str(), std::ofstream::out);
        if (! os.is_open()) {
            std::cerr << "could not open gaps file " << fn << std::endl;
            exit(1);
        }
        create_gaps_bed(os);
        os.close();
    }
    void create_gaps_bed(std::ostream& os = std::cout) const
    {
        // header
        os << "track name=\"gaps_" << filename << "\" ";
        os << "description=\"N-gaps (minimum length " <<
            FastaSequenceStats::min_gap_length << ") for Fasta file " <<
            filename << "\"";
        os << std::endl;
        // BED entries for each gap, from each contig in order
        for (size_t i = 0; i < seqs.size(); ++i)
            seqs[i].gaps_bed(os);
    }
};


int main(int argc, char *argv[])
{
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}

    // This reads the Fasta file in argv[1] and calculates stats for
    // all sequences.
    FastaFileStats fastastats(argv[1]);

    // Produce a BED file describing all gaps
    fastastats.create_gaps_bed("gaps.bed");

	return 0;
}
