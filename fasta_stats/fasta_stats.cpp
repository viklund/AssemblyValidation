// TODO: get working! :-)
// TODO: add CpG islands... make a transition matrix?

#define DEBUG 1
static int _DEBUG = DEBUG;

extern "C" {
#include <zlib.h>
#include "kseq.h"
}

#include "SimpleOpt.h"

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


// Class to hold sequence intervals, we use it for gaps.  Intervals have a
// sequence name, a start, a length, and an index that is the base (0 or 1)
// used to represent the first position of the sequence.  Note that this is
// *not* what is used for output... for BED and GFF formats, for example, the
// index is used to calculate the correct start position for either.

class SequenceInterval {
  public:
    std::string seq;
    uint64_t start;
    uint64_t length;
  private:
    uint64_t index;
  public:
    SequenceInterval(const std::string& s, const uint64_t st,
                     const uint64_t l, const uint64_t i = 1)
        : seq(s), start(st), length(l), index(i)
    {
        if (index != 0 && index != 1) {
            std::cerr << "SequenceInterval: index must be 0 or 1" << std::endl;
            exit(1);
        }
    }
    uint64_t start_bed() const { return start - index; }
    uint64_t end_bed()   const { return start - index + length; }
    uint64_t start_gff() const { return start - index + 1; }
    uint64_t end_gff()   const { return start - index + length; }
    const std::string bed_record() const
    {
        std::stringstream s;
        s << seq << '\t' << start_bed() << '\t' << end_bed();
        return s.str();
    }
};

class SizeStats {
  public:
    uint64_t num;
    uint64_t total_length;
    uint64_t min_size;
    uint64_t num_min_size;
    uint64_t max_size;
    double   mean_size;
    double   median_size;
    SizeStats()
        : num(0), total_length(0), min_size(0), num_min_size(0),
          max_size(0), mean_size(0.0), median_size(0.0)
    { }
    SizeStats(std::vector<uint64_t>& g) { fill(g); }
    SizeStats(uint64_t l) { fill(l); }
    void fill(std::vector<uint64_t>& g)
    {
        num = uint64_t(g.size());
        if (! num)
            return;
        size_t i;
        total_length = 0;
        for (i = 0; i < num; ++i)
            total_length += g[i];
        std::sort(g.begin(), g.end());
        min_size = g.front();
        for (i = 1; i < num && g[i] == min_size; ++i);
        num_min_size = uint64_t(i);
        max_size = g.back();
        mean_size = double(total_length) / double(num);
        if ((num % 2) == 0)
            median_size = (g[num / 2] + g[(num / 2) - 1]) / 2.0;
        else
            median_size = g[num / 2];
    }
    void fill(uint64_t l)
    {
        num = num_min_size = 1;
        total_length = min_size = max_size = mean_size = median_size = l;
    }
    void dump(std::ostream& os = std::cout) const
    {
        os << "SizeStats: num " << num;
        os << "  total_length " << total_length;
        os << "  min_size " << min_size;
        os << "  num_min_size " << num_min_size;
        os << "  max_size " << max_size;
        os << "  mean_size " << mean_size;
        os << "  median_size " << median_size;
        os << std::endl;
    }
};

typedef std::map<char, uint64_t>       char_count_map;
typedef char_count_map::iterator       char_count_map_I;
typedef char_count_map::const_iterator char_count_map_cI;

static void dump_composition(const char_count_map& c,
                             std::ostream& os = std::cout)
{
    os << "composition: ";
    for (char_count_map_cI it = c.begin(); it != c.end(); ++it)
        os << (it == c.begin() ? "" : "  ") << it->first << ":" << it->second;
    os << std::endl;
}

class SequenceStats {
  public:
    std::string         name;
    std::string         comment;
    unsigned long       num;
    uint64_t            length;
    unsigned long       file_index;
    char_count_map      composition;
    SizeStats           sequences;
    SizeStats           gaps;
  public:
    void dump(std::ostream& os = std::cout) const
    {
        os << name << " :" << comment << ": file_index " << file_index <<
            " num " << num << " len " << length << std::endl;
        os << name << " ";
        dump_composition(composition, os);
        os << name << " .sequences" << std::endl;
        sequences.dump(os);
        os << name << " .gaps" << std::endl;
        gaps.dump(os);
    }
};

class FastaSequenceStats {
  public:

    // stats.name     name of the Fasta sequence
    // stats.comment  sequence description
    // stats.length   length of the sequence

    SequenceStats                 stats;

    // N-gap sizes
    std::vector<uint64_t>         g;

    // full N-gap descriptions
    std::vector<SequenceInterval> _gaps;

    // const because we don't want anyone to change these values yet
    static const bool           track_case = false;
    static const bool           track_all_chars = false;
    static const uint64_t       min_gap_length = 1;

  public:

    // -------- c-tor, d-tor
    //
    FastaSequenceStats(const kseq_t* k, unsigned long file_index = 0)
    {
        stats.name.assign(k->name.s);
        if (k->comment.s)
            stats.comment.assign(k->comment.s);
        stats.num = 1;
        stats.length = k->seq.l;
        stats.file_index = file_index;
        if (! stats.length)  // length-0 sequence, because the kseq lib can return it
            return;
        calc_composition(k->seq.s);
        stats.sequences.fill(stats.length);
    }

    // -------- c-tor accessory functions
    //
    void calc_composition(const char* s)
    {
        const char* const start = s;
        uint64_t current_gap_start = 0;
        uint64_t current_gap_length = 0;
        unsigned char c;
        while ((c = *s++)) {
            // position with index 0 is s - start - 1
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
        if (uint64_t(s - start - 1) != stats.length) {
            std::cerr << stats.name << ": inconsistency between stated length " <<
                stats.length << " and calculated length " << s - start - 1 << std::endl;
            exit(1);
        }
        if (current_gap_start && current_gap_length >= min_gap_length)
            _gaps.push_back(SequenceInterval(stats.name,
                                             current_gap_start,
                                             current_gap_length));

        // Done reading the sequence for its composition.
        // Summarise gaps: fill size vector g and call stats.gaps.fill()
        if (_DEBUG > 1) std::cout << "_gaps.size() = " << _gaps.size() << std::endl;
        g.resize(_gaps.size());
        for (size_t i = 0; i < _gaps.size(); ++i) {
            g[i] = _gaps[i].length;
            if (_DEBUG > 1) std::cout << "adding g[" << i << "] = " << g[i] << std::endl;
        }
        if (_DEBUG > 1) std::cout << "g.size() = " << g.size() << std::endl;
        stats.gaps.fill(g);
    }

    // -------- extract info
    //
    void gaps_bed(std::ostream& os = std::cout) const
    {
        for (size_t i = 0; i < _gaps.size(); ++i)
           os << _gaps[i].bed_record() << std::endl;
    }
    void dump(std::ostream& os = std::cout) const
    {
        os << "* ";
        stats.dump(os);
        dump_composition(stats.composition, os);
        os << "* " << stats.name << " gaps" << std::endl;
        gaps_bed(os);
    }

};

class FastaFileStats {
  public:
    std::string                     filename;

    // one entry for each sequence
    typedef std::vector<FastaSequenceStats>  seqs_t;
    typedef seqs_t::iterator                 seqs_I;
    typedef seqs_t::const_iterator           seqs_cI;
    seqs_t                                   seqs;

    typedef std::map<std::string, seqs_I>    seqs_by_name_t;
    typedef seqs_by_name_t::iterator         seqs_by_name_I;
    typedef seqs_by_name_t::const_iterator   seqs_by_name_cI;
    seqs_by_name_t                           seqs_by_name;

    // Summary stats for the whole file, filled after done with individual
    // sequences.
    SequenceStats                            stats;

    // -------- c-tor, d-tor
    //
    FastaFileStats() { }
    FastaFileStats(const char* fn)
    {
        run(fn);
    }

    void run(const char* fn)
    {
        if (! fn) {
            std::cerr << "must supply filename" << std::endl;
            exit(1);
        }
        filename.assign(fn);
        gzFile fp = gzopen(fn, "r");
        if (! fp) {
            std::cerr << "could not open fasta file " << filename << std::endl;
            exit(1);
        }
        kseq_t *seq = kseq_init(fp);
        int64_t l;
        uint64_t file_index = 0;
        while ((l = kseq_read(seq)) >= 0) {

            ++file_index;
            if (! seq->name.l) {
                std::cerr << "must supply sequence name" << std::endl;
                exit(1);
            }
            if (_DEBUG > 0) {
                printf("name: %s\n", seq->name.s);
                if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
                printf("seq: %s\n", seq->seq.s);
                if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            }

            FastaSequenceStats s(seq, file_index);
            seqs.push_back(s);
            if (seqs_by_name.find(s.stats.name) != seqs_by_name.end()) {
                std::cerr << "duplicate sequence name: " << s.stats.name << std::endl;
                exit(1);
            }
            seqs_by_name[s.stats.name] = seqs.end() - 1;  // it to back element

        }
        kseq_destroy(seq);
        gzclose(fp);

        // done reading file and calculating sequence-specific stats
        // now calculate summary stats for the whole genome
        stats.name.assign(fn);
        stats.num = uint64_t(seqs.size());
        stats.length = 0;
        std::vector<uint64_t>  s(seqs.size());  // for sequences
        std::vector<uint64_t>  g;  // for gaps
        for (size_t i = 0; i < seqs.size(); ++i) {
            stats.length += s[i] = seqs[i].stats.length; 
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

    // DONE: Answer composition query for given sequence name, including gap summary
    //
    // Statistics for a particular sequence
    SequenceStats stats_for_sequence(const std::string& s) const
    {
        // map interator -> second is iterator into seqs[]
        seqs_by_name_cI it = seqs_by_name.find(s);
        if (it != seqs_by_name.end())
            return it->second->stats;
        std::cerr << "could not find sequence named " << s << std::endl;
        exit(1);
    }

    void dump(std::ostream& os = std::cout) const
    {
        os << "*** FastaFileStats:  ";
        stats.dump(os);
        os << std::endl;
        for (size_t i = 0; i < seqs.size(); ++i) {
            seqs[i].dump(os);
            os << std::endl;
        }
    }
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

     //----------------- Command-line options

    enum { OPT_debug,
           OPT_help };

    CSimpleOpt::SOption options[] = {
#ifdef DEBUG
        { OPT_debug,         "--debug",           SO_REQ_SEP },
#endif
        { OPT_help,          "--help",            SO_NONE },
        { OPT_help,          "-?",                SO_NONE }, 
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            std::cerr << "invalid argument " << args.OptionText() << std::endl;
            exit(1);
        }
        if (args.OptionId() == OPT_help) {
            std::cerr << "sorry, can't help you" << std::endl;
            exit(1);
#ifdef DEBUG
        } else if (args.OptionId() == OPT_debug) {
            _DEBUG = args.OptionArg() ? atoi(args.OptionArg()) : 1;
#endif
        } else {
            std::cerr << "unknown argument " << args.OptionText() << std::endl;
            exit(1);
        }
    }
    if (args.FileCount() != 1) {
        std::cerr << "at most one sequence file as input" << std::endl;
        exit(1);
    }

    // This reads the Fasta file in argv[1] and calculates stats for
    // all sequences.
    //FastaFileStats fastastats(argv[1]);
    FastaFileStats fastastats(args.File(0));

    //fastastats.run(argv[1]);

    fastastats.dump(std::cout);

    fastastats.create_gaps_bed("gaps.bed");

	return 0;
}
