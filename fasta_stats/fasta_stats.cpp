// TODO: add CpG islands... make a transition matrix?
// TODO: make sure can const the returned stats, must use .at() instead of []
// TODO: error if not a Fasta sequence
// TODO: adjust default output filenames
// DONE: get working! :-)

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

const std::string tab = "\t";
const std::string comma = ",";

#define DEBUG 0
static int         opt_debug = DEBUG;

static std::string opt_query;  // restrict output to only this sequence
static bool        opt_total = true;  // produce totals output
static bool        opt_sequences = true;  // produce individual sequence output
static bool        opt_header = true;  // add header to table output
static std::string opt_sep = comma;  // table column separator
static std::string opt_output;  // output file
static bool        opt_stdin = false;
static bool        from_stdin = false;
static std::string opt_gaps_bed;  // gaps BED file
static bool        opt_create_gaps_bed = true;


// Remove path from filename
std::string localBasename(const std::string& pathname, bool remove_extension = false) {
    size_t slash = pathname.find_last_of('/');
    std::string filename(pathname);
    if (slash < pathname.length())  // found '/'
        filename = pathname.substr(slash + 1);
    if (remove_extension) {
        size_t dot = filename.find_last_of('.');
        if (slash < filename.length() && slash > 0)  // found '.' not at start
            filename = filename.substr(0, dot);
    }
    return filename;
}


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

class LenStats {
  public:
    uint64_t num;
    uint64_t total_len;
    uint64_t min_len;
    uint64_t num_min_len;
    uint64_t max_len;
    double   mean_len;
    double   median_len;
    LenStats()
        : num(0), total_len(0), min_len(0), num_min_len(0),
          max_len(0), mean_len(0.0), median_len(0.0)
    { }
    LenStats(std::vector<uint64_t>& g) { fill(g); }
    LenStats(uint64_t l) { fill(l); }
    void fill(std::vector<uint64_t>& g)
    {
        num = uint64_t(g.size());
        if (! num)
            return;
        size_t i;
        total_len = 0;
        for (i = 0; i < num; ++i)
            total_len += g[i];
        std::sort(g.begin(), g.end());
        min_len = g.front();
        for (i = 1; i < num && g[i] == min_len; ++i);
        num_min_len = uint64_t(i);
        max_len = g.back();
        mean_len = double(total_len) / double(num);
        if ((num % 2) == 0)
            median_len = (g[num / 2] + g[(num / 2) - 1]) / 2.0;
        else
            median_len = g[num / 2];
    }
    void fill(uint64_t l)
    {
        num = num_min_len = 1;
        total_len = min_len = max_len = mean_len = median_len = l;
    }
    void dump(std::ostream& os = std::cout) const
    {
        os << "LenStats: num " << num;
        os << "  total_len " << total_len;
        os << "  min_len " << min_len;
        os << "  num_min_len " << num_min_len;
        os << "  max_len " << max_len;
        os << "  mean_len " << mean_len;
        os << "  median_len " << median_len;
        os << std::endl;
    }
};

typedef std::map<char, uint64_t>       char_count_map;
typedef char_count_map::iterator       char_count_map_I;
typedef char_count_map::const_iterator char_count_map_cI;
static void char_count_map_CTOR(char_count_map& m) {
    m['A'] = 0;
    m['C'] = 0;
    m['G'] = 0;
    m['T'] = 0;
    m['N'] = 0;
    m['*'] = 0;
}

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
    uint64_t            CpG;
    LenStats            sequences;
    LenStats            gaps;
  public:
    SequenceStats() {
        CpG = 0;
        char_count_map_CTOR(composition);
    }
  public:
    void dump(std::ostream& os = std::cout) const {
        os << name << " :" << comment << ": file_index " << file_index <<
            " num " << num << " len " << length << std::endl;
        os << name << " ";
        dump_composition(composition, os);
        os << name << " .sequences" << std::endl;
        sequences.dump(os);
        os << name << " .gaps" << std::endl;
        gaps.dump(os);
    }
    static std::string table_header(const std::string sep = opt_sep) {
        std::stringstream s;
        s << "name";
        s << sep << "num";
        s << sep << "len";
        s << sep << "A";
        s << sep << "C";
        s << sep << "G";
        s << sep << "T";
        s << sep << "N";
        s << sep << "O";
        s << sep << "GC";
        s << sep << "CpG";
        s << sep << "num";
        s << sep << "totlen";
        s << sep << "minlen";
        s << sep << "maxlen";
        s << sep << "meanlen";
        s << sep << "medlen";
        s << sep << "gnum";
        s << sep << "gaptotlen";
        s << sep << "gapminlen";
        s << sep << "gapmaxlen";
        s << sep << "gapmeanlen";
        s << sep << "gapmedlen";
        return s.str();
    }
    std::string table_line(const std::string sep = opt_sep) {
        std::stringstream s;
        s << name;
        s << sep << num;
        s << sep << length;
        s << sep << composition['A'];
        s << sep << composition['C'];
        s << sep << composition['G'];
        s << sep << composition['T'];
        s << sep << composition['N'];
        s << sep << composition['*'];
        s << sep << double(composition['C'] + composition['G']) / double(length);
        s << sep << CpG;
        s << sep << sequences.num;
        s << sep << sequences.total_len;
        s << sep << sequences.min_len;
        s << sep << sequences.max_len;
        s << sep << sequences.mean_len;
        s << sep << sequences.median_len;
        s << sep << gaps.num;
        s << sep << gaps.total_len;
        s << sep << gaps.min_len;
        s << sep << gaps.max_len;
        s << sep << gaps.mean_len;
        s << sep << gaps.median_len;
        return s.str();
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
    static const uint64_t       min_gap_len = 1;

  public:

    // -------- c-tor, d-tor
    //
    FastaSequenceStats(const kseq_t* k, unsigned long file_index = 0) {
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
    void calc_composition(const char* s) {
        const char* const start = s;
        uint64_t current_gap_start = 0;
        uint64_t current_gap_len = 0;
        unsigned char c;
        int CpG_state = 0; // 1: if C was just seen 0: otherwise
        while ((c = *s++)) {
            // position with index 0 is s - start - 1
            if (! track_case) c = toupper(c);
            if (c == 'N' || c == 'n') {
                stats.composition[c]++;
                if (! current_gap_start) {
                    current_gap_start = s - start;
                    current_gap_len = 1;
                } else ++current_gap_len;
                CpG_state = 0;
                continue;
            } else if (current_gap_start) {
                if (current_gap_len >= min_gap_len)
                    _gaps.push_back(SequenceInterval(stats.name,
                                                     current_gap_start,
                                                     current_gap_len));
                current_gap_start = 0;
                current_gap_len = 0;
            }
            switch(c) {
                case 'A': case 'a':
                case 'T': case 't':
                    stats.composition[c]++;
                    CpG_state = 0;
                    break;
                case 'C': case 'c':
                    stats.composition[c]++;
                    CpG_state = 1;
                    break;
                case 'G': case 'g':
                    stats.composition[c]++;
                    if (CpG_state == 1) {
                        stats.CpG += 2; // forward and reverse
                        CpG_state = 0;
                    }
                    break;
                default:
                    stats.composition[track_all_chars ? c : '*']++;
                    CpG_state = 0;
                    break;
            }
        }
        if (uint64_t(s - start - 1) != stats.length) {
            std::cerr << stats.name << ": inconsistency between stated length " <<
                stats.length << " and calculated length " << s - start - 1 << std::endl;
            exit(1);
        }
        if (current_gap_start && current_gap_len >= min_gap_len)
            _gaps.push_back(SequenceInterval(stats.name,
                                             current_gap_start,
                                             current_gap_len));

        // Done reading the sequence for its composition.
        // Summarise gaps: fill size vector g and call stats.gaps.fill()
        if (opt_debug > 2) std::cout << "_gaps.size() = " << _gaps.size() << std::endl;
        g.resize(_gaps.size());
        for (size_t i = 0; i < _gaps.size(); ++i) {
            g[i] = _gaps[i].length;
            if (opt_debug > 2) std::cout << "adding g[" << i << "] = " << g[i] << std::endl;
        }
        if (opt_debug > 2) std::cout << "g.size() = " << g.size() << std::endl;
        stats.gaps.fill(g);
    }

    // -------- extract info
    //
    void gaps_bed(std::ostream& os = std::cout) const {
        for (size_t i = 0; i < _gaps.size(); ++i)
           os << _gaps[i].bed_record() << std::endl;
    }
    void dump(std::ostream& os = std::cout) const {
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

    typedef std::map<std::string, size_t>    seqs_by_name_t;
    typedef seqs_by_name_t::iterator         seqs_by_name_I;
    typedef seqs_by_name_t::const_iterator   seqs_by_name_cI;
    seqs_by_name_t                           seqs_by_name;

    // Summary stats for the whole file, filled after done with individual
    // sequences.
    SequenceStats                            stats;

    // -------- c-tor, d-tor
    //
    FastaFileStats() { }
    FastaFileStats(const char* fn) {
        run(fn);
    }

    void run(const char* fn) {
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
            if (opt_debug > 0) {
                printf("name: %s\n", seq->name.s);
                if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
                printf("seq: %s\n", seq->seq.s);
                if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            }

            FastaSequenceStats s(seq, file_index);
            seqs.push_back(s);

            size_t seqs_idx = seqs.size() - 1;

            if (seqs_by_name.find(s.stats.name) != seqs_by_name.end()) {
                std::cerr << "duplicate sequence name: " << s.stats.name << std::endl;
                exit(1);
            }
            seqs_by_name[s.stats.name] = seqs_idx;  // index of back element
            if (opt_debug > 2) {
                fprintf(stderr, "inserting %s : %ld\n", s.stats.name.c_str(), seqs_idx);
                seqs_by_name_cI it;
                if ((it = seqs_by_name.find(s.stats.name)) != seqs_by_name.end()) {
                    std::cerr << "found just-inserted sequence (" << s.stats.name <<
                        "): " << seqs[it->second].stats.name << std::endl;
                    std::cerr << seqs[it->second].stats.table_line(":") << std::endl;
                } else {
                    std::cerr << "could not find just-inserted sequence: " <<
                        s.stats.name << std::endl;
                    exit(1);
                }
            }
        }
        kseq_destroy(seq);
        gzclose(fp);

        // done reading file and calculating sequence-specific stats
        // now calculate summary stats for the whole genome
        stats.name.assign(fn);
        stats.num = uint64_t(seqs.size());
        stats.length = 0;
        stats.CpG = 0;
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
            stats.CpG += seqs[i].stats.CpG;
        }
        // do we want e.g., N50 when we generate stats for sequences?
        stats.sequences.fill(s);
        stats.gaps.fill(g);
    }

    // -------- extract info
    //
    // Statistics for the entire file
    SequenceStats stats_for_total() const {
        return stats;
    }

    // DONE: Answer composition query for given sequence name, including gap summary
    //
    // Statistics for a particular sequence
    SequenceStats& stats_for_sequence(const std::string& s) {
        if (opt_debug > 2)
            std::cerr << "stats_for_sequence: looking for sequence " << s << std::endl;
        seqs_by_name_cI it = seqs_by_name.find(s);
        if (it != seqs_by_name.end()) {
            if (opt_debug > 2) {
                std::cerr << "found sequence " << it->first << " at " <<
                    it->second << "  dumping contents... " << std::endl;
                seqs[it->second].dump(std::cerr);
            }
            return seqs[it->second].stats;
        }
        std::cerr << "could not find sequence " << s << std::endl;
        exit(1);
    }

    void dump(std::ostream& os = std::cout) const {
        os << "*** FastaFileStats:  ";
        stats.dump(os);
        os << std::endl;
        for (seqs_cI it = seqs.begin(); it != seqs.end(); ++it) {
            std::cerr << "*it = " << &(*it) << "  " << it->stats.name << std::endl;
            it->dump(os);
            os << std::endl;
        }
    }

    void create_sequence_table(std::ostream& os = std::cout,
                               const std::string sep = opt_sep) {
        for (seqs_I it = seqs.begin(); it != seqs.end(); ++it)
            os << it->stats.table_line(sep) << std::endl;
    }

    // -------- produce BED file describing observed N-gaps
    //
    void create_gaps_bed(const std::string& fn,
                         const std::string& query = "") const {
        std::ofstream os(fn.c_str(), std::ofstream::out);
        if (! os.is_open()) {
            std::cerr << "could not open gaps file " << fn << std::endl;
            exit(1);
        }
        create_gaps_bed(os, query);
        os.close();
    }
    void create_gaps_bed(std::ostream& os = std::cout, 
                         const std::string& query = "") const {
        std::string nm = (query.empty() ? filename : query);
        // header
        os << "track name=\"gaps_" << nm << "\" ";
        os << "description=\"N-gaps (minimum length " <<
            FastaSequenceStats::min_gap_len << ") for " <<
            (query.empty() ? "filename " : "sequence ") <<
            nm << "\"";
        os << std::endl;
        // BED entries for each gap, from each contig in order
        for (seqs_cI it = seqs.begin(); it != seqs.end(); ++it)
            if (query.empty() || query == it->stats.name)
                it->gaps_bed(os);
    }
};


int main(int argc, char *argv[]) {
    enum { OPT_output,
           OPT_stdin,
           OPT_query,
           OPT_gaps_bed,   OPT_no_gaps_bed,
           OPT_tab,        OPT_comma, 
           OPT_header,     OPT_no_header, 
           OPT_total,      OPT_no_total, 
           OPT_sequences,  OPT_no_sequences, 
           OPT_debug,
           OPT_help };

    CSimpleOpt::SOption options[] = {
        { OPT_query,         "--query",           SO_REQ_SEP },
        { OPT_query,         "-q",                SO_REQ_SEP },
        { OPT_output,        "--output",          SO_REQ_SEP },
        { OPT_output,        "-o",                SO_REQ_SEP },
        { OPT_stdin,         "-",                 SO_NONE },
        { OPT_stdin,         "--stdin",           SO_NONE },
        { OPT_gaps_bed,      "--gaps-bed",        SO_REQ_SEP },
        { OPT_gaps_bed,      "-g",                SO_REQ_SEP },
        { OPT_no_gaps_bed,   "--no-gaps-bed",     SO_NONE },
        { OPT_no_gaps_bed,   "-G",                SO_NONE },
        { OPT_header,        "--header",          SO_NONE },
        { OPT_header,        "-d",                SO_NONE },
        { OPT_no_header,     "--no-header",       SO_NONE },
        { OPT_no_header,     "-D",                SO_NONE },
        { OPT_total,         "--total",           SO_NONE },
        { OPT_total,         "-t",                SO_NONE },
        { OPT_no_total,      "--no-total",        SO_NONE },
        { OPT_no_total,      "-T",                SO_NONE },
        { OPT_sequences,     "--sequences",       SO_NONE },
        { OPT_sequences,     "-s",                SO_NONE },
        { OPT_no_sequences,  "--no-sequences",    SO_NONE },
        { OPT_no_sequences,  "-S",                SO_NONE },
        { OPT_comma,         "--comma",           SO_NONE },
        { OPT_tab,           "--tab",             SO_NONE },
#ifdef DEBUG
        { OPT_debug,         "--debug",           SO_REQ_SEP },
#endif
        { OPT_help,          "--help",            SO_NONE },
        { OPT_help,          "-h",                SO_NONE }, 
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
        } else if (args.OptionId() == OPT_query) {
            opt_query = args.OptionArg();
        } else if (args.OptionId() == OPT_output) {
            opt_output = args.OptionArg();
        } else if (args.OptionId() == OPT_stdin) {
            opt_stdin = true;
        } else if (args.OptionId() == OPT_gaps_bed) {
            opt_gaps_bed = args.OptionArg();
            opt_create_gaps_bed = true;
        } else if (args.OptionId() == OPT_no_gaps_bed) {
            opt_create_gaps_bed = false;
        } else if (args.OptionId() == OPT_header) {
            opt_header = true;
        } else if (args.OptionId() == OPT_no_header) {
            opt_header = false;
        } else if (args.OptionId() == OPT_total) {
            opt_total = true;
        } else if (args.OptionId() == OPT_no_total) {
            opt_total = false;
        } else if (args.OptionId() == OPT_sequences) {
            opt_sequences = true;
        } else if (args.OptionId() == OPT_no_sequences) {
            opt_sequences = false;
        } else if (args.OptionId() == OPT_comma) {
            opt_sep = comma;
        } else if (args.OptionId() == OPT_tab) {
            opt_sep = tab;
#ifdef DEBUG
        } else if (args.OptionId() == OPT_debug) {
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : 1;
#endif
        } else {
            std::cerr << "unknown argument " << args.OptionText() << std::endl;
            exit(1);
        }
    }

    // Input and output files
    std::string input;
    if (args.FileCount() > 1) {
        std::cerr << "at most one sequence file as input" << std::endl;
        exit(1);
    } else if (args.FileCount() == 1) {
        input.assign(args.File(0));
    } else if (opt_stdin) {
        input = "/dev/stdin";
        from_stdin = true;
    } else {
        std::cerr << "at least one sequence file as input, or stdin with -/--stdin" << std::endl;
        exit(1);
    }
    if (opt_output.empty()) {
        if (opt_stdin) {
            opt_output = "/dev/stdout";
        } else {
            // assemble filename from input filename
            opt_output = localBasename(input, false) + ".stats.txt";
        }
    }
    if (opt_create_gaps_bed && opt_gaps_bed.empty())
        opt_gaps_bed = opt_stdin ? "gaps.bed" : localBasename(input, false) + ".gaps.bed";
    if (opt_debug)
        std::cerr << "input:" << input << ":   output:" << opt_output << 
            ":  gaps.bed:" << opt_gaps_bed << ":" << std::endl;
    std::ofstream output(opt_output.c_str(), std::ofstream::out);
    FastaFileStats fastastats(input.c_str());

    // produce output
    if (opt_debug > 1)
        fastastats.dump(std::cerr);
    if (opt_header)
        output << SequenceStats::table_header(opt_sep) << std::endl;
    if (! opt_query.empty()) {
        SequenceStats& s = fastastats.stats_for_sequence(opt_query);
        output << s.table_line(opt_sep) << std::endl;
        if (! opt_gaps_bed.empty())
            fastastats.create_gaps_bed(opt_gaps_bed, opt_query);
        return 0;
    }
    if (opt_total)
        output << fastastats.stats.table_line(opt_sep) << std::endl;
    if (opt_sequences)
        fastastats.create_sequence_table(output, opt_sep);
    if (opt_create_gaps_bed)
        fastastats.create_gaps_bed(opt_gaps_bed);

	return 0;
}
