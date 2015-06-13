// TODO: separate all-sequences and single-sequences output
// TODO: make sure can const the returned stats, must use .at() instead of []
// TODO: error if not a Fasta sequence
// TODO: adjust default output filenames
// TODO: something is up with the assembly stats, maps need to be checked
// DONE: add CpG islands... make a transition matrix?
// DONE: get working! :-)

#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cctype>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace FastaStats {

extern "C" {
#include <zlib.h>
#include "kseq.h"
}

int         fs_debug = 0;
std::string fs_separator = ",";

// These are classes returned from queries to FastaFileStats.  They are
// filled when needed but unused otherwise.
typedef std::map<double, uint64_t>       QuantStats_t;
typedef QuantStats_t::iterator           QuantStats_t_I;
typedef QuantStats_t::const_iterator     QuantStats_t_cI;

struct SummarySequenceStats {
    std::string filename;
    uint64_t genome_size, num, total_len, min_len, num_min_len,
             short_len, num_short_len, max_len;
    double mean_len, median_len;
    uint64_t A, C, G, T, N, O, CpG;
    double GC, GC_with_N;
    uint64_t gap_num, gap_total_len, gap_min_len, gap_num_min_len,
             gap_max_len;
    double gap_mean_len, gap_median_len;
    QuantStats_t Nq, Lq, NGq, LGq;
};
struct SingleSequenceStats {
    std::string seq_name;
    std::string comment;
    uint64_t file_index;
    uint64_t length;
    uint64_t A, C, G, T, N, O, CpG;
    double GC, GC_with_N;
    uint64_t gap_num, gap_total_len, gap_min_len, gap_num_min_len,
             gap_max_len;
    double gap_mean_len, gap_median_len;
    void dump(std::ostream& os, bool header = true, std::string sep = ",") const {
        if (header)
            os << "name" << sep << "comment" << sep << "len" << sep << "idx" <<
                sep << "A" << sep << "C" << sep << "G" << sep << "T" << sep << "N" <<
                sep << "O" << sep << "GC" << sep << "GCwN" << sep << "CpG" <<
                sep << "gapnum" << sep << "gaptotlen" << sep << "gapminlen" << 
                sep << "gapnminlen" << sep << "gapmaxlen" << sep << "gapmeanlen" <<
                sep << "gapmedlen" <<
                std::endl;
        os << seq_name << sep << comment << sep << length << sep << file_index << 
            sep << A << sep << C << sep << G << sep << T << sep << N <<
            sep << O << sep << GC << sep << GC_with_N << sep << CpG <<
            sep << gap_num << sep << gap_total_len << sep << gap_min_len <<
            sep << gap_num_min_len << sep << gap_max_len << sep << gap_mean_len <<
            sep << gap_median_len <<
            std::endl;
    }
};

class FastaFileStats {
    // variables used for calculating 'assembly' stats for all sequences
    static const uint64_t     min_gap_len = 1;
    std::vector<double>       _quantiles;  // filled by set_quantiles()

  public:

    // class GapComposition --------------------------
    //
    // fill()                  catalogue gaps in sequence
    // fill_length_vector()    fill a vector with gap lengths
    // calc_stats()            calc summary stats from length vector
    // dump()                  debugging dump of contents
    // print_bed()             print BED intervals for each sequence in g
    //
    // Interval                inner class holding start and length of gap
    //         .bed_record()   create a string holding a BED record for a gap
    //
    class GapComposition {
      public:
        std::string         seq_name;
        uint64_t            num, total_len, min_len, num_min_len, max_len;
        double              mean_len, median_len;

        struct Interval {  
            std::string seq_name;
            uint64_t start, length;  // 0-based
            Interval(const std::string& s, const uint64_t st, const uint64_t len) 
                : seq_name(s), start(st), length(len)
            { }
            uint64_t start_bed() const { return start; }
            uint64_t end_bed()   const { return start + length; }
            uint64_t start_gff() const { return start + 1; }
            uint64_t end_gff()   const { return start + length; }
            const std::string bed_record() const {
                std::stringstream s;
                s << seq_name << '\t' << start_bed() << '\t' << end_bed();
                return s.str();
            }
        };

        std::vector<Interval>       g;

        GapComposition() { }

        void fill_length_vector(std::vector<uint64_t>& v) {
            v.resize(g.size());
            for (size_t i = 0; i < g.size(); ++i)
                v[i] = g[i].length;
        }
        void calc_stats() {
            num = total_len = min_len = num_min_len = max_len = 0;
            mean_len = median_len = 0.0;
            std::vector<uint64_t> lv(g.size());
            fill_length_vector(lv);
            num = uint64_t(lv.size());
            if (! num) return;
            size_t i;
            total_len = 0;
            for (i = 0; i < num; ++i)
                total_len += lv[i];
            std::sort(lv.begin(), lv.end());
            min_len = lv.front();
            for (i = 1; i < num && g[i].length == min_len; ++i);
            num_min_len = uint64_t(i);
            max_len = lv.back();
            mean_len = double(total_len) / double(num);
            if ((num % 2) == 0)
                median_len = (g[num / 2].length + g[(num / 2) - 1].length) / 2.0;
            else
                median_len = g[num / 2].length;
        }
        void fill(const std::string& n, const char* s) {
            seq_name = n;
            const char* const start = s;
            uint64_t current_gap_start = 0;
            uint64_t current_gap_len = 0;
            unsigned char c;
            while ((c = toupper(*s++))) {
                if (c == 'N') {
                    if (! current_gap_start) {
                        current_gap_start = s - start - 1;
                        current_gap_len = 1;
                    } else ++current_gap_len;
                    continue;
                } else if (current_gap_start) {
                    if (current_gap_len >= min_gap_len)
                        g.push_back(Interval(n, current_gap_start,
                                             current_gap_len));
                    current_gap_start = 0;
                    current_gap_len = 0;
                }
            }
            if (current_gap_start && current_gap_len >= min_gap_len)
                g.push_back(Interval(n, current_gap_start,
                                     current_gap_len));
            calc_stats();
        }
        void dump(std::ostream& os, bool do_all = true) const {
                os << "Gap: name=" << seq_name << " n=" << num << " tot=" <<
                    total_len << " min=" << min_len << " nmin=" <<
                    num_min_len << " max=" << max_len;
                os << " mean=" << mean_len << " med=" << median_len << std::endl;
            for (size_t i = 0; do_all && i < g.size(); ++i)
                os << "Gap: " << g[i].seq_name << " start:len " << g[i].start <<
                    ":" << g[i].length << std::endl;
        }
        void print_bed(std::ostream& os = std::cout) const {
            for (size_t i = 0; i < g.size(); ++i)
               os << g[i].bed_record() << std::endl;
        }
    };  // end class GapComposition



    // class SequenceComposition --------------------------
    //
    // .m                      std::map<char, uint64_t> of base counts
    // .CpG                    number of CpG
    // fill()                  catalogue base composition of sequence
    // dump()                  debugging dump of contents
    //
    class SequenceComposition {
        static const bool                  track_case = false;
        static const bool                  track_all_chars = false;
      public:
        typedef std::map<char, uint64_t>   comp;
        typedef comp::iterator             comp_I;
        typedef comp::const_iterator       comp_cI;

        std::string                        seq_name;
        comp                               m;
        uint64_t                           CpG;

        SequenceComposition () {
            m['A'] = 0; m['C'] = 0; m['G'] = 0; m['T'] = 0;
            m['N'] = 0; m['*'] = 0;
            CpG = 0;
        }
        double calc_GC(uint64_t sz) {
            return(double(m['G'] + m['C']) / double(sz));
        }
        void fill(const std::string& n, const char* s) {
            seq_name = n;
            unsigned char c;
            int CpG_state = 0; // 1: if C was just seen 0: otherwise
            while ((c = *s++)) {
                if (! track_case) 
                    c = toupper(c);
                switch(c) {
                    case 'N': case 'n':
                    case 'A': case 'a':
                    case 'T': case 't':
                        m[c]++;
                        CpG_state = 0;
                        break;
                    case 'C': case 'c':
                        m[c]++;
                        CpG_state = 1;
                        break;
                    case 'G': case 'g':
                        m[c]++;
                        if (CpG_state == 1) {
                            CpG += 2; // forward and reverse
                            CpG_state = 0;
                        }
                        break;
                    default:
                        m[track_all_chars ? c : '*']++;
                        CpG_state = 0;
                        break;
                }
            }
        }
        void dump(std::ostream& os = std::cout) const {
            os << "SequenceComposition: ";
            for (comp_cI it = m.begin(); it != m.end(); ++it)
                os << (it == m.begin() ? "" : "  ") << it->first <<
                    ":" << it->second;
            os << std::endl;
        }
    };  // end class SequenceComposition



    class SingleSequence {
      public:
        std::string         name;
        std::string         comment;
        uint64_t            length;
        unsigned long       file_index;
        SequenceComposition composition;
        double              GC;
        GapComposition      gaps;

        // -------- c-tor, d-tor
        //
        SingleSequence()
            : name(""), comment(""), length(0), file_index(0), GC(0.0)
        { }

        SingleSequence(const kseq_t* k, unsigned long fi = 0) {
            name.assign(k->name.s);
            if (k->comment.s)
                comment.assign(k->comment.s);
            file_index = fi;
            length = k->seq.l;
            if (! length)  // length-0 sequence, because the kseq lib can return it
                return;
            gaps.fill(name, k->seq.s);
            composition.fill(name, k->seq.s);
            GC = composition.calc_GC(length - gaps.total_len);
        }
        SingleSequenceStats sequence_stats() {
            SingleSequenceStats ans;
            ans.seq_name = name;
            ans.comment = comment;
            ans.length = length;
            ans.file_index = file_index;
            ans.A = composition.m['A'];
            ans.C = composition.m['C'];
            ans.G = composition.m['G'];
            ans.T = composition.m['T'];
            ans.N = composition.m['N'];
            ans.O = composition.m['*'];
            ans.CpG = composition.CpG;
            ans.GC = GC;
            ans.GC_with_N = composition.calc_GC(length);
            ans.gap_num = gaps.num;
            ans.gap_total_len = gaps.total_len;
            ans.gap_min_len = gaps.min_len;
            ans.gap_num_min_len = gaps.num_min_len;
            ans.gap_max_len = gaps.max_len;
            ans.gap_mean_len = gaps.mean_len;
            ans.gap_median_len = gaps.median_len;
            return ans;
        }
        void dump(std::ostream& os = std::cout) const {
            os << name << " :" << comment << ": file_index " << file_index <<
                " len " << length << std::endl;
            os << name << " ";
            composition.dump(os);
            os << name << " .gaps" << std::endl;
            gaps.dump(os);
        }
        static std::string table_header(const std::string sep = fs_separator) {
            std::stringstream s;
            s << "name";
            s << sep << "comment";
            s << sep << "len";
            s << sep << "idx";
            s << sep << "A";
            s << sep << "C";
            s << sep << "G";
            s << sep << "T";
            s << sep << "N";
            s << sep << "O";
            s << sep << "GC";
            s << sep << "GCwN";
            s << sep << "CpG";
            s << sep << "gapnum";
            s << sep << "gaptotlen";
            s << sep << "gapminlen";
            s << sep << "gapnminlen";
            s << sep << "gapmaxlen";
            s << sep << "gapmeanlen";
            s << sep << "gapmedlen";
            return s.str();
        }
        std::string table_line(const std::string sep = fs_separator) {
            std::stringstream s;
            s << name;
            s << sep << comment;
            s << sep << length;
            s << sep << file_index;
            s << sep << composition.m['A'];
            s << sep << composition.m['C'];
            s << sep << composition.m['G'];
            s << sep << composition.m['T'];
            s << sep << composition.m['N'];
            s << sep << composition.m['*'];
            s << sep << GC;
            s << sep << composition.calc_GC(length);
            s << sep << composition.CpG;
            s << sep << gaps.num;
            s << sep << gaps.total_len;
            s << sep << gaps.min_len;
            s << sep << gaps.num_min_len;
            s << sep << gaps.max_len;
            s << sep << gaps.mean_len;
            s << sep << gaps.median_len;
            return s.str();
        }
    };  // end class SingleSequence



    std::string           filename;
    uint64_t              genome_size, num, total_len,
                          min_len, num_min_len;
    static const uint64_t short_len = 500;
    uint64_t              num_short_len, max_len;
    double                mean_len, median_len;
    SequenceComposition   summary_composition;
    double                GC;
    GapComposition        summary_gaps;

    QuantStats_t     N_stats,  L_stats;   // based on total_len
    QuantStats_t     NG_stats, LG_stats;  // based on genome_size

    // one entry for each sequence
    typedef std::vector<SingleSequence>      seqs_t;
    typedef seqs_t::iterator                 seqs_I;
    typedef seqs_t::const_iterator           seqs_cI;
    seqs_t                                   seqs;
    std::vector<uint64_t>                    lengths;

    typedef std::map<std::string, size_t>    seqs_by_name_t;
    typedef seqs_by_name_t::iterator         seqs_by_name_I;
    typedef seqs_by_name_t::const_iterator   seqs_by_name_cI;
    seqs_by_name_t                           seqs_by_name;

    // -------- c-tor
    //
    FastaFileStats() { }
    FastaFileStats(const char* fn, uint64_t gsz = 0)
        : genome_size(gsz)
    {
        run(fn);
    }

    // ---- methods
    //
    // set up quantile locations
    void set_quantiles(const double from = 0.0,
                       const double to   = 1.0,
                       const double by   = 0.1) {
        _quantiles.clear();
        for (double q = from; q <= to; q += by)
            _quantiles.push_back(q);
    }
    // add a quantile and sort
    void add_quantile(const double a) {
        _quantiles.push_back(a);
        std::sort(_quantiles.begin(), _quantiles.end());
    }
    void fill_lengths() {
        lengths.resize(seqs.size());
        for (size_t i = 0; i < seqs.size(); ++i)
            lengths[i] = seqs[i].length;
    }
    // determine N50 and L50 from quantiles given a base length
    void fill_quants(QuantStats_t& N, QuantStats_t& L, uint64_t sz) {
        N.clear();
        L.clear();
        std::sort(lengths.begin(), lengths.end());
        std::reverse(lengths.begin(), lengths.end());
        size_t qi = 0;
        uint64_t cumlen = 0;
        double dcumlen;
        for (size_t si = 0; si < lengths.size() && qi < _quantiles.size(); ++si) {
            cumlen += lengths[si];
            dcumlen = double(cumlen) / double(sz);
            while (_quantiles[qi] <= dcumlen) {
                N[_quantiles[qi]] = lengths[si];
                L[_quantiles[qi]] = si + 1;
                ++qi;
            }
        }
    }
    // fill in sequence stats including quantiles etc.
    void fill() {
        fill_lengths();
        num = uint64_t(lengths.size());
        if (! num) return;
        size_t i;
        total_len = 0;
        for (i = 0; i < num; ++i)
            total_len += lengths[i];
        std::sort(lengths.begin(), lengths.end());
        min_len = lengths.front();
        for (i = 1; i < num && lengths[i] == min_len; ++i);
        num_min_len = uint64_t(i);
        for (i = 0; i < num && lengths[i] <= short_len; ++i);
        num_short_len = uint64_t(i);
        max_len = lengths.back();
        mean_len = double(total_len) / double(num);
        if ((num % 2) == 0)
            median_len = (lengths[num / 2] + lengths[(num / 2) - 1]) / 2.0;
        else
            median_len = lengths[num / 2];
        // genomic quantiles: NG20, NG50, etc.
        set_quantiles();
        fill_quants(N_stats, L_stats, total_len);
        if (genome_size > 0)
            fill_quants(NG_stats, LG_stats, genome_size);
    }
    void gaps_fill() {
        for (size_t i = 0; i < seqs.size(); ++i) {
            summary_gaps.g.insert(summary_gaps.g.end(), seqs[i].gaps.g.begin(),
                                                        seqs[i].gaps.g.end());
        }
        summary_gaps.calc_stats();
    }
    void composition_fill() {
        for (size_t i = 0; i < seqs.size(); ++i) {
            for (SequenceComposition::comp_cI it = seqs[i].composition.m.begin();
                 it != seqs[i].composition.m.end();
                 ++it) {
                summary_composition.m[it->first] += it->second;
            }
            summary_composition.CpG += seqs[i].composition.CpG;
        }
        // subtract N-gaps from total length
        GC = summary_composition.calc_GC(total_len - summary_gaps.total_len);
    }
    void dump(std::ostream& os = std::cout) const
    {
        os << "FastaFileStats: genome_size " << genome_size;
        os << "  num " << num;
        os << "  total_len " << total_len;
        os << "  min_len " << min_len;
        os << "  num_min_len " << num_min_len;
        os << "  short_len " << min_len;
        os << "  num_short_len " << num_min_len;
        os << "  max_len " << max_len;
        os << "  mean_len " << mean_len;
        os << "  median_len " << median_len;
        os << std::endl;
        os << "FastaFileStats: ";
        summary_gaps.dump(os, false);
        os << "FastaFileStats: N/L";
        for (size_t qi = 0; qi < _quantiles.size(); ++qi) {
            double q = _quantiles[qi];
            os << " " << (q * 100) << " " << N_stats.at(q) << "/" << L_stats.at(q);
        }
        os << std::endl;
        if (genome_size > 0) {
            os << "LengthStats: NG/LG";
            for (size_t qi = 0; qi < _quantiles.size(); ++qi) {
                double q = _quantiles[qi];
                os << " " << (q * 100) << " " << NG_stats.at(q) << "/" << LG_stats.at(q);
            }
            os << std::endl;
        }
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
            if (fs_debug > 0) {
                printf("name: %s\n", seq->name.s);
                if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
                printf("seq: %s\n", seq->seq.s);
                if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            }

            SingleSequence s(seq, file_index);
            seqs.push_back(s);

            size_t seqs_idx = seqs.size() - 1;

            if (seqs_by_name.find(s.name) != seqs_by_name.end()) {
                std::cerr << "duplicate sequence name: " << s.name << std::endl;
                exit(1);
            }
            seqs_by_name[s.name] = seqs_idx;  // index of back element
            if (fs_debug > 2) {
                fprintf(stderr, "inserting %s : %ld\n", s.name.c_str(), seqs_idx);
                seqs_by_name_cI it;
                if ((it = seqs_by_name.find(s.name)) != seqs_by_name.end()) {
                    std::cerr << "found just-inserted sequence (" << s.name <<
                        "): " << seqs[it->second].name << std::endl;
                    std::cerr << seqs[it->second].table_line(":") << std::endl;
                } else {
                    std::cerr << "could not find just-inserted sequence: " <<
                        s.name << std::endl;
                    exit(1);
                }
            }
        }
        kseq_destroy(seq);
        gzclose(fp);

        // done reading file and calculating sequence-specific stats
        // now calculate summary stats for the whole set of sequences
        fill();  // sequences
        gaps_fill(); // gaps
        composition_fill();  // sequence composition, gaps must come first
    }

    // -------- extract info
    //
    SummarySequenceStats summary_stats() {
        if (seqs.empty()) {
            std::cerr << "SummarySequenceStats: must call run() first" << std::endl;
            exit(1);
        }
        SummarySequenceStats ans;
        ans.filename = filename;
        ans.genome_size = genome_size;
        ans.num = num;
        ans.total_len = total_len;
        ans.min_len = min_len;
        ans.num_min_len = num_min_len;
        ans.short_len = short_len;
        ans.num_short_len = num_short_len;
        ans.max_len = max_len;
        ans.mean_len = mean_len;
        ans.median_len = median_len;
        ans.A = summary_composition.m['A'];
        ans.C = summary_composition.m['C'];
        ans.G = summary_composition.m['G'];
        ans.T = summary_composition.m['T'];
        ans.N = summary_composition.m['N'];
        ans.O = summary_composition.m['*'];
        ans.CpG = summary_composition.CpG;
        ans.GC = summary_composition.calc_GC(total_len - summary_gaps.total_len);
        ans.GC_with_N = summary_composition.calc_GC(total_len);
        ans.gap_num = summary_gaps.num;
        ans.gap_total_len = summary_gaps.total_len;
        ans.gap_min_len = summary_gaps.min_len;
        ans.gap_num_min_len = summary_gaps.num_min_len;
        ans.gap_max_len = summary_gaps.max_len;
        ans.gap_mean_len = summary_gaps.mean_len;
        ans.gap_median_len = summary_gaps.median_len;
        ans.Nq = N_stats;
        ans.Lq = L_stats;
        ans.NGq = NG_stats;
        ans.LGq = LG_stats;
        return ans;
    }

    // DONE: Answer composition query for given sequence name, including gap summary
    //
    // Statistics for a particular sequence
    SingleSequenceStats sequence_stats(const std::string& s) {
        if (fs_debug > 2)
            std::cerr << "sequence_stats: looking for sequence " << s << std::endl;
        seqs_by_name_cI it = seqs_by_name.find(s);
        if (it != seqs_by_name.end()) {
            if (fs_debug > 2) {
                std::cerr << "found sequence " << it->first << " at " <<
                    it->second << "  dumping contents... " << std::endl;
                seqs[it->second].dump(std::cerr);
            }
            return seqs[it->second].sequence_stats();
        }
        std::cerr << "could not find sequence " << s << std::endl;
        exit(1);
    }
    static std::string table_header(const std::string sep = fs_separator) {
        std::stringstream s;
        s << "filename";
        s << sep << "genomesz";
        s << sep << "num";
        s << sep << "totlen";
        s << sep << "minlen";
        s << sep << "maxlen";
        s << sep << "meanlen";
        s << sep << "medlen";
        s << sep << "A";
        s << sep << "C";
        s << sep << "G";
        s << sep << "T";
        s << sep << "N";
        s << sep << "O";
        s << sep << "GC";
        s << sep << "GCwN";
        s << sep << "CpG";
        s << sep << "gapnum";
        s << sep << "gaptotlen";
        s << sep << "gapminlen";
        s << sep << "gapmaxlen";
        s << sep << "gapmeanlen";
        s << sep << "gapmedlen";
        return s.str();
    }
    std::string table_line(const std::string sep = fs_separator) {
        std::stringstream s;
        s << filename;
        s << sep << genome_size;
        s << sep << num;
        s << sep << total_len;
        s << sep << min_len;
        s << sep << max_len;
        s << sep << mean_len;
        s << sep << median_len;
        s << sep << summary_composition.m['A'];
        s << sep << summary_composition.m['C'];
        s << sep << summary_composition.m['G'];
        s << sep << summary_composition.m['T'];
        s << sep << summary_composition.m['N'];
        s << sep << summary_composition.m['*'];
        s << sep << GC;
        s << sep << summary_composition.calc_GC(total_len);
        s << sep << summary_composition.CpG;
        s << sep << summary_gaps.num;
        s << sep << summary_gaps.total_len;
        s << sep << summary_gaps.min_len;
        s << sep << summary_gaps.max_len;
        s << sep << summary_gaps.mean_len;
        s << sep << summary_gaps.median_len;
        return s.str();
    }

    void create_summary_table(std::ostream& os = std::cout,
                              const std::string sep = fs_separator,
                              bool header = true) {
        if (header)
            os << table_header(sep) << std::endl;
        os << table_line(sep) << std::endl;
    }

    void create_sequence_table(std::ostream& os = std::cout,
                               const std::string sep = fs_separator,
                               bool header = true) {
        if (header)
            os << SingleSequence::table_header(sep) << std::endl;
        for (seqs_I it = seqs.begin(); it != seqs.end(); ++it)
            os << it->table_line(sep) << std::endl;
    }

    // -------- produce BED file describing observed N-gaps
    //
    void create_gaps_bed(const std::string& fn,
                         const std::string& query = "",
                         bool header = true) const {
        std::ofstream os(fn.c_str(), std::ofstream::out);
        if (! os.is_open()) {
            std::cerr << "could not open gaps file " << fn << std::endl;
            exit(1);
        }
        create_gaps_bed(os, query, header);
        os.close();
    }
    void create_gaps_bed(std::ostream& os = std::cout, 
                         const std::string& query = "",
                         bool header = true) const {
        if (header) {
            std::string nm = (query.empty() ? filename : query);
            os << "track name=\"gaps_" << nm << "\" ";
            os << "description=\"N-gaps (minimum length " <<
                min_gap_len << ") for " <<
                (query.empty() ? "filename " : "sequence ") <<
                nm << "\"";
            os << std::endl;
        }
        // BED entries for each gap, from each contig in order
        for (seqs_cI it = seqs.begin(); it != seqs.end(); ++it)
            if (query.empty() || query == it->name)
                it->gaps.print_bed(os);
    }
};  // end class FastaFileStats

}  // namespace FastaStats
