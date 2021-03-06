// TODO: separate all-sequences and single-sequences output
// TODO: make sure can const the returned stats, must use .at() instead of []
// TODO: error if not a Fasta sequence
// TODO: adjust default output filenames
//
// TODO: incorporate CG gaps into standard gaps code
// TODO: include gaptype in calc_stats()?? or somewhere else??
// DONE: set min_gap_CG_len (now part of constructor)
// DONE: include kseq.h source in here
// DONE: something is up with the assembly stats, maps need to be checked
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

// Heng Li's kseq.h included complete.  I was thinking to remove the FastQ
// recognition for streamlining but decided against that as it is so
// lightweight already.
//
// http://lh3lh3.users.sourceforge.net/kseq.shtml

/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/* Last Modified: 12APR2009 */

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_MAX   1

typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

typedef struct __kstream_t {
    char *buf;
    int begin, end, is_eof;
    gzFile f;
} kstream_t;

typedef struct {
    kstring_t name, comment, seq, qual;
    int last_char;
    kstream_t *f;
} kseq_t;


#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

static inline kstream_t *ks_init(gzFile f)
{
    kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));
    ks->f = f;
    ks->buf = (char*)malloc(4096);
    return ks;
}
static inline void ks_destroy(kstream_t *ks)
{
    if (ks) {
        free(ks->buf);
        free(ks);
    }
}

static inline int ks_getc(kstream_t *ks)
{
    if (ks->is_eof && ks->begin >= ks->end) return -1;
    if (ks->begin >= ks->end) {
        ks->begin = 0;
        ks->end = gzread(ks->f, ks->buf, 4096);
        if (ks->end < 4096) ks->is_eof = 1;
        if (ks->end == 0) return -1;
    }
    return (int)ks->buf[ks->begin++];
}

static int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret)
{
    if (dret) *dret = 0;
    str->l = 0;
    if (ks->begin >= ks->end && ks->is_eof) return -1;
    for (;;) {
        int i;
        if (ks->begin >= ks->end) {
            if (!ks->is_eof) {
                ks->begin = 0;
                ks->end = gzread(ks->f, ks->buf, 4096);
                if (ks->end < 4096) ks->is_eof = 1;
                if (ks->end == 0) break;
            } else break;
        }
        if (delimiter > KS_SEP_MAX) {
            for (i = ks->begin; i < ks->end; ++i)
                if (ks->buf[i] == delimiter) break;
        } else if (delimiter == KS_SEP_SPACE) {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i])) break;
        } else if (delimiter == KS_SEP_TAB) {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break;
        } else i = 0; /* never come to here! */
        if (str->m - str->l < size_t(i - ks->begin + 1)) {
            str->m = str->l + (i - ks->begin) + 1;
            kroundup32(str->m);
            str->s = (char*)realloc(str->s, str->m);
        }
        memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
        str->l = str->l + (i - ks->begin);
        ks->begin = i + 1;
        if (i < ks->end) {
            if (dret) *dret = ks->buf[i];
            break;
        }
    }
    if (str->l == 0) {
        str->m = 1;
        str->s = (char*)calloc(1, 1);
    }
    str->s[str->l] = '\0';
    return str->l;
}

static inline kseq_t *kseq_init(gzFile fd)
{
    kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));
    s->f = ks_init(fd);
    return s;
}
static inline void kseq_rewind(kseq_t *ks)
{
    ks->last_char = 0;
    ks->f->is_eof = ks->f->begin = ks->f->end = 0;
}
static inline void kseq_destroy(kseq_t *ks)
{
    if (!ks) return;
    free(ks->name.s); free(ks->comment.s); free(ks->seq.s);    free(ks->qual.s);
    ks_destroy(ks->f);
    free(ks);
}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
static int64_t kseq_read(kseq_t *seq)
{
    int c;
    kstream_t *ks = seq->f;
    if (seq->last_char == 0) { /* then jump to the next header line */
        while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
        if (c == -1) return -1; /* end of file */
        seq->last_char = c;
    } /* the first header char has been read */
    seq->comment.l = seq->seq.l = seq->qual.l = 0;
    if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1;
    if (c != '\n') ks_getuntil(ks, '\n', &seq->comment, 0);
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        if (isgraph(c)) { /* printable non-space character */
            if (seq->seq.l + 1 >= seq->seq.m) { /* double the memory */
                seq->seq.m = seq->seq.l + 2;
                kroundup32(seq->seq.m); /* rounded to next closest 2^k */
                seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m);
            }
            seq->seq.s[seq->seq.l++] = (char)c;
        }
    }
    if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */
    seq->seq.s[seq->seq.l] = 0;    /* null terminated string */
    if (c != '+') return seq->seq.l; /* FASTA */
    if (seq->qual.m < seq->seq.m) {    /* allocate enough memory */
        seq->qual.m = seq->seq.m;
        seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m);
    }
    while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
    if (c == -1) return -2; /* we should not stop here */
    while ((c = ks_getc(ks)) != -1 && seq->qual.l < seq->seq.l)
        if (c >= 33 && c <= 127) seq->qual.s[seq->qual.l++] = (unsigned char)c;
    seq->qual.s[seq->qual.l] = 0; /* null terminated string */
    seq->last_char = 0;    /* we have not come to the next header line */
    if (seq->seq.l != seq->qual.l) return -2; /* qual string is shorter than seq string */
    return seq->seq.l;
}

#endif

// end of kseq.h

}  // extern "C"

///////////////////////////////////////////////////////////


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
    void dump(std::ostream& os, bool header = true, std::string sep = ",") const {
        if (header) {
            os << "filename" << sep << "genome_size" << sep << "num" <<
                sep << "total_len" <<
                sep << "min_len" << sep << "num_min_len" <<
                sep << "short_len" << sep << "num_short_len" <<
                sep << "max_len" << sep << "mean_len" << sep << "median_len" <<
                sep << "A" << sep << "C" << sep << "G" << sep << "T" <<
                sep << "N" << sep << "O" <<
                sep << "GC" << sep << "GCwN" << sep << "CpG" <<
                sep << "gapnum" << sep << "gaptotlen" << sep << "gapminlen" <<
                sep << "gapnminlen" << sep << "gapmaxlen" <<
                sep << "gapmeanlen" << sep << "gapmedlen" <<
                sep << "QuantsNotIncluded" <<
                std::endl;
        }
        os << filename << sep << genome_size << sep << num <<
            sep << total_len <<
            sep << min_len << sep << num_min_len <<
            sep << short_len << sep << num_short_len <<
            sep << max_len << sep << mean_len << sep << median_len <<
            sep << A << sep << C << sep << G << sep << T <<
            sep << N << sep << O <<
            sep << GC << sep << GC_with_N << sep << CpG <<
            sep << gap_num << sep << gap_total_len << sep << gap_min_len <<
            sep << gap_num_min_len << sep << gap_max_len <<
            sep << gap_mean_len << sep << gap_median_len <<
            sep << "QuantsNotIncluded" <<
            std::endl;
    }
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
        if (header) {
            os << "name" << sep << "comment" << sep << "len" << sep << "idx" <<
                sep << "A" << sep << "C" << sep << "G" << sep << "T" << sep << "N" <<
                sep << "O" << sep << "GC" << sep << "GCwN" << sep << "CpG" <<
                sep << "gapnum" << sep << "gaptotlen" << sep << "gapminlen" <<
                sep << "gapnminlen" << sep << "gapmaxlen" << sep << "gapmeanlen" <<
                sep << "gapmedlen" <<
                std::endl;
        }
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
    static const bool         track_gaps_CG = true;  // set in constructor
    static const uint64_t     min_gap_CG_len = 10;    // set in constructor
    std::vector<double>       _quantiles;  // filled by set_quantiles()

    static const bool         debug_gaps = false;

  public:

    enum  gaptype_t  { gap_U = 0, gap_N, gap_C, gap_G };
    static std::string gaptype_string(gaptype_t gt) {
        switch(gt) {
            case gap_N: return "gap_N"; break;
            case gap_C: return "gap_C"; break;
            case gap_G: return "gap_G"; break;
            default:    return "gap_U"; break;
        }
    }
    static std::string gaptype_string_short(gaptype_t gt) {
        switch(gt) {
            case gap_N: return "N"; break;
            case gap_C: return "C"; break;
            case gap_G: return "G"; break;
            default:    return "U"; break;
        }
    }

    // class GapComposition --------------------------
    //
    // fill_N()                catalogue N-gaps in sequence
    // fill_NCG()              catalogue N- C- and G-gaps in sequence
    // fill_length_vector()    fill a vector with gap lengths
    // calc_stats()            calc summary stats from length vector
    // dump()                  debugging dump of contents
    // print_bed()             print BED intervals for each sequence in g
    //
    // Interval                inner class holding start, length and type of gap
    //         .bed_record()   create a string holding a BED record for a gap
    //
    class GapComposition {
      public:
        std::string         seq_name;
        uint64_t            num, total_len, min_len, num_min_len, max_len;
        double              mean_len, median_len;

        struct Interval {
            std::string seq_name;
            uint64_t    start, length;  // 0-based
            gaptype_t   gaptype;
            Interval(const std::string& s, const uint64_t st, const uint64_t len,
                     const gaptype_t gt = gap_N)
                : seq_name(s), start(st), length(len), gaptype(gt)
            { 
                if (debug_gaps)
                    std::cerr << "Interval() ctor: " << s << " s:" << start << 
                        " l:" << length << " t:" << gaptype << std::endl;
            }
            uint64_t start_bed() const { return start; }
            uint64_t end_bed()   const { return start + length; }
            uint64_t start_gff() const { return start + 1; }
            uint64_t end_gff()   const { return start + length; }
            std::string get_gaptype() const {
                return gaptype_string(gaptype);
            }
            std::string get_gaptype_short() const {
                return gaptype_string_short(gaptype);
            }
            const std::string bed_record() const {
                std::stringstream s;
                s << seq_name << '\t' << start_bed() << '\t' << end_bed() << '\t' << get_gaptype();
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
        // TODO: should this be where the gaps_CG options are set, to restrict their scope??
        // If so the output routines queue off of gaptype as set in the Interval so would
        // handle output correctly regardless
        void fill_NCG(const std::string& n, const char* s) {
            seq_name = n;
            const char* const start = s;
            uint64_t          current_gap_start = 0;
            uint64_t          current_gap_len = 0;
            gaptype_t         current_gaptype = gap_U;
            unsigned char c;
            if (debug_gaps) std::cerr << "fill_NCG: n:" << n << ":  s:" << s << ":" << std::endl;

#define CHECK_GAPLENGTH(_l, _gt) (((_gt == gap_N && _l >= min_gap_len) || ((_gt == gap_C || _gt == gap_G) && _l >= min_gap_CG_len)) ? true : false)

            while ((c = toupper(*s++))) {
                if (debug_gaps) std::cerr << "fill_NCG loop: c:" << c << "  s:" << s << 
                    ":  " << gaptype_string(current_gaptype) <<
                    " @" << current_gap_start << " ." << current_gap_len << std::endl;
                if (c == 'N') {
                    switch(current_gaptype) {
                        case gap_N:  // continue gap
                            ++current_gap_len;
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " continue @" << current_gap_start << " ." << current_gap_len << std::endl;
                            break;
                        default:  // end current gap, start a gap
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " ." << current_gap_len << " end" << std::endl;
                            if (CHECK_GAPLENGTH(current_gap_len, current_gaptype)) {
                                if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " >= min len, stashing" << std::endl;
                                g.push_back(Interval(n, current_gap_start, current_gap_len, current_gaptype));
                            }
                            // >>>>>>> FALLTHROUGH
                        case gap_U:  // start a gap
                            current_gaptype = gap_N;
                            current_gap_start = s - start - 1;
                            current_gap_len = 1;
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " started @" << current_gap_start << " ." << current_gap_len << std::endl;
                            break;
                    }
                } else if (c == 'C') {
                    switch(current_gaptype) {
                        case gap_C:  // continue gap
                            ++current_gap_len;
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " continue @" << current_gap_start << " ." << current_gap_len << std::endl;
                            break;
                        default:  // end current gap, start a gap
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " ." << current_gap_len << " end" << std::endl;
                            if (CHECK_GAPLENGTH(current_gap_len, current_gaptype)) {
                                if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " >= min len, stashing" << std::endl;
                                g.push_back(Interval(n, current_gap_start, current_gap_len, current_gaptype));
                            }
                            // >>>>>>> FALLTHROUGH
                        case gap_U:  // start a gap
                            current_gaptype = gap_C;
                            current_gap_start = s - start - 1;
                            current_gap_len = 1;
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " started @" << current_gap_start << " ." << current_gap_len << std::endl;
                            break;
                    }
                } else if (c == 'G') {
                    switch(current_gaptype) {
                        case gap_G:  // continue gap
                            ++current_gap_len;
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " continue @" << current_gap_start << " ." << current_gap_len << std::endl;
                            break;
                        default:  // end current gap, start a gap
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " ." << current_gap_len << " end" << std::endl;
                            if (CHECK_GAPLENGTH(current_gap_len, current_gaptype)) {
                                if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " >= min len, stashing" << std::endl;
                                g.push_back(Interval(n, current_gap_start, current_gap_len, current_gaptype));
                            }
                            // >>>>>>> FALLTHROUGH
                        case gap_U:  // start a gap
                            current_gaptype = gap_G;
                            current_gap_start = s - start - 1;
                            current_gap_len = 1;
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " started @" << current_gap_start << " ." << current_gap_len << std::endl;
                            break;
                    }
                } else {  // A or T or something else
                    if (current_gaptype != gap_U) {
                        if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " ." << current_gap_len << " end" << std::endl;
                        if (CHECK_GAPLENGTH(current_gap_len, current_gaptype)) {
                            if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " end @" << current_gap_start << " >= min len, stashing" << std::endl;
                            g.push_back(Interval(n, current_gap_start, current_gap_len, current_gaptype));
                        }
                        current_gap_start = 0;
                        current_gap_len = 0;
                        current_gaptype = gap_U;
                        if (debug_gaps) std::cerr << "gap reset" << std::endl;
                    }
                }
            }
            if (current_gaptype != gap_U) {
                if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " FINAL @" << current_gap_start << " ." << current_gap_len << " end" << std::endl;
                if (CHECK_GAPLENGTH(current_gap_len, current_gaptype)) {
                    if (debug_gaps) std::cerr << gaptype_string(current_gaptype) << " FINAL @" << current_gap_start << " >= min len, stashing" << std::endl;
                    g.push_back(Interval(n, current_gap_start, current_gap_len, current_gaptype));
                }
            }
            calc_stats();
        }
        void fill_N(const std::string& n, const char* s) {
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
                        g.push_back(Interval(n, current_gap_start, current_gap_len, gap_N));
                    current_gap_start = 0;
                    current_gap_len = 0;
                }
            }
            if (current_gap_start && current_gap_len >= min_gap_len)
                g.push_back(Interval(n, current_gap_start, current_gap_len, gap_N));
            calc_stats();
        }
        void dump(std::ostream& os, bool do_all = true) const {
                os << "Gap: name=" << seq_name << " n=" << num << " tot=" <<
                    total_len << " min=" << min_len << " nmin=" <<
                    num_min_len << " max=" << max_len;
                os << " mean=" << mean_len << " med=" << median_len << std::endl;
            for (size_t i = 0; do_all && i < g.size(); ++i)
                os << "Gap: " << g[i].seq_name << " start:len type " << g[i].start <<
                    ":" << g[i].length << " " << g[i].get_gaptype_short() << std::endl;
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
            gaps.fill_NCG(name, k->seq.s); //TODO
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
            s << sep << "gapnum"; //TODO add gaptype stuff later
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
            s << sep << gaps.num; //TODO add gaptype stuff later
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
    FastaFileStats(const char* fn, uint64_t genomesz = 0, uint64_t gCGmin = 0ULL)
        : genome_size(genomesz)
    {
        //debug_gaps = true;
        //min_gap_CG_len = gCGmin;
        //track_gaps_CG = min_gap_CG_len > 0ULL;
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
            os << "description=\"N-gaps (minimum length " << min_gap_len << ")";
            if (track_gaps_CG)
               os << " and C- and G-gaps (minimum length " << min_gap_CG_len << ")";
            os << " for " << (query.empty() ? "filename " : "sequence ") << nm << "\"";
            os << std::endl;
        }
        // BED entries for each gap, from each contig in order
        for (seqs_cI it = seqs.begin(); it != seqs.end(); ++it)
            if (query.empty() || query == it->name)
                it->gaps.print_bed(os);
    }
};  // end class FastaFileStats

}  // namespace FastaStats
