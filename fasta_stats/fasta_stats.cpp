// TODO: add detection of C- and G-gaps
// TODO: separate all-sequences and single-sequences output
// TODO: make sure can const the returned stats, must use .at() instead of []
// TODO: error if not a Fasta sequence
// DONE: use namespace std
// DONE: adjust default output filenames
// DONE: something is up with the assembly stats, maps need to be checked
// DONE: add CpG islands... make a transition matrix?
// DONE: get working! :-)

#include "SimpleOpt.h"
// the FastaStats namespace and classes
#include "FastaStats.h"
using namespace FastaStats;

using namespace std;

#include <vector>
#include <map>
#include <cstdlib>
#include <string>
#include <cctype>
#include <climits>
#include <iostream>
#include <fstream>

const string tab = "\t";
const string comma = ",";

#define DEBUG 0
static int              opt_debug = DEBUG;

static string           opt_query;                 // restrict output to only this sequence
static bool             opt_summary = true;        // produce summary table for all sequences
static bool             opt_sequences = true;      // produce output for each sequence
static bool             opt_header = true;         // add headers to tables output
static string           opt_sep = comma;           // table column separator
static string           opt_assembly;              // summary assembly stats file
const static string     assembly_default_file = "assembly.txt";
const static string     assembly_default_suffix = ".assembly.txt";
static string           opt_output;                // output file
const static string     output_default_suffix = ".stats.txt";
static bool             opt_stdin = false;
static bool             from_stdin = false;

static bool             opt_create_gaps_bed = true;
static string           opt_gaps_bed;              // gaps BED file
const static string     gaps_bed_default_file = "gaps.bed";
const static string     gaps_bed_default_suffix = ".gaps.bed";

static bool             opt_gaps_CG = false;
static uint64_t         opt_gaps_CG_min_minimum = 10ULL;
static uint64_t         opt_gaps_CG_min = opt_gaps_CG_min_minimum;

static bool             opt_assembly_stats = false;
static uint64_t         opt_genome_size = 0ULL;


// Remove path from filename
string localBasename(const string& pathname, bool remove_extension = false) {
    size_t slash = pathname.find_last_of('/');
    string filename(pathname);
    if (slash < pathname.length())  // found '/'
        filename = pathname.substr(slash + 1);
    if (remove_extension) {
        size_t dot = filename.find_last_of('.');
        if (slash < filename.length() && slash > 0)  // found '.' not at start
            filename = filename.substr(0, dot);
    }
    return filename;
}

void usage (int arg = 0)
{
    cerr << 
"USAGE:" << endl <<
endl <<
"     fasta_stats [ options ] [ fasta-file.fa ]" << endl <<
endl <<
"OPTIONS:" << endl <<
endl <<
"    -q/--query NAME     Produce stats for sequence NAME.  Default is to produce" << endl <<
"                        stats for all sequences." << endl <<
"    -/--stdin           Expect input on STDIN, write output to STDOUT" << endl <<
"                        and gaps to '" << gaps_bed_default_file << "' unless --output and/or" << endl <<
"                        --gaps-bed are specified" << endl <<
"    -o/--output FILE    Write stats output to FILE.  Default is to write to" << endl <<
"                        output file with suffix " << output_default_suffix << endl <<
"    -s/--sequences      Produce stats for individual sequences" << endl <<
"    -S/--no-sequences   Do NOT produce stats for individual sequences" << endl <<
"    -g/--gaps-bed FILE  Write BED file containing intervals for N-gaps," << endl <<
"                        if not specified defaults to output file with" << endl <<
"                        suffix " << gaps_bed_default_suffix << endl <<
"    -G/--no-gaps-bed    Do NOT write a BED file of gaps" << endl <<
"    --gaps-CG           Include intervals for C- and G-gaps in gaps BED file." << endl <<
"                        Setting --gaps-CG-min implies this option." << endl <<
"                        With this option, summary gap statistics include C- and" << endl <<
"                        G-gaps, and the 4th column of the gaps BED file indicates" << endl <<
"                        the type of gap the entry refers to" << endl <<
"    --gaps-CG-min INT   Minimum size (bp) to consider a homopolymer run of C or G" << endl <<
"                        to be a C- or G-gap; default is " << opt_gaps_CG_min_minimum << endl <<
endl <<
"    -u/--summary        Include summary stats for all sequences" << endl <<
"    -U/--no-summary     Do NOT include summary stats for all sequences" << endl <<
"    --assembly-stats    Include assembly stats in summary stats" << endl <<
"    --genome-size INT   Expected genome size, required for NG and LG" << endl <<
"                        assembly stats" << endl <<
endl <<
"    -d/--header         Write a header of column names to the output file" << endl <<
"    -D/--no-header      Do NOT write a header to the output file" << endl <<
endl <<
"    --comma             Separate output columns with commas ','" << endl <<
"    --tab               Separate output columns with tabs '\\t'" << endl <<
"    -h/-?/--help        Produce this help message" << endl <<
"    --debug INT         Debug output level INT" << endl <<
endl;
    exit(arg);
}

int main(int argc, char *argv[]) {

    enum { OPT_output,
           OPT_stdin,
           OPT_query,
           OPT_gaps_bed,       OPT_no_gaps_bed,
           OPT_gaps_CG,        OPT_gaps_CG_min,
           OPT_header,         OPT_no_header, 
           OPT_summary,        OPT_no_summary, 
           OPT_sequences,      OPT_no_sequences, 
           OPT_assembly_stats, OPT_genome_size,
           OPT_tab,            OPT_comma, 
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
        { OPT_gaps_CG,       "--gaps-CG",         SO_NONE },
        { OPT_gaps_CG_min,   "--gaps-CG-min",     SO_REQ_SEP },
        { OPT_header,        "--header",          SO_NONE },
        { OPT_header,        "-d",                SO_NONE },
        { OPT_no_header,     "--no-header",       SO_NONE },
        { OPT_no_header,     "-D",                SO_NONE },
        { OPT_summary,       "--summary",           SO_NONE },
        { OPT_summary,       "-u",                SO_NONE },
        { OPT_no_summary,    "--no-summary",        SO_NONE },
        { OPT_no_summary,    "-U",                SO_NONE },
        { OPT_sequences,     "--sequences",       SO_NONE },
        { OPT_sequences,     "-s",                SO_NONE },
        { OPT_no_sequences,  "--no-sequences",    SO_NONE },
        { OPT_no_sequences,  "-S",                SO_NONE },
        { OPT_assembly_stats,"--assembly-stats",  SO_NONE },
        { OPT_genome_size,   "--genome-size",     SO_REQ_SEP },
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
            cerr << "invalid argument " << args.OptionText() << endl;
            exit(1);
        }
        switch(args.OptionId()) {
        case OPT_help:
            usage();
            break;
        case OPT_query:
            opt_query = args.OptionArg(); break;
        case OPT_output:
            opt_output = args.OptionArg(); break;
        case OPT_stdin:
            opt_stdin = true; break;
        case OPT_gaps_CG:
            opt_gaps_CG = true; break;
        case OPT_gaps_CG_min:
            opt_gaps_CG_min = strtoull(args.OptionArg(), NULL, 10);
            if (opt_gaps_CG_min < opt_gaps_CG_min_minimum) {
                cerr << "--gaps-CG-min argument too small " << args.OptionArg() << endl;
                exit(1);
            }
            opt_gaps_CG = true;
            break;
        case OPT_header:
            opt_header = true; break;
        case OPT_no_header:
            opt_header = false; break;
        case OPT_summary:
            opt_summary = true; break;
        case OPT_no_summary:
            opt_summary = false; break;
        case OPT_sequences:
            opt_sequences = true; break;
        case OPT_no_sequences:
            opt_sequences = false; break;
        case OPT_assembly_stats:
            opt_assembly_stats = true; break;
        case OPT_genome_size:
            opt_genome_size = strtoull(args.OptionArg(), NULL, 10);
            if (opt_genome_size == 0ULL || opt_genome_size == ULLONG_MAX) {
                cerr << "--genome-size argument invalid " << args.OptionArg() << endl;
                exit(1);
            }
            break;
        case OPT_comma:
            opt_sep = comma; break;
        case OPT_tab:
            opt_sep = tab; break;
#ifdef DEBUG
        case OPT_debug:
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : 1;
            fs_debug = opt_debug; // set debugging in FastaStats namespace as well
            break;
#endif
        default:
            cerr << "unknown argument " << args.OptionText() << endl;
            exit(1);
        }
    }

    // Input and output files
    string input;
    if (args.FileCount() > 1) {
        cerr << "at most one sequence file as input" << endl;
        exit(1);
    } else if (args.FileCount() == 1) {
        input.assign(args.File(0));
    } else if (opt_stdin) {
        input = "/dev/stdin";
        from_stdin = true;
    } else {
        cerr << "at least one sequence file as input, or stdin with -/--stdin" << endl;
        exit(1);
    }
    if (opt_output.empty()) {
        if (opt_stdin) {
            opt_output = "/dev/stdout";
        } else {
            // assemble filename from input filename
            opt_output = localBasename(input, false) + output_default_suffix;
        }
    }
    if (opt_create_gaps_bed && opt_gaps_bed.empty())
        opt_gaps_bed = opt_stdin ? gaps_bed_default_file : localBasename(input, false) + gaps_bed_default_suffix;
    if (opt_debug) {
        cerr << "input:" << input << ":   output:" << opt_output << 
            ":  gaps_bed:" << opt_gaps_bed << ":" << 
            ":  gaps_CG:" << opt_gaps_CG << ":" << opt_gaps_CG_min << endl;
    }

    ofstream output(opt_output.c_str(), ofstream::out);


    ////
    ////  Start object to gather/serve Fasta file stats
    ////
    FastaFileStats fastastats(input.c_str(),
                              opt_genome_size,
                              opt_gaps_CG ? opt_gaps_CG_min : 0ULL);
    ////


    // produce output
    if (opt_debug > 1)
        fastastats.dump(cerr);

    if (! opt_query.empty()) {    // for a single query sequence

        SingleSequenceStats s = fastastats.sequence_stats(opt_query);

        s.dump(output);

        if (! opt_gaps_bed.empty())
            fastastats.create_gaps_bed(opt_gaps_bed, opt_query, opt_header);

        return 0;
    }

    if (opt_summary)
        fastastats.create_summary_table(output, opt_sep, opt_header);

    if (opt_sequences)
        fastastats.create_sequence_table(output, opt_sep, opt_header);

    if (opt_create_gaps_bed)
        fastastats.create_gaps_bed(opt_gaps_bed, "", opt_header);

	return 0;
}
