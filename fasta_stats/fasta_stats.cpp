// TODO: separate all-sequences and single-sequences output
// TODO: make sure can const the returned stats, must use .at() instead of []
// TODO: error if not a Fasta sequence
// DONE: adjust default output filenames
// DONE: something is up with the assembly stats, maps need to be checked
// DONE: add CpG islands... make a transition matrix?
// DONE: get working! :-)

#include "SimpleOpt.h"
// the FastaStats namespace and classes
#include "FastaStats.h"
using namespace FastaStats;

#include <vector>
#include <map>
#include <cstdlib>
#include <string>
#include <cctype>
#include <climits>
#include <iostream>
#include <fstream>

const std::string tab = "\t";
const std::string comma = ",";

#define DEBUG 0
static int         opt_debug = DEBUG;

static std::string opt_query;                 // restrict output to only this sequence
static bool        opt_summary = true;        // produce summary table for all sequences
static bool        opt_sequences = true;      // produce output for each sequence
static bool        opt_header = true;         // add headers to tables output
static std::string opt_sep = comma;           // table column separator
static std::string opt_assembly;              // summary assembly stats file
const static std::string assembly_default_file = "assembly.txt";
const static std::string assembly_default_suffix = ".assembly.txt";
static std::string opt_output;                // output file
const static std::string output_default_suffix = ".stats.txt";
static bool        opt_stdin = false;
static bool        from_stdin = false;
static std::string opt_gaps_bed;              // gaps BED file
const static std::string gaps_bed_default_file = "gaps.bed";
const static std::string gaps_bed_default_suffix = ".gaps.bed";
static bool        opt_create_gaps_bed = true;
static bool        opt_assembly_stats = false;
static uint64_t    opt_genome_size = 0ULL;


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

void usage (int arg = 0)
{
    std::cerr << 
"USAGE:" << std::endl <<
"" << std::endl <<
"     fasta_stats [ options ] [ fasta-file.fa ]" << std::endl <<
"" << std::endl <<
"OPTIONS:" << std::endl <<
"" << std::endl <<
"    -q/--query NAME     Produce stats for sequence NAME.  Default is to produce" << std::endl <<
"                        stats for all sequences." << std::endl <<
"    -/--stdin           Expect input on STDIN, write output to STDOUT" << std::endl <<
"                        and gaps to '" << gaps_bed_default_file << "' unless --output and/or" << std::endl <<
"                        --gaps-bed are specified" << std::endl <<
"    -o/--output FILE    Write stats output to FILE.  Default is to write to" << std::endl <<
"                        'inputfilename.stats.fa'." << std::endl <<
"    -s/--sequences      Produce stats for individual sequences" << std::endl <<
"    -S/--no-sequences   Do NOT produce stats for individual sequences" << std::endl <<
"    -g/--gaps-bed FILE  Write BED file containing intervals for N-gaps," << std::endl <<
"                        if not specified defaults to output file with" << std::endl <<
"                        suffix " << gaps_bed_default_suffix << std::endl <<
"    -G/--no-gaps-bed    Do NOT write a BED file of N-gaps" << std::endl <<
"" << std::endl <<
"    -u/--summary        Include summary stats for all sequences" << std::endl <<
"    -U/--no-summary     Do NOT include summary stats for all sequences" << std::endl <<
"    --assembly-stats    Include assembly stats in summary stats" << std::endl <<
"    --genome-size INT   Expected genome size, required for NG and LG" << std::endl <<
"                        assembly stats" << std::endl <<
"" << std::endl <<
"    -d/--header         Write a header of column names to the output file" << std::endl <<
"    -D/--no-header      Do NOT write a header to the output file" << std::endl <<
"" << std::endl <<
"    --comma             Separate output columns with commas ','" << std::endl <<
"    --tab               Separate output columns with tabs '\\t'" << std::endl <<
"    -h/-?/--help        Produce this help message" << std::endl <<
"    --debug INT         Debug output level INT" << std::endl <<
"" << std::endl;
    exit(arg);
}

int main(int argc, char *argv[]) {

    enum { OPT_output,
           OPT_stdin,
           OPT_query,
           OPT_gaps_bed,       OPT_no_gaps_bed,
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
            std::cerr << "invalid argument " << args.OptionText() << std::endl;
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
        case OPT_gaps_bed:
            opt_gaps_bed = args.OptionArg();
            opt_create_gaps_bed = true;
            break;
        case OPT_no_gaps_bed:
            opt_create_gaps_bed = false; break;
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
                std::cerr << "--genome-size argument invalid " << args.OptionArg() << std::endl;
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
            opt_output = localBasename(input, false) + output_default_suffix;
        }
    }
    if (opt_create_gaps_bed && opt_gaps_bed.empty())
        opt_gaps_bed = opt_stdin ? gaps_bed_default_file : localBasename(input, false) + gaps_bed_default_suffix;
    if (opt_debug)
        std::cerr << "input:" << input << ":   output:" << opt_output << 
            ":  gaps.bed:" << opt_gaps_bed << ":" << std::endl;

    std::ofstream output(opt_output.c_str(), std::ofstream::out);

    FastaFileStats fastastats(input.c_str(), opt_genome_size);

    // produce output
    if (opt_debug > 1)
        fastastats.dump(std::cerr);
    if (! opt_query.empty()) {
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
