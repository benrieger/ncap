#define MAIN

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <sstream>
#include <stdlib.h> // for exit values
#include <sys/stat.h>
#include <signal.h>
using namespace std;

// include boost's program_options
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp> 
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
//using namespace boost;
namespace po = boost::program_options;

// include OpenMP to be able to set the number of threads globally
#include <omp.h>

// include own software parts
#include "../include/ngs_version.hpp"
#include "../include/ngs_utils.hpp"
#include "../include/ngs_interface.hpp"
#include "../include/ngs_filter.hpp"
#include "../include/ngs_errors.hpp"
#include "../include/ngs_global.hpp"
#include "../include/ngs_io.hpp"
#include "../include/ngs_read.hpp"
#include "../include/ngs_aligner.hpp"
#include "../include/ngs_confread.hpp"

using namespace ngs;

// no namespace, only for simple timings
#include "../include/timer.hpp"

typedef void (*signal_handler)(int);

//! internal use, only.
bool file_exists(const std::string &fname) {
    struct stat finfo;
    int fstat;
    fstat = stat(fname.c_str(), &finfo);
    if (fstat == 0) return true;
    return false;
}

void userbreak(int) {
   cerr << endl << endl << "NGSPREPROC has been stopped by you." << endl;
   exit(EXIT_SUCCESS);
}

void genericerror(int err) {
   cerr << endl << endl
        << "NGSPREPROC has crashed for an unkown reason." << endl
        << "Perhaps there is a memory problem." << endl
        << "To help improve this program, please send a bug" << endl
        << "report to meesters@uni-mainz.de" << endl
        << "Please include all data in your report, which are" << endl
        << "necessary to reproduce the crash." << endl;
   exit(err);
}

void parsingerror(const string msg) {
   cerr << endl << endl
        << "NGSPREPROC was unable to parse your commandline." << endl
        << "The parser says: " << endl
        << msg << endl;
   exit(EXIT_FAILURE);
}

void setupcrashhandlers() {
   signal(SIGINT,  (signal_handler) userbreak);
   signal(SIGSEGV, (signal_handler) genericerror);
   signal(SIGABRT, (signal_handler) genericerror);
}

void helpmessage(po::options_description &options) {
   cerr << "*------------------------------------------------*" << endl
        << "| NCAP |    ";
   cerr << tuples::set_open(' ') << tuples::set_close(' ')
        << tuples::set_delimiter('.') << VERSION;
   // consider lenght of the version tuple
   unsigned short int space = 5;
   if (VERSION.get<0>() >= 10) space--;
   if (VERSION.get<1>() >= 10) space--;
   if (VERSION.get<2>() >= 10) space--;
   for (size_t i = 0; i <= space; i++) cerr << " ";
   cerr << "|  " << DATE << "  |" << endl;
   cerr << "*------------------------------------------------*" << endl;
   cerr << "| (c) Christian Meesters                         |" << endl
        << "| Gnu Public License v3                          |" << endl
        << "*------------------------------------------------*" << endl << endl;
   cerr << options << endl;
}

struct encoding {
  encoding(std::string const& val):
    value(val)
  { }
  std::string value;
};

void validate(boost::any& v, 
              std::vector<std::string> const& values,
              encoding* /* target_type */,
              int) {
    using namespace boost::program_options;

    // Make sure no previous assignment to 'v' was made.
    validators::check_first_occurrence(v);

    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    std::string const &s = validators::get_single_string(values);
    // we need to 'un-const' to alter the case
    std::string* ps;
    ps = (std::string*)(&s);

    // convert encoding to upper case
    std::transform(ps->begin(), ps->end(), ps->begin(), ::toupper);

    if (s == "SANGER" || s == "SOLEXA" || s == "ILLUMINA") {
        v = boost::any(encoding(s));
    } 
    else {
        parsingerror("the argument for option '--encoding' is invalid");
    }
    if ( s == "SANGER" || s.size() == 0 ) ud.offset = 33;
    else ud.offset = 64;
}

int main(int args, char* argv[]) {
    #ifdef TIMER
    {
        Timer total("total_timer");
    #endif
    #ifndef DEBUG
    // before anything else, set up the crash handlers
    setupcrashhandlers();
    #endif // DEBUG
    
    bool version = false;
    string uadaptor, fadaptor, madaptor, primers, input, output, forward, reverse, configfname;
    unsigned int lboundary = 0, rboundary = 0;

    /* WARNING: Do not include duplicate option names for the shell
       and the configuration file.
    */
    // options available only on the commandline
    po::options_description general("General Command Line Options");
    general.add_options()
            ("help,h",      "produce help message")
            ("verbose,v",   po::value<bool>(&ud.verbose),
                            "if set message will be written (to stderr)")
            ("threads,t",   po::value<unsigned int>(&ud.nthreads),
                            "set number of threads to use")
            ("input,i",     po::value<string>(&input),
                            "- if omitted the programs reads from stdin (single source)"
                            " - else, a full qualified path needs to be given"
                            " - to specifiy paired read files, separate with a comma:"
                            "   <file1,file2>"
                            " - files compressed with gzip (usually ending on .gz) can be read")
            ("if",          po::value<string>(&forward),
		                    "can be used to specify the forward sequence file of a paired"
		                    " reads tuple"
		                    " (mutually exclusive with '-i,--input', requires --rf)")
            ("rf",          po::value<string>(&reverse),
		                    "can be used to specify the 'reverse' sequence file of a paired"
		                    " reads tuple"
		                    " (mutually exclusive with '-i,--input', requires --rf)")
            ("configfile,c", po::value<string>(&configfname),
                             "the path to the configfile (may be omitted)")
            ("output,o",    po::value<string>(&output),
                            "if set output will not be written to stdout, but to the designated file/path"
                            " in case of paired end reads output files should be given as comma separated strings"
                            " e.g.: \"<prefix>_1.fastq,<prefix>_2.fastq\"")
            ("minimum,m",   po::value<unsigned int>(&ud.minimum),
                            "minimum length an adaptor should fit to (defaults to 15)")
            ("length,L",    po::value<unsigned int>(&ud.length),
                            "discard remaining read(s), if their length is < l (default is 5)")
            ("uadaptors",   po::value<string>(&uadaptor),
                            "the adaptor string (universal = applied on all sequences); several adaptors may be given (comma separated)")
            ("fadaptors",   po::value<string>(&fadaptor),
                            "the adaptor string (forward = applied on all forward sequences in mate pair sequencing); several adaptors may be given (comma separated)")
            ("madaptors",   po::value<string>(&madaptor),
                            "the adaptor string (universal = applied on all sequences); several adaptors may be given (comma separated)")
            ("e_adaptor",   po::value<short unsigned int>(&ud.a_errors),
                            "number of errors an alignment may show")
            ("primers,p",   po::value<string>(&primers),
                            "the primer string, several primers may be given (comma separated)")
            ("quality,q",   po::value<unsigned int>(&ud.quality),
                            "minimum quality threshold score the remaining read should pass (defaults to 0)")
	    ("fraction,f",  po::value<double>(&ud.fraction),
	                    "minimum percentage of bases which must have -q quality")
	    ("lboundary,l", po::value<unsigned int>(&lboundary),
                            "l-value of a tuple in the form of (l,r), "
                            " where 'l' is the first base to keep (default =1)"
                            " and 'r' is the last base to keep (default = no limit)"
                            " ; see 'rboundary,r' options")
            ("rboundary,r", po::value<unsigned int>(&rboundary),
	                    "r-value of a tuple in the form of (l,r), "
                            " where 'l' is the first base to keep (default =1)"
                            " and 'r' is the last base to keep (default = no limit)"
                            " ; see 'rboundary,r' options")
            ("encoding",    po::value<encoding>(),
                            "FASTQ format encoding, defaults to \"sanger\" resulting in a +33 ASCII offset for the quality string. Allowed values are: 'sanger', 'solexa' and 'illumina' with 33, 64 and 64 offsets, respectively")
            ("allow-N",     po::value<bool>(&ud.allowN)->zero_tokens(),
                            "if set, nucleotides labelled as 'N' are not counted as errors")
            ("rc",          po::value<bool>(&ud.reverse_complement)->zero_tokens(),
                            "is set, the reverse complement will be taken for each sequence"
                            "if a paired-end file tuple is given, only the second file will"
                            "be 'reverse-complemented'")
            ("keep_comments", po::value<bool>(&ud.keep_comment)->zero_tokens(),
                            "if set, comments will be kept (takes more RAM)")
            ("collapse",    po::value<bool>(&ud.collapse)->zero_tokens(),
	                    "will discard all double reads"
                            "in case of paired end mates, the mate will also be invalided")
            ("fasta",       po::value<bool>(&ud.fasta)->zero_tokens(),
                            "if set, output will be written FASTA format")
            ("version",     po::value<bool>(&version)->zero_tokens(),
                            "if set, NGSPREPROC will print its version and exit")
           ;

    // the file to contain adaptor information
    po::options_description config("ConfigFile");
    config.add_options()
        ("forward"  ,po::value<string>(), "")
        ("mate"     ,po::value<string>(), "")
        ("universal",po::value<string>(), "")
        ;

    // add options 
    po::options_description cmdline_options;
    cmdline_options.add(general);

    po::options_description config_file_options;
    config_file_options.add(config);

    // parse command line
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(args, argv, cmdline_options), vm);
        po::notify(vm);
    }
    // catching unknown options
    catch (exception_detail::clone_impl<exception_detail::error_info_injector<program_options::unknown_option> > &msg) {
        parsingerror(msg.what());
    }
    // catching missing argument
    catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::invalid_command_line_syntax> > &msg) {
        error(msg.what());
    }
    // invalid option values
    catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::invalid_option_value> > &msg) {
	parsingerror(msg.what());
    }
    // catch multiple occurences of a command line option
    catch (exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::multiple_occurrences> > &msg) {
        string message = " as a command line option";
        parsingerror(msg.what() + message);
    }
    // before we do anything else: check whether we should report the version:
    if (version) {
        cerr << "version number: " << tuples::set_open(' ') << tuples::set_close(' ')
             << tuples::set_delimiter('.') << VERSION << endl;

        return EXIT_SUCCESS;
    }

    if (ud.a_errors >= ud.minimum) {
        parsingerror("'a_error' >= min. alignment length");
    }
    if (ud.p_errors >= ud.minimum) {
        parsingerror("'p_error' >= min. alignment length");
    }

    // check, whether there is a file name given by the adaptor option
    ifstream configfile(configfname.c_str());
    if (!configfile) {
       // we now need to extract the necessary run time info
        
         boost::algorithm::split(ud.uadaptors, uadaptor, boost::is_any_of(",;"), boost::token_compress_on);
       boost::algorithm::split(ud.fadaptors, fadaptor, boost::is_any_of(",;"), boost::token_compress_on);
       boost::algorithm::split(ud.madaptors, madaptor, boost::is_any_of(",;"), boost::token_compress_on);
       boost::algorithm::split(ud.primers, primers, boost::is_any_of(",;"), boost::token_compress_on);
    } else {
        po::store(po::parse_config_file(configfile, config_file_options), vm);
        po::notify(vm);
        ngs::conffile::conf_extractor(vm, 
                                      ud.fadaptors,
                                      ud.madaptors,
                                      ud.uadaptors);

    }
    // try dealing witht he input string
    if ( (vm.count("input") || vm.count("i")) &&  (vm.count("if") || vm.count("rf"))) {
        commandline_error("-i/--input cannot be used together with --if and/or --rf");
    }
    else if ( vm.count("input") || vm.count("i") ) {
        std::vector<std::string> tmp;
        boost::algorithm::split(tmp, input, boost::is_any_of(",;"), boost::token_compress_on);
        ud.input.get<0>() = tmp[0];
        if ( tmp.size() == 2 ) {
            ud.input.get<1>() = tmp[1];
          
        }
        else if ( tmp.size() > 2 ) {
            parsingerror("unable to process more than 2 input files");
        }
    }
    else if (vm.count("if")) {
		ud.input.get<0>() = forward;
		if (!vm.count("rf")) {
			commandline_error("--if needs --rf as its counterpard");
		}
		else ud.input.get<1>() = reverse;
    }

    void (*collapsefunc) (vector<Read> &rv, vector<Read> &ov) = NULL;
    if (vm.count("collapse")) {
        if (ud.input.get<1>().size() > 0 ) {
             collapsefunc = &collapse;
        } else {
             collapsefunc = &collapse_paired_end;
        }
    }

    if ( vm.count("output") || vm.count("o") ) {
       std::vector<string> tmp;
       boost::algorithm::split(tmp, output, boost::is_any_of(",;"), boost::token_compress_on);
       if ( tmp.size() > 0 ) {
           ud.output.get<0>() = tmp[0];
       }
       if ( tmp.size() == 2 ) {
           ud.output.get<1>() = tmp[1];
       }
       else if ( tmp.size() > 2 ) {
          parsingerror("unable to process more than 2 output files");
       }
    }
    else if ( ud.input.get<1>().size() > 0 ) { // multiple inputs, but no ouput string?
        boost::tuple<string, string> items = split_filename( ud.input.get<0>() );
        stringstream ss;
        ss << items.get<0>() << "_filtered." << items.get<1>();
        ud.output.get<0>() = ss.str();
        if (ud.input.get<0>().size() == 2) { // 2nd input string given?
            boost::tuple<string, string> items = split_filename( ud.input.get<1>() );
            stringstream ss;
            ss << items.get<0>() << "_filtered." << items.get<1>();
            ud.output.get<0>() = ss.str();
        }
    }

    // checking whether we are dealing with 0, 1, or more adaptors
    // but first we distinguish between forward and reverse reads, 
    // if applicable:
    ud.fadaptors.insert(ud.fadaptors.end(), ud.uadaptors.begin(), ud.uadaptors.end());
    ud.madaptors.insert(ud.madaptors.end(), ud.uadaptors.begin(), ud.uadaptors.end());
    
    void (*fadaptorfunc) (Read &r, const std::vector<std::string> &a) = NULL;
    void (*madaptorfunc) (Read &r, const std::vector<std::string> &a) = NULL;
    // nothing given?
    if (ud.fadaptors[0].size() == 0 && ud.fadaptors.size() == 1 &&
        ud.madaptors[0].size() == 0 && ud.madaptors.size() == 1) {
	fadaptorfunc = &null_func_vec;
        madaptorfunc = &null_func_vec;
    }
    // just forward adaptor(s) given?
    else if (ud.fadaptors[0].size() > 0 &&
             ud.madaptors[0].size() == 0 && ud.madaptors.size() == 1) {
	fadaptorfunc = &trim_adaptors;
	madaptorfunc = &null_func_vec;
    }
    // just one mate adaptor(s) given?
    else if (ud.fadaptors[0].size() == 0 && ud.fadaptors.size() == 1 &&
             ud.madaptors[0].size() > 0 ) {
        fadaptorfunc = &null_func_vec;
        madaptorfunc = &trim_adaptors;
    }
    else {
        fadaptorfunc = &trim_adaptors;
        madaptorfunc = &trim_adaptors;
    } 
    // repeating for primers
    void (*primerfunc) (Read &r) = NULL;
    if (ud.primers[0].size() == 0 && ud.primers.size() == 1) {
	primerfunc = &null_func;
    }
    else if (ud.primers[0].size() > 0 && ud.primers.size() == 1) {
        ud._primer = ud.primers[0];
	primerfunc = &trim_primer;
    }
    else {
        primerfunc = &trim_primers;
    } 

    // setting boundaries
    void (*trimfunc) (Read &r) = NULL;
    if ( (vm.count("lboundary") || vm.count("l")) && 
         (vm.count("rboundary") || vm.count("r")) ) { 
        ud.boundaries.get<0>() = lboundary; 
        ud.boundaries.get<1>() = rboundary;
        trimfunc = &trimlr;}
    else if ( vm.count("lboundary") || vm.count("l") ) { 
        ud.boundaries.get<0>() = lboundary;
        trimfunc = &trimr;
    }
    else if ( vm.count("rboundary") || vm.count("r") ) { 
        ud.boundaries.get<1>() = lboundary;
        trimfunc = &trimr;
    }
    else trimfunc = &null_func;

    if (vm.count("help") || vm.count("h")) {
        helpmessage(general);
        exit(EXIT_SUCCESS);
    }

    // before everything is implemented, inmplementation checks are placed, here
    if (ud.allowN) not_implemented_error("Sorry, considering N's is not yet implemented.");

    // set thread information
    omp_set_num_threads(ud.nthreads);
    omp_set_nested(true);

    // set the output formatter
    void (*outfunc)(const vector<Read> &r, ostream &o, ostream &e) = NULL;
    if ( ud.fasta ) {
        outfunc = &write_fasta;
    } else {
        outfunc = &write_reads;
    }

    if (!ud.collapse) read_filter_parallelio(fadaptorfunc, madaptorfunc, primerfunc, outfunc, trimfunc);
    else
    read_filter(fadaptorfunc, madaptorfunc, primerfunc, outfunc, trimfunc, collapsefunc);
    #ifdef TIMER
    }
    #endif

    return EXIT_SUCCESS;
}  

