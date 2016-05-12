#include <vector>
#include <string>
#include <ostream>
#include <fstream>

#include <boost/tuple/tuple.hpp>

#ifndef NGS_INTERFACE
#define NGS_INTERFACE

namespace ngs {
    //! default value for left hand boundary
    const unsigned int left = 1;
    //! default value for right hand boundary
    const unsigned int right = 10000000;
    /*! This struct with public members, only, is to provide easy
        accessible flags globally.
    */
    struct UsageData {
        //! verbose on/off
        bool verbose;
        //! number of threads
        unsigned int nthreads;
        //! file name(s) for input (2 in case of paired ends)
        boost::tuple<std::string, std::string> input;
        //! strings for output - if omitted std::cout will be used
        boost::tuple<std::string, std::string> output;
        //! minimum overlap length for an adaptor
        unsigned int minimum;
        //! minimum length of the remaining read
        unsigned int length;
        //! univsersal adaptors
        std::vector<std::string> uadaptors;
        //! forward adaptors
        std::vector<std::string> fadaptors;
        //! reverse adaptors
        std::vector<std::string> madaptors;
        //! number of allowed errors for a given alignment
        short unsigned int a_errors;
        //! primers 5'
        std::vector<std::string> primers;
        //! internal use: for setting only one primer and avoiding loops
        std::string _primer;
        //! number of allowed errors for a given alignment
        short unsigned int p_errors;
        //! quality threshold
        unsigned int quality;
	//! quality fraction threshold
	double fraction;
	//! boundaries for trimming
	boost::tuple<unsigned int, unsigned int> boundaries;
        //! whether or not Ns should be allowed in the sequence
        bool allowN;
        //! whether or not the reverse complement should be calculated
        bool reverse_complement;
        //! whether or not the output shall be in FASTA format
        bool fasta;
        //! whether or not the original comments should be kept
        bool keep_comment;
	//! whether or not identical reads (seq identical) should be kept
	bool collapse;
        //! ASCII encoding offset
        short int offset;

        UsageData(): verbose(false),
                     nthreads(1),
                     input("", ""),
                     output("", ""),
                     minimum(15),
                     length(5),
                     uadaptors(std::vector<std::string>()),
                     fadaptors(std::vector<std::string>()),
                     madaptors(std::vector<std::string>()),
                     a_errors(0),
                     primers(std::vector<std::string>()),
                     quality(0),
                     fraction(0),
                     boundaries(left, right),
                     allowN(false),
                     reverse_complement(false),
                     fasta(false),
                     keep_comment(false),
                     collapse(false),
                     offset(33)
                     {};
    };
} // end namespace ngs
    
#endif // NGS_INTERFACE
