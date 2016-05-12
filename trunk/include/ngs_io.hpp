#include <vector>
#include <string>
#include <memory>

#include <boost/iostreams/filtering_stream.hpp>
using namespace boost;

#ifndef NGS_IO
#define NGS_IO

namespace ngs {
    // forward declaration
    class Read;

    //! short hand for input filtering streams
    typedef boost::iostreams::filtering_stream<boost::iostreams::input> infiletype;

    //! a function to read a given number of lines
    bool read_lines(infiletype &forwardfile,
                    infiletype &reversefile,
                    std::vector<std::shared_ptr<std::string> > &forward,
                    std::vector<std::shared_ptr<std::string> > &reverse,
                    size_t nlines);

    //! a function to read a given number of lines -- only forward
    bool read_lines(infiletype &forwardfile,
                    std::vector<std::shared_ptr<std::string> > &forward,
                    size_t nlines);

    //! a function to read Reads
    void read_reads(std::vector<Read> &forward,
                    std::vector<Read> &reverse);
    //! a function to write Reads
    void write_reads(const std::vector<Read> &reads,
                     std::ostream &out, 
                     std::ostream& err);
    //! a function to write in FASTA format
    void write_fasta(const std::vector<Read> &reads, 
                     std::ostream &out, 
                     std::ostream& err);
} //~ngs

#endif
