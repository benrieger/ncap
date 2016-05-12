#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <ctype.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/regex.hpp>

#include "ngs_errors.hpp"

#ifndef NGS_IO_UTILS
#define NGS_IO_UTILS

namespace io_utils {
    boost::regex valid_qual("\\S");
 
    unsigned long int fastq_sniff(const std::string &fname) {
        struct stat fileInfo;
        unsigned long int line_count = 0;
        unsigned long int nbytes     = 0;

        std::ifstream istream;
	boost::iostreams::filtering_stream<boost::iostreams::input> inputfile;

	// dealing with a gzipped file?
	if (fname.find("gz") != std::string::npos) {
	    istream.open(fname.c_str());
	    inputfile.push(boost::iostreams::gzip_decompressor());
	    inputfile.push(istream);
            inputfile.unsetf(std::ios_base::skipws);
            line_count = std::count(std::istream_iterator<char>(inputfile), std::istream_iterator<char>(), '\n');
	}
	else {
            // combine the stat object and filename
	    if (stat(fname.c_str(), &fileInfo) != 0 ) {
                ngs::ioerror(fname);
            }
	    istream.open(fname.c_str());
	    inputfile.push(istream);
            nbytes = fileInfo.st_size;
	}

	// start the sniffing
	std::string id_line, seq_line, comment_line, quality_line;

	std::getline(inputfile, id_line);
	std::getline(inputfile, seq_line);
	std::getline(inputfile, comment_line);
	std::getline(inputfile, quality_line);

	// checking whether sequence and quality characters equal
	if (seq_line.size() != quality_line.size()) 
	    ngs::formating_error(fname, "sequence length != quality length (read 1)");
	// TODO: check for only ATCG or ATCG/N in seq
	// TODO: check for @ in seq start and + in comment start
	// checking whether all characters in quality string are alphanumeric
	boost::cmatch what;
        if (boost::regex_match(quality_line.c_str(), what, valid_qual)) {
	   ngs::formating_error(fname, "whitespace character in quality string");
	}

	size_t all = id_line.size() + seq_line.size() + comment_line.size() + quality_line.size();
        // catch the case that a file did not exist or contained empty lines
        if ( fileInfo.st_size == 0 || all == 0 ) {
            ngs::error("Problem with the input file. Aborting.\nEmpty lines?");
        }

	// dealing with a gzipped file?
	if (fname.find("gz") != std::string::npos) {
            return line_count;
        }
	return nbytes / all;
    }

} //~ io_utils

#endif 
