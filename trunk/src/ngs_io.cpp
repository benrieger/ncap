#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace boost;

#include <omp.h>

// include own software parts
#include "../include/ngs_utils.hpp"
#include "../include/ngs_io.hpp"
#include "../include/ngs_io_utils.hpp"
#include "../include/ngs_read.hpp"
#include "../include/ngs_errors.hpp"
#include "../include/ngs_interface.hpp"
#include "../include/ngs_global.hpp"
using namespace ngs;

#include <omp.h>

EXTERN ngs::UsageData ud;

namespace ngs {
    bool read_lines(infiletype &forwardfile,
                    infiletype &reversefile,
                    std::vector<std::shared_ptr<std::string> > &forward,
                    std::vector<std::shared_ptr<std::string> > &reverse,
                    size_t nlines) {
        if (nlines % 4 != 0) ngs::error("nlines must be divisible by 4");
        size_t nforward = 0, nreverse = 0;
        std::string fline, rline;
        while (nforward < nlines && nreverse < nlines) {
            std::getline(forwardfile, fline);
            std::getline(reversefile, rline);
            if (__builtin_expect(fline.empty() || rline.empty(), 0)) break;
            forward.push_back(std::make_shared<std::string>(std::move(fline)));
            reverse.push_back(std::make_shared<std::string>(std::move(rline)));
            nforward++;
            nreverse++;
        }
        
        return forwardfile.eof() && reversefile.eof();
    }

    bool read_lines(infiletype &forwardfile,
                    std::vector<std::shared_ptr<std::string> >&forward,
                    size_t nlines) {
        if (nlines % 4 != 0) ngs::error("nlines must be divisible by 4");
        size_t nforward = 0;
        std::string fline;
        bool eof_reached = false;
        while (nforward < nlines) {
            std::getline(forwardfile, fline);
            if (__builtin_expect(fline.empty(), 0)) {
                eof_reached = true;
                break;
            }
            forward.push_back(std::make_shared<std::string>(std::move(fline)));
            nforward++;
        }

        return eof_reached;
    }

    void read_reads(std::vector<ngs::Read> &forward,
                    std::vector<ngs::Read> &reverse) {
        // open input stream depending on the name suffix
        //std::ifstream istream; // note: stream must be defined before(!) filtering_stream!!!
        //boost::iostreams::filtering_stream<boost::iostreams::input> inputfile;
        // declare lines associated with reads:
        std::string id_line;
        std::string seq_line;
        std::string comment_line ; //= "+";
        std::string quality_line;

        // temporary vector for reads of this file
        std::vector<ngs::Read> vr;
        // internal counter
        unsigned int nreads;

        std::vector<std::string> temp_input;
        temp_input.push_back(ud.input.get<0>());
        if ( ud.input.get<1>().size() > 0 )
            temp_input.push_back(ud.input.get<1>());
		
        #pragma omp parallel for private(id_line,seq_line,comment_line,quality_line,nreads,vr)
        for (size_t i = 0; i < temp_input.size(); i++) {
            // guess number of reads to be read and validate format
            unsigned long int reads_in_file = 50000000; // just a guess
            // can be empty, if reading from stdin - but than sniffing is not possible
            if (! temp_input[i].empty() )
                reads_in_file = io_utils::fastq_sniff(temp_input[i]);

            std::ifstream istream; // note: stream must be defined before(!) filtering_stream!!!
            boost::iostreams::filtering_stream<boost::iostreams::input> inputfile;
            //#pragma omp critical(file_opening)
            {
            if (temp_input[i].find(".gz") != std::string::npos) { // dealing with a gzipped file?
                istream.open(temp_input[i].c_str());
                if (!istream) ngs::ioerror(temp_input[i]);
                inputfile.push(boost::iostreams::gzip_decompressor());
                inputfile.push(istream);
            }
            else if ( temp_input[i].empty()  ) {
                inputfile.push(std::cin);
            }
            else {
                istream.open(temp_input[i].c_str());
                if (!istream) ngs::ioerror(temp_input[i]);
                inputfile.push(istream);
            }
            } //~file_opening
            // we count the number of lines in this file
			//size_t lcount = std::count(std::istreambuf_iterator<char>(inputfile), std::istreambuf_iterator<char>(), '\n');
			// jump back to the beginning
			//inputfile.seekg(0, std::ios_base::beg);
			
            // we start creating space for the temporary vector
            //vr.reserve(lcount); // 5e6 
            vr.reserve(reads_in_file);

            nreads = 0;

            // walk over the file
            while (true) {
                std::getline(inputfile, id_line);
                std::getline(inputfile, seq_line);
                //if (ud.keep_comment) std::getline(inputfile, comment_line);
                std::getline(inputfile, comment_line);
                std::getline(inputfile, quality_line);
                if (inputfile.eof()) break;
		// TODO: test whether correct or not
                if (nreads++ >= vr.capacity()) vr.resize(vr.size()+5000000);
                vr.push_back( ngs::Read(id_line, 
                                        seq_line, 
                                        comment_line, 
                                        quality_line));


            }
            vr.shrink_to_fit();
            if (i == 0 ) forward.swap(vr);
            else if (i == 1) reverse.swap(vr);
        } //for input
    } //~read_reads

    void write_reads(const std::vector<ngs::Read> &reads, std::ostream &out, std::ostream &err) {
        std::vector<ngs::Read>::const_iterator rit;
        for (rit = reads.begin(); rit != reads.end(); rit++) {
            if ((*rit).is_valid()) {
                out << (*rit);
            }
            else {
                err << (*rit);
            }
        }
        out.flush();
        err.flush();
    }

    void write_fasta(const std::vector<ngs::Read> &reads, std::ostream &out, std::ostream &err) {
        std::vector<ngs::Read>::const_iterator rit;
        for (rit = reads.begin(); rit != reads.end(); rit++) {
            if ( (*rit).is_valid() )
                out << "> " << (*rit).get_id() << '\n' 
                    << (*rit).get_seq() << '\n';
            else
                err << "> " << (*rit).get_id() << ": " 
                    << (*rit).get_comment() << '\n'
                    << (*rit).get_seq() << '\n';

        }
        out.flush();
        err.flush();
    }
} //~ngs
