#include <vector>
#include <string>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace boost;

#include <omp.h>

// include own software parts
#include "../include/ngs_global.hpp"
#include "../include/ngs_read.hpp"
#include "../include/ngs_interface.hpp"
#include "../include/ngs_io.hpp"
#include "../include/ngs_errors.hpp"
#include "../include/ngs_statistics.hpp"

#include "../include/timer.hpp"

namespace ngs {
    const size_t READ_AHEAD = 100000;
    void read_filter_parallelio(
                     void (*fadaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*madaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*primer_func) (Read &r),
                     void (*out_func)(const std::vector<Read>& r, 
                                     std::ostream &o, std::ostream &e),
                     void (*trim_func)(Read &r)
                     ) {
        // 1st: as it is not possible to copy ostreams, we need to 
        //      prepare them locally:
        std::ofstream out1, out2;
        stringstream serr1, serr2;
        std::ofstream err1, err2;
        if (ud.output.get<1>().size() == 0 && ud.input.get<1>().size() > 0 ) {
            std::stringstream msg;
            msg << "Unable to write to unspecified output, while two "
                << "input files were assigned." << std::endl;
            error(msg.str());
        }
        else if ( ud.output.get<0>().size() > 0 ) { // a filename?
            out1.open( ud.output.get<0>(), std::ofstream::out );
            // we need to suppress flushing at every line
            //out1.put(out1.widen('\n'));
            serr1 << ud.output.get<0>() << "_invalid.fastq";
            err1.open( serr1.str(), std::ofstream::out );
            //err1.put(err1.widen('\n'));
        } else {
            out1.copyfmt(std::cout);
            out1.clear(std::cout.rdstate());
            out1.basic_ios<char>::rdbuf(std::cout.rdbuf());
        }
        if ( ud.output.get<1>().size() > 0 ) { // a filename?
            out2.open( ud.output.get<1>(), std::ofstream::out );
            // we need to suppress flushing at every line
            //out2.put(out1.widen('\n'));
            serr2 << ud.output.get<1>() << "_invalid.fastq";
            err2.open( serr2.str(), std::ofstream::out );
            //err2.put(err2.widen('\n'));
        } else {
            out2.copyfmt(std::cout);
            out2.clear(std::cout.rdstate());
            out2.basic_ios<char>::rdbuf(std::cout.rdbuf());
        }

        // counter pointers: 1 per thread
        std::vector<Counter*> counter;
        counter.reserve(omp_get_num_threads());
        Counter* cit;
        for (int i = 0; i < omp_get_max_threads(); i++) { 
            //cit = Counter::getCounter(forward[0].size());
            cit = Counter::getCounter();
            counter.push_back(cit);
        };
        // global stats model
        GlobalStats* sp = GlobalStats::get_Statistics();

        bool   eof_reached  = false; // did we hit the end of the
                                     // input files?
        bool   paired_input = false; // dealing with two input files?
        std::ifstream forward_istream, reverse_istream;
        // boost::iostreams::filtering_stream<boost::iostreams::input>
        infiletype forwardfile, reversefile;
        if (ud.input.get<0>().find(".gz") != std::string::npos) { // dealing with a gzipped file?
            forward_istream.open(ud.input.get<0>().c_str());
            if (!forward_istream) ngs::ioerror(ud.input.get<0>());
            forwardfile.push(boost::iostreams::gzip_decompressor());
            forwardfile.push(forward_istream);
        } else {
            forward_istream.open(ud.input.get<0>().c_str());
            if (!forward_istream) ngs::ioerror(ud.input.get<0>());
            forwardfile.push(forward_istream);
        }

        if (! ud.input.get<1>().empty() and ud.input.get<1>().find(".gz") != std::string::npos) { // dealing with a gzipped file?
            reverse_istream.open(ud.input.get<1>().c_str());
            if (!reverse_istream) ngs::ioerror(ud.input.get<1>());
            reversefile.push(boost::iostreams::gzip_decompressor());
            reversefile.push(reverse_istream);
            paired_input = true;
        } else if (! ud.input.get<1>().empty() ) {
            reverse_istream.open(ud.input.get<1>().c_str());
            if (!reverse_istream) ngs::ioerror(ud.input.get<1>());
            reversefile.push(reverse_istream);
            paired_input = true;
        }
        
        // we bind to the appropriate filter function
        typedef void (Counter::*quality_filter) (Read&, std::vector<Read>*, size_t) ;
        quality_filter q = NULL;
        if ( paired_input ) {
            q = &Counter::qual_fil_paired_end;
        } else {
            q = &Counter::qual_fil;
        }
        
        // we bind an appropriate function for reverse complement calculations
        // it only makes sense, if at all, for the mate
        typedef void (Read::*reverse_complement_func)();
        reverse_complement_func rvcfunc_reverse = NULL;
        reverse_complement_func rvcfunc_forward = NULL;
        if (ud.reverse_complement && paired_input) {
            rvcfunc_reverse = &Read::reverse_complement;
            rvcfunc_forward = &Read::nullfunc;
        }
        else if (ud.reverse_complement && not paired_input) {
            rvcfunc_forward = &Read::reverse_complement;
            rvcfunc_reverse = &Read::nullfunc;
        }
        else {
            rvcfunc_forward = &Read::nullfunc;
            rvcfunc_reverse = &Read::nullfunc;
        }

        #pragma omp parallel shared(forwardfile, reversefile, eof_reached) private(cit)
        {
            int my_thread = omp_get_thread_num();
            // construct container to hold reads
            std::vector<Read> forward, reverse;
            std::vector<std::shared_ptr<std::string> > forward_lines, reverse_lines;
            while (! eof_reached ) {
                forward_lines.reserve(READ_AHEAD);
                forward.reserve(READ_AHEAD);

                #pragma omp critical(reading)
                {
                if ( paired_input ) {
                    reverse_lines.reserve(READ_AHEAD);
                    reverse.reserve(READ_AHEAD);

                    eof_reached = read_lines(forwardfile,
                                             reversefile,
                                             forward_lines,
                                             reverse_lines,
                                             READ_AHEAD);
                } else {
                    eof_reached = read_lines(forwardfile,
                                             forward_lines,
                                             READ_AHEAD);
                }
                }
                // produce reads from lines
                for (auto lit = forward_lines.begin(); lit != forward_lines.end(); lit+=4 ) {
                    {
                    forward.push_back(Read(*(lit),
                                           *(lit+1),
                                           *(lit+2),
                                           *(lit+3)));
                    }
                }
                forward_lines.erase(forward_lines.begin(), forward_lines.end());
                
                size_t position;
                for (auto rit = forward.begin(); 
                     rit < forward.end(); rit++) {
                     // calculate the position within the vector
                     // to be used to access the other vector, if necessary for invalidation
                     position = rit - forward.begin();
                     // get a 'threadprivate' counter
                     cit = counter[my_thread];
                     // clip adaptors, if applicable
                     fadaptor_func(*rit, ud.fadaptors);
                     // clip primers, if applicable
                     primer_func(*rit); 
                     // trimm, if asked for
                     trim_func(*rit);
                     // reverse complement, if applicable
                     (*rit.*rvcfunc_forward)();
                     // adapt gathered limits
                     rit->adapt_limits();
                     // gather stats and filter
                     (cit->*q)(*rit, &reverse, position);
                }
                    
                // produce reads from lines
                for (auto lit = reverse_lines.begin(); lit != reverse_lines.end(); lit+=4 ) {
                    reverse.push_back(Read(*lit,
                                           *(lit+1),
                                           *(lit+2),
                                           *(lit+3)));
                }
                reverse_lines.erase(reverse_lines.begin(), reverse_lines.end());

                for (auto rit = reverse.begin(); 
                     rit < reverse.end(); rit++) {
                     // calculate the position within the vector
                     // to be used to access the other vector, if necessary for invalidation
                     position = rit - reverse.begin();
                     // get a 'threadprivate' counter
                     cit = counter[my_thread];
                     // clip adaptors, if applicable
                     madaptor_func(*rit, ud.madaptors);
                     // clip primers, if applicable
                     primer_func(*rit);
                     // reverse complement, if applicable
                     (*rit.*rvcfunc_reverse)();
		     // trimm, if asked for
		     trim_func(*rit);
                     // adapt gathered limits
                     rit->adapt_limits();
                     // gather stats and filter
                     (cit->*q)(*rit, &forward, position);
                }

                #pragma omp critical(writing_output)
                {
                //std::cout << forward.size() << std::endl;
                for (size_t i = 0; i < 2; i++) {
                    if (i == 0) {
                        // write output for forward mate
                        out_func(forward, out1, err1);
                    }
                    if (i == 1 && paired_input ) {
                        // write ouput for reverse mate
                        out_func(reverse, out2, err2);
                    }
                }
                } //~ critical:writing_output
                
                // save memory
                forward.erase(forward.begin(), forward.end());
                if ( paired_input ) {
                   reverse.erase(reverse.begin(), reverse.end());
                }
            } //~ if not eof_reached
        } // ~ end of parallel section
            
        // update and delete counter resources
        for (auto it = counter.begin(); it != counter.end(); it++) {
            (*sp) += *(*it);
            delete *it;
        }
        counter.clear();
        
   try {
        out1.close();
    } catch (...) {};
    try {
        out2.close();
    } catch (...) {};
    try {
      err1.close();
    } catch (...) {};
    try {
      err2.close();
    } catch (...) {};

    } //~read_filter_parallelio

    void read_filter(
                     void (*fadaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*madaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*primer_func) (Read &r),
                     void (*out_func)(const std::vector<Read>& r, 
                                     std::ostream &o, std::ostream &e),
                     void (*trim_func)(Read &r),
                     void (*collapse_func) (std::vector<Read> &rv,
                                            std::vector<Read> &ov) 
                     ) {
        // 1st: as it is not possible to copy ostreams, we need to 
        //      prepare them locally:
        std::ofstream out1, out2;
        stringstream serr1, serr2;
        std::ofstream err1, err2;
        if (ud.output.get<1>().size() == 0 && ud.input.get<1>().size() > 0 ) {
            std::stringstream msg;
            msg << "Unable to write to unspecified output, while two "
                << "input files were assigned." << std::endl;
            error(msg.str());
        }
        else if ( ud.output.get<0>().size() > 0 ) { // a filename?
            out1.open( ud.output.get<0>(), std::ofstream::out );
            // we need to suppress flushing at every line
            //out1.put(out1.widen('\n'));
            serr1 << ud.output.get<0>() << "_invalid.fastq";
            err1.open( serr1.str(), std::ofstream::out );
            //err1.put(err1.widen('\n'));
        } else {
            out1.copyfmt(std::cout);
            out1.clear(std::cout.rdstate());
            out1.basic_ios<char>::rdbuf(std::cout.rdbuf());
        }
        if ( ud.output.get<1>().size() > 0 ) { // a filename?
            out2.open( ud.output.get<1>(), std::ofstream::out );
            // we need to suppress flushing at every line
            //out2.put(out1.widen('\n'));
            serr2 << ud.output.get<1>() << "_invalid.fastq";
            err2.open( serr2.str(), std::ofstream::out );
            //err2.put(err2.widen('\n'));
        } else {
            out2.copyfmt(std::cout);
            out2.clear(std::cout.rdstate());
            out2.basic_ios<char>::rdbuf(std::cout.rdbuf());
        }
        // 1st: construct container to hold reads
        std::vector<Read> forward, reverse;
        // 2nd: read in the files
        #ifdef TIMER
        {
            Timer io_timer("read_timer");
        #endif
        read_reads(forward, reverse);
        #ifdef TIMER
        }
        #endif

        if (forward.size() != reverse.size() && ud.input.get<1>().size() > 0 ) {
            ngs::error("Input files for paired reads are of different size.");
        }

        #ifdef TIMER
        {
            Timer process_timer("process_timer");
        #endif

        // we bind to the appropriate filter function
        typedef void (Counter::*quality_filter) (Read&, std::vector<Read>*, size_t) ;
        quality_filter q = NULL;
        if ( ud.input.get<1>().size() > 0 ) {
            q = &Counter::qual_fil_paired_end;
        } else {
            q = &Counter::qual_fil;
        }

        // we bind an appropriate function for reverse complement calculations
        typedef void (Read::*reverse_complement_func)();
        reverse_complement_func rvcfunc_forward = NULL;
        reverse_complement_func rvcfunc_reverse = NULL;
        if (ud.reverse_complement && ud.input.get<1>().size() > 0) {
            rvcfunc_reverse = &Read::reverse_complement;
            rvcfunc_forward = &Read::nullfunc;
        }
        else if (ud.reverse_complement && ud.input.get<1>().size() == 0 ) {
            rvcfunc_forward = &Read::reverse_complement;
            rvcfunc_reverse = &Read::nullfunc;
        }
        else {
            rvcfunc_forward = &Read::nullfunc;
            rvcfunc_reverse = &Read::nullfunc;
        }

        // counter pointers: 1 per thread
        std::vector<Counter*> counter;
        counter.reserve(omp_get_num_threads());
        Counter* cit;
        for (int i = 0; i < omp_get_max_threads(); i++) { 
            //cit = Counter::getCounter(forward[0].size());
            cit = Counter::getCounter();
            counter.push_back(cit);
        };
        // global stats model
        GlobalStats* sp = GlobalStats::get_Statistics();

        // 3rd: in case of paired reads - open sections
        #pragma omp parallel
        {
            // we need to know the position for quality checks
            size_t position;
                
            #pragma omp for schedule(static) private(position, cit) nowait
            for (auto rit = forward.begin(); 
                 rit < forward.end(); rit++) {
                 // calculate the position within the vector
                 // to be used to access the other vector, if necessary for invalidation
                 position = rit - forward.begin();
                 // get a 'threadprivate' counter
                 cit = counter[omp_get_thread_num()];
                 // reverse complement, if applicable
                 (*rit.*rvcfunc_forward)();
                 // clip adaptors, if applicable
                 fadaptor_func(*rit, ud.fadaptors);
                 // clip primers, if applicable
                 primer_func(*rit); 
                 // trimm, if asked for
                 trim_func(*rit);
                 // adapt gathered limits
                 rit->adapt_limits();
                 // gather stats and filter
                 (cit->*q)(*rit, &reverse, position);
            }

            #pragma omp for schedule(static) private(position, cit)  
            for (auto rit = reverse.begin(); 
                 rit < reverse.end(); rit++) {
                 // calculate the position within the vector
                 // to be used to access the other vector, if necessary for invalidation
                 position = rit - reverse.begin();
                 // get a 'threadprivate' counter
                 cit = counter[omp_get_thread_num()];
                 // reverse complement, if applicable
                 (*rit.*rvcfunc_reverse)();
                 // clip adaptors, if applicable
                 madaptor_func(*rit, ud.madaptors);
                 // clip primers, if applicable
                 primer_func(*rit);
		 // trimm, if asked for
		 trim_func(*rit);
                 // adapt gathered limits
                 rit->adapt_limits();
                 // gather stats and filter
                 (cit->*q)(*rit, &forward, position);
            }
	}
        // update and delete counter resources
        for (auto it = counter.begin(); it != counter.end(); it++) {
            (*sp) += *(*it);
            delete *it;
        }
        counter.clear();
        #ifdef TIMER
        }
        #endif
        
        // start collapsing, if requested
	if (ud.collapse) {
            collapse_func(forward, reverse);
	}

        #ifdef TIMER
        {
            Timer writeout_timer("writeout_timer");
        #endif


        #pragma omp parallel for
        for (size_t i = 0; i < 2; i++) {
            if (i == 0) {
                // write output for forward mate
                out_func(forward, out1, err1);
            }
            if (i == 1 && reverse.size() > 1) {
                // write ouput for reverse mate
                out_func(reverse, out2, err2);
            }
        }
        #ifdef TIMER
        }
        #endif

        try {
            out1.close();
        } catch (...) {};
        try {
            out2.close();
        }
        catch (...) {};
        try {
            err1.close();
        } catch (...) {};
        try {
            err2.close();
        } catch (...) {};
    }

    void triml(Read &r) {
        r.set_left(r.get_l() + ud.boundaries.get<0>());
    }
    void trimr(Read &r) {
	r.set_right(r.get_r() - ud.boundaries.get<1>());
    }
    void trimlr(Read &r) {
        r.set_left(r.get_l() + ud.boundaries.get<0>());
	r.set_right(r.get_r() - ud.boundaries.get<1>());
    }

    void collapse(std::vector<Read> &rv, std::vector<Read> &ov) {
        // second iterator for permutation testing
        std::vector<Read>::iterator sit;
        #pragma omp parallel for private(sit)
        for (size_t i = 0; i < rv.size(); i++ ) {
            Read lhs = rv[i];
            for (sit = rv.begin() + i; sit != rv.end(); sit++) {
                // are they equal and both valid?
                if (lhs == *sit && lhs.is_valid() && (*sit).is_valid()) {
                    if (lhs >= *sit) {
                        (*sit).set_validity(false, "invalidated upon collapsing");
                    } else {
                        lhs.set_validity(false, "invalidated upon collapsing");
                    }
                }

            }
        }
    }

    void collapse_paired_end(std::vector<Read> &rv,
                             std::vector<Read> &ov) {
        // second iterator for permutation testing
        std::vector<Read>::iterator sit;
        // position of the iterator, for lookup
        size_t position;
        // various Read instances
        Read lhs, rhs, rhs2;
        #pragma omp parallel for private(sit, position, lhs, rhs, rhs2)
        for (size_t i = 0; i < rv.size(); i++ ) {
            lhs = rv[i];
            for (sit = rv.begin() + i; sit != rv.end(); sit++) {
                // are they equal and both valid?
                if (lhs == *sit && lhs.is_valid() && (*sit).is_valid()) {
                    position = sit - rv.begin();
                    rhs      = ov[i];
                    rhs2     = ov[position];
                    if (rhs == rhs2 && rhs.is_valid() && rhs2.is_valid()) {
                        if (lhs >= *sit) {
                            (*sit).set_validity(false, "invalidated upon collapsing");
                        } else {
                            lhs.set_validity(false, "invalidated upon collapsing");
                        }
                        if (rhs >= rhs2) {
                            rhs2.set_validity(false, "invalidated upon collapsing");
                        } else {
                            rhs.set_validity(false, "invalidated upon collapsing");
                        }
                    }
                }
            }
        }
    }

   

}
