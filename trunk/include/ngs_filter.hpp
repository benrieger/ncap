#include <string>
#include <vector>
#include <iostream>

#include <boost/tuple/tuple.hpp>

#include "ngs_read.hpp"
#include "ngs_statistics.hpp"

#ifndef NGS_FILTER
#define NGS_FILTER

namespace ngs {
    //! a limit indicator: start to trim, overlapp length, number of errors
    typedef boost::tuple<unsigned short int, 
                         unsigned short int, 
                         unsigned short int> limits;
    
    //! simple filter (adaptor trimming, clipping, quality, etc.) -- no collapse feature
    void read_filter_parallelio(
                     void (*fadaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*madaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*primer_func) (Read &r),
                     void (*out_func)(const std::vector<Read>& r, 
                                     std::ostream &o, std::ostream &e),
                     void (*trim_func)(Read &r)
                     );

    //! simple filter (adaptor trimming, clipping, quality, etc.)
    void read_filter(void (*fadaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*madaptor_func) (Read &r, const std::vector<std::string> &a),
                     void (*primer_func) (Read &r),
                     void (*out_func)(const std::vector<Read> &r, 
                                     std::ostream &o, std::ostream &e),
                     void (*trim_func)(Read &r),
                     void (*collapse_func) (std::vector<Read> &rv,
                                            std::vector<Read> &ov)
                     );

    //! trim read with set boundaries
    void triml(Read &r);
    //! trim read with set boundaries
    void trimr(Read &r);
    //! trim read with set boundaries
    void trimlr(Read &r);
    //! collapse functionality for single end reads
    void collapse(std::vector<Read> &rv, 
                  std::vector<Read> &ov);
    //! collapse functionality for paired end reads
    void collapse_paired_end(std::vector<Read> &rv,
		             std::vector<Read> &ov);
 
} //~ngs
#endif //NGS_FILTER
