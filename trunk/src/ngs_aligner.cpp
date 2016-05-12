#include <vector>
#include <string>
#include <iostream>

#include <omp.h>

#include "../include/ngs_read.hpp"
#include "../include/ngs_interface.hpp"
#include "../include/ngs_global.hpp"

namespace ngs {
    /* The following function was the original test dummy function.
     * It is not in use, any longer.
    void trim_adaptor(Read &r, const std::string &adaptor) {
	size_t        rlength      = r.get_length();
	size_t        offset       = rlength - ud.minimum;
        std::string::const_iterator seq_it, a_it; 
        for (size_t i = offset; i-- > 0; ) {
            unsigned short int correct = 0;
            unsigned short int errors  = 0;
            for (seq_it = r.s_begin() + i, a_it = adaptor.begin(); 
                 seq_it != r.s_end() && a_it != adaptor.end(); 
                 ++seq_it, ++a_it) {
                // REMARK: tested and not better at all times
                //if (__builtin_expect(*seq_it == *a_it, 0)) correct += 1;
                if (*seq_it == *a_it) correct += 1;
                else errors += 1;
                if (errors > ud.a_errors) break;
            }
            // if the number of errors is below our minimum, we consider it a hit
            if (errors <= ud.a_errors && (rlength - (errors + correct)) >= ud.length) {
                r.set_right(i);
                return;
            }
        }
    }
    */

    void trim_adaptors(Read &r, const std::vector<std::string> &adaptors) {
	size_t        rlength      = r.get_length();
	size_t        offset       = rlength - ud.minimum;
        std::string::const_iterator seq_it, a_it;
        std::vector<std::string>::const_iterator va_it;

        // loop over given adaptors
        for (va_it = adaptors.begin(); va_it != adaptors.end(); va_it++) {
            //loop over all possible offsets
            for (size_t i = offset; i == 0; i--) {
                unsigned short int correct = 0;
                unsigned short int errors  = 0;
                for (seq_it = r.s_begin() + i + r.get_l(), a_it = (*va_it).begin(); 
                     seq_it != r.s_begin() + r.get_r() && a_it != (*va_it).end(); 
                     ++seq_it, ++a_it) {
                     if (*seq_it == *a_it) correct += 1;
                     else errors += 1;
                     if (errors > ud.a_errors) break;
                }
                // if the number of errors is below our minimum, we consider it a hit
                if (errors <= ud.a_errors && (rlength - (errors + correct)) >= ud.length) {
                    r.set_right(i);
                    return;
                }
            }
        }
    }

    void trim_primer(Read &r) {
	size_t        rlength      = r.get_length();
        std::string::const_iterator seq_it, p_it;

        for (size_t i = 0; i < (ud.minimum > ud._primer.size() - ud.minimum ? ud.minimum : ud._primer.size()-ud.minimum);
            ++i) {
            unsigned short int correct = 0;
            unsigned short int errors  = 0;
            for (seq_it = r.s_begin(), p_it = ud._primer.begin() + i;
                 seq_it != r.s_end(), p_it != ud._primer.end();
                 ++seq_it, ++p_it ) {
                 if (*seq_it == *p_it) correct += 1;
                 else errors += 1;
                 if (errors > ud.a_errors) break;
            }
            // if the number of errors is below our minimum, we consider it a hit
            if (errors <= ud.a_errors && (rlength - (errors + correct)) >= ud.length) {
                r.set_left(p_it - ud._primer.begin());
                return;
            }
        }
    }

    void trim_primers(Read &r) {
	size_t        rlength      = r.get_length();
        std::string::const_iterator seq_it, p_it;
        std::vector<std::string>::const_iterator vp_it;

        // loop over given primers
        for (vp_it = ud.primers.begin(); vp_it != ud.primers.end(); vp_it++) {
            for (size_t i = 0; i < (ud.minimum > ud._primer.size() - ud.minimum ? ud.minimum : ud._primer.size()-ud.minimum);
                ++i) {
                unsigned short int correct = 0;
                unsigned short int errors  = 0;
                for (seq_it = r.s_begin() + r.get_l(), p_it = (*vp_it).begin() + i;
                     seq_it != r.s_begin() + r.get_r(), p_it != (*vp_it).end();
                     ++seq_it, ++p_it ) {
                    if (*seq_it == *p_it) correct += 1;
                    else errors += 1;
                    if (errors > ud.a_errors) break;
                }
                // if the number of errors is below our minimum, we consider it a hit
                if (errors <= ud.a_errors && (rlength - (errors + correct)) >= ud.length) {
                    r.set_left(p_it - ud._primer.begin());
                    return;

                }
            }
        }
    }
    
} //~ngs
