#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <functional>

#include "../include/ngs_statistics.hpp"
#include "../include/ngs_interface.hpp"
#include "../include/ngs_read.hpp"
#include "../include/ngs_global.hpp"

#include <omp.h>

// defining singleton pointer
ngs::GlobalStats* ngs::GlobalStats::inst_ = NULL;

ngs::GlobalStats::GlobalStats() {
}

ngs::Counter::Counter() : full_(),
                          trimmed_(),
                          readlength_() {}

/*
 * ngs::Counter::Counter(size_t readlength) {
    // the size is the readlength given
    // as it can be bigger, we add a grace value
    size = readlength + 50;
    for (size_t i = 0; i < size; i++) {
        for (char a = 0; a <= 255; a++) {
            full_[i][a]    = 0;
            trimmed_[i][a] = 0;
        }
    }
}
*/

namespace ngs {
    void Counter::qual_fil(Read &r,
                  std::vector<Read> *o_r,
                  size_t position) {
        const char* cq          = r.get_quality().c_str();
        size_t    size          = r.get_quality().size();

        // internal represetation of qua:>lity
        unsigned int* iquality     = new unsigned int[size]; 
        for (size_t i = 0; i < size; i++ ) {
            iquality[i] = (unsigned int)cq[i] + ud.offset;
        }

        // update the readlength
        readlength_[r.r - r.l] += 1;

        // overall quality
        for (size_t i = 0; i < size; i++) {
            full_[i][cq[i]] += iquality[i];
        }

        // trimmed quality
        for (size_t i = r.get_l(); i < size - (size - r.get_r()); i++ ) {
            trimmed_[i][cq[i]] += iquality[i];
            if (iquality[i] < ud.quality) {
                r.r = r.r > i ? i : r.r;
            }
        }
        if ((r.r - r.l) < ud.length) {
            r.set_validity(false, "remaining length is too short");
        }
        delete []iquality;
    }

    void Counter::qual_fil_paired_end(Read &r,
                  std::vector<Read> *o_r,
                  size_t position) {
        ngs::ReadStat qrstat, trstat;
	const char* cq          = r.get_quality().c_str();
	size_t    size          = r.get_quality().size();
        // internal represetation of quality
        unsigned int* iquality     = new unsigned int[size]; 
//        #pragma omp parallel for simd
        for (size_t i = 0; i < size; i++ ) {
            iquality[i] = (unsigned int)cq[i] + ud.offset;
        }
        
        // update the readlength
        readlength_[r.r - r.l] += 1;

        // overall quality
//        #pragma omp parallel for simd
        for (size_t i = 0; i < size; i++) {
            full_[i][cq[i]] += iquality[i];
        }

        // trimmed quality
        for (size_t i = r.get_l(); i < size - (size - r.get_r()); i++ ) {
            trimmed_[i][cq[i]] += iquality[i];
            if (iquality[i] < ud.quality) {
                r.r = r.r > i ? i : r.r;
            }
        }
        if ((r.r - r.l) < ud.length) {
            r.set_validity(false, "remaining length is too short");
            // take care of the mate
            std::stringstream other_comment;
            other_comment << "mate (" 
                          << r.get_id()
                          << ") got too short.";
            #pragma omp critical (quality_invalidate_other)
            {
                o_r->at(position).set_validity(false, other_comment.str());
            }
        }
        delete []iquality;
    }

    ngs::GlobalStats* GlobalStats::get_Statistics() {
        if (inst_ == NULL) {
            inst_ = new GlobalStats();
        }
        return inst_;
    }

    ngs::GlobalStats& GlobalStats::operator+=(const Counter c) {
        for (auto& mc: c._get_full()) {
            for (auto& mcc: mc.second) {
                this->full_[mc.first][mcc.first] += mcc.second;
            }
        }
        for (auto& mc: c._get_trimmed()) {
            for (auto& mcc: mc.second) {
                this->full_[mc.first][mcc.first] += mcc.second;
            }
        }
        return *this;
    }
} //~ngs
